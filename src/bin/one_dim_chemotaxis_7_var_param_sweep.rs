use ks_rs::adr::one_dim::{DomainParams, Problem1D};
use ks_rs::models::chemotaxis_7_var::*;
use ks_rs::timestepping::{ExplicitTimeStepper, SspRungeKutta33};
use ks_rs::utilities::lhsu;

use ndarray::parallel::prelude::*;
use ndarray::Axis;
use std::fs;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use structopt::StructOpt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

#[derive(StructOpt, Debug)]
#[structopt(name = "chemotaxis", rename_all = "verbatim")]
struct Opt {
    #[structopt(long, default_value = "101")]
    n_cell: usize,
    #[structopt(long, default_value = "100.0")]
    t_max: f64,
    #[structopt(long, default_value = "0.01")]
    output_time_interval: f64,
    #[structopt(long, default_value = "res")]
    dir: String,
    #[structopt(long = "config")]
    config_path: String,
    #[structopt(long)]
    j_phi_i_right_i_min: f64,
    #[structopt(long)]
    j_phi_i_right_i_max: f64,
    #[structopt(long)]
    j_phi_c_b_left_i_min: f64,
    #[structopt(long)]
    j_phi_c_b_left_i_max: f64,
    #[structopt(long)]
    m_i_min: f64,
    #[structopt(long)]
    m_i_max: f64,
    #[structopt(long)]
    n_parameter_sample: usize,
}

fn set_initial_conditions(problem: &mut Problem1D<Chemotaxis>) {
    fn cos_ramp(x: f64, n: f64) -> f64 {
        use std::f64::consts::PI;
        if x < 1.0 / n {
            0.5 * (1.0 + (n * PI * x).cos())
        } else {
            0.0
        }
    }

    for cell in problem.interior_cells() {
        let x = problem.x(cell);

        *problem.var_mut(C_U, cell) = cos_ramp(x, 10.0);
        *problem.var_mut(C_B, cell) = 0.0;
        *problem.var_mut(C_S, cell) = 0.0;
        *problem.var_mut(PHI_I, cell) = problem.functions.phi_i_init;
        *problem.var_mut(PHI_M, cell) = problem.functions.phi_m_init;
        *problem.var_mut(PHI_C_U, cell) = 0.0;
        *problem.var_mut(PHI_C_B, cell) = 0.0;
    }

    problem.update_ghost_cells();
}

fn inflammation_status(problem: &Problem1D<Chemotaxis>) -> f64 {
    let t_1 = problem.functions.t_1;
    let t_2 = problem.functions.t_2;
    let time = problem.time;

    // Simplest piecewise constant inflammation status
    // TODO make this flexible/generic to allow different ramp functions?
    if time < t_1 || time > t_2 {
        0.0
    } else {
        1.0
    }
}

fn update_params(problem: &mut Problem1D<Chemotaxis>) {
    let i_s = inflammation_status(problem);
    problem.functions.m = (1.0 - i_s) * problem.functions.m_h + i_s * problem.functions.m_i;
    problem.functions.j_phi_c_b_left =
        (1.0 - i_s) * problem.functions.j_phi_c_b_left_h + i_s * problem.functions.j_phi_c_b_left_i;
    problem.functions.j_phi_i_right =
        (1.0 - i_s) * problem.functions.j_phi_i_right_h + i_s * problem.functions.j_phi_i_right_i;
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let chemotaxis = {
        let config_file = fs::File::open(&opt.config_path)?;
        let reader = BufReader::new(config_file);
        serde_json::from_reader(reader)?
    };

    println!("{:#?}", opt);
    println!("{:#?}", chemotaxis);

    let n_cell = opt.n_cell;
    let t_max = opt.t_max;
    let output_time_interval = opt.output_time_interval;
    let dir = opt.dir;

    let domain = DomainParams { n_cell, width: 1.0 };

    let mut problem = Problem1D::new(ChemotaxisVariable::N_VARIABLE, domain, chemotaxis);
    problem.set_variable_names(&[
        "$C_u$",
        "$C_b$",
        "$C_s$",
        "$\\\\phi_i$",
        "$\\\\phi_m$",
        "$\\\\phi_{C_u}$",
        "$\\\\phi_{C_b}$",
    ])?;

    set_initial_conditions(&mut problem);

    let homeostasis_path: PathBuf = [&dir, "homeostasis"].iter().collect();
    fs::create_dir_all(&homeostasis_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let file = fs::File::create(&homeostasis_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let mut i = 1;
    let mut outputs = 1;

    // timestep until just before the the parameters are switched to their inflammation values
    while problem.time < t_max && problem.time < problem.functions.t_1 {
        update_params(&mut problem);

        let dt = problem.calculate_dt();
        ssp_rk33.step(&mut problem, dt);

        if problem.time >= outputs as f64 * output_time_interval {
            let file =
                fs::File::create(homeostasis_path.join(format!("output_{:05}.csv", outputs)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;
            outputs += 1;
        }
        i += 1;
    }

    // get a latin hypercube sample of the parameter space for the unknown inflammation parameters
    let param_min = &[
        opt.j_phi_i_right_i_min,
        opt.j_phi_c_b_left_i_min,
        opt.m_i_min,
    ];
    let param_max = &[
        opt.j_phi_i_right_i_max,
        opt.j_phi_c_b_left_i_max,
        opt.m_i_max,
    ];

    let samples = lhsu(param_min, param_max, opt.n_parameter_sample);

    println!("{}", samples);

    // loop over the inflammation parameter value samples
    samples
        .axis_iter(Axis(0))
        .into_par_iter()
        .enumerate()
        .try_for_each(|(sample_idx, sample)| -> Result<()> {
            println!("{}: {}", sample_idx, sample);

            // clone the problem state - it will be independent of the inflammation parameter values
            let mut problem = problem.clone();
            let mut ssp_rk33 = ssp_rk33.clone();
            let mut i = i;
            let mut outputs = outputs;

            // set the relevant parameter values from the current sample
            problem.functions.j_phi_i_right_i = sample[0];
            problem.functions.j_phi_c_b_left_i = sample[1];
            problem.functions.m_i = sample[2];

            let inflammation_path: PathBuf = [&dir, &sample_idx.to_string()].iter().collect();
            fs::create_dir_all(&inflammation_path)?;

            // continue timestepping until t_max is reached
            while problem.time < t_max {
                update_params(&mut problem);

                let dt = problem.calculate_dt();
                ssp_rk33.step(&mut problem, dt);

                if problem.time >= outputs as f64 * output_time_interval {
                    let file = fs::File::create(
                        inflammation_path.join(format!("output_{:05}.csv", outputs)),
                    )?;
                    let mut buf_writer = BufWriter::new(file);
                    println!("Outputting at time = {}, i = {}", problem.time, i);
                    problem.output(&mut buf_writer)?;
                    outputs += 1;
                }
                i += 1;
            }
            Ok(())
        })?;

    Ok(())
}