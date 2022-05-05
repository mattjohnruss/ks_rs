use ks_rs::adr::one_dim::{DomainParams, Problem1D, Variable};
use ks_rs::models::chemotaxis_7_var::*;
use ks_rs::timestepping::{ExplicitTimeStepper, SspRungeKutta33};

use std::f64::consts::PI;
use std::fs;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use structopt::StructOpt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

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
}

fn set_initial_conditions(problem: &mut Problem1D<Chemotaxis>) {
    // Set maturation and cell influx rates to homeostasis values as in this test problem we're not
    // simulating inflammation
    let p = &mut problem.functions.p;
    p.m = p.m_h;
    p.j_phi_i = p.j_phi_i_h;

    for cell in problem.interior_cells() {
        let x = problem.x(cell);

        let p = &problem.functions.p;
        // Ignore phi_i_init and set phi_i to it's exact solution
        let phi_i = 0.5 * x * x * (p.j_phi_i / p.d_phi_i);

        *problem.var_mut(C_U, cell) = 1.0 - x;
        *problem.var_mut(C_B, cell) = x;
        *problem.var_mut(C_S, cell) = x * (1.0 - x);
        *problem.var_mut(PHI_I, cell) = phi_i;
        *problem.var_mut(PHI_M, cell) = 1.0;
        *problem.var_mut(PHI_C_U, cell) = 1.0;
        *problem.var_mut(PHI_C_B, cell) = 0.5 * (1.0 - (2.0 * PI * x).cos());
    }

    problem.update_ghost_cells();
}

fn forcing(var: Variable, x: f64, _t: f64, p: &ChemotaxisParameters) -> f64 {
    match var.into() {
        C_U => {
            -p.pe
                + p.alpha_plus
                + p.n_ccr7 * p.beta_plus
                - 0.5 * (p.j_phi_i / p.d_phi_i) * x * x * (x - 1.0) * (p.q_u + p.gamma_ui)
                + p.gamma_um
                - x * (p.alpha_minus + p.alpha_plus + p.n_ccr7 * p.beta_plus + p.gamma_um)
                - p.phi_bar_over_c_bar * p.n_ccr7 * p.beta_minus
        }
        C_B => {
            -p.alpha_plus
                + 0.5 * (p.j_phi_i / p.d_phi_i) * x.powi(3) * (p.q_b + p.gamma_bi)
                + x * (p.alpha_minus + p.alpha_plus + p.n_ccr7 * p.beta_plus + p.gamma_bm)
                - p.phi_bar_over_c_bar
                    * p.n_ccr7
                    * p.beta_minus
                    * (PI * x).sin().powi(2)
        }
        C_S => {
            2.0 * p.d_c_s
                + p.pe * (1.0 - 2.0 * x)
                - x * p.gamma_bm
                + 0.5
                    * (p.j_phi_i / p.d_phi_i)
                    * x
                    * x
                    * (p.q_s * x * (x - 1.0) - p.gamma_bi * x + p.gamma_ui * (x - 1.0))
                + p.gamma_um * (x - 1.0)
        }
        PHI_I => -p.j_phi_i + 0.5 * (p.j_phi_i / p.d_phi_i) * x * x * p.m,
        PHI_M => {
            -0.5 * (p.j_phi_i / p.d_phi_i) * p.m * x * x
                - 1.5 * p.beta_minus
                + p.phi_bar_over_c_bar.recip() * p.beta_plus
                + 0.5 * p.beta_minus * (2.0 * PI * x).cos()
        }
        PHI_C_U => {
            p.alpha_plus
                + p.beta_minus
                + p.phi_bar_over_c_bar.recip() * p.beta_plus * (x - 1.0)
                - p.alpha_minus * (PI * x).sin().powi(2)
        }
        PHI_C_B => {
            -p.alpha_plus
                + (p.alpha_minus + p.beta_minus) * (PI * x).sin().powi(2)
                - p.phi_bar_over_c_bar.recip() * p.beta_plus * x
                + PI * (p.chi_b * (2.0 * PI * x).sin() - 2.0 * p.d_phi_c_b * PI * (2.0 * PI * x).cos())
        }
    }
}

fn exact_solution(var: Variable, x: f64, _time: f64, p: &ChemotaxisParameters) -> f64 {
    match var.into() {
        C_U => 1.0 - x,
        C_B => x,
        C_S => x * (1.0 - x),
        PHI_I => 0.5 * x.powi(2) * (p.j_phi_i / p.d_phi_i),
        PHI_M => 1.0,
        PHI_C_U => 1.0,
        PHI_C_B => 0.5 * (1.0 - (2.0 * PI * x).cos()),
    }
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let chemotaxis_params: ChemotaxisParameters = {
        let config_file = fs::File::open(&opt.config_path)?;
        let reader = BufReader::new(config_file);
        serde_json::from_reader(reader)?
    };

    println!("{:#?}", opt);
    println!("{:#?}", chemotaxis_params);

    let chemotaxis = Chemotaxis::with_forcing(chemotaxis_params, forcing);

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

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    // output ICs
    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    // output exact ICs
    let file = fs::File::create(dir_path.join(format!("output_exact_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output_fn(&mut buf_writer, |var, x, time| {
        exact_solution(var, x, time, &problem.functions.p)
    })?;
    buf_writer.flush()?;

    let mut i = 1;
    let mut outputs = 1;

    let start = std::time::Instant::now();

    while problem.time < t_max {
        let dt = problem.calculate_dt();
        ssp_rk33.step(&mut problem, dt);

        if problem.time >= outputs as f64 * output_time_interval {
            // output solution
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", outputs)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;

            // output exact solution
            let file = fs::File::create(dir_path.join(format!("output_exact_{:05}.csv", outputs)))?;
            let mut buf_writer = BufWriter::new(file);
            problem.output_fn(&mut buf_writer, |var, x, time| {
                exact_solution(var, x, time, &problem.functions.p)
            })?;
            buf_writer.flush()?;

            outputs += 1;
        }
        i += 1;
    }

    let duration = std::time::Instant::now() - start;
    println!("Time for {i} timesteps: {} ({:.2} timesteps / s)", duration.as_secs_f64(), i as f64 / duration.as_secs_f64());

    Ok(())
}
