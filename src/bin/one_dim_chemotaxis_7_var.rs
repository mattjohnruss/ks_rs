use ks_rs::adr::one_dim::{DomainParams, Problem1D};
use ks_rs::models::chemotaxis_7_var::*;
use ks_rs::steady_state::SteadyStateDetector;
use ks_rs::timestepping::{ExplicitTimeStepper, SspRungeKutta33};
use ks_rs::utilities::cos_ramp;

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
    #[structopt(long)]
    ssd_threshold: f64,
}

fn set_initial_conditions(problem: &mut Problem1D<Chemotaxis>) {
    // Set relevant parameters to their homestatic values
    problem.functions.p.m = problem.functions.p.m_h;
    problem.functions.p.j_phi_i = problem.functions.p.j_phi_i_h;

    for cell in problem.interior_cells() {
        let x = problem.x(cell);

        let p = &problem.functions.p;

        let phi_i_steady = p.j_phi_i * ((p.m_h / p.d_phi_i).sqrt() * x).cosh() / ((p.m_h / p.d_phi_i).sqrt().sinh() * (p.d_phi_i * p.m_h).sqrt());

        *problem.var_mut(C_U, cell) = cos_ramp(x, 10.0);
        *problem.var_mut(C_B, cell) = 0.0;
        *problem.var_mut(C_S, cell) = 0.0;
        *problem.var_mut(PHI_I, cell) = phi_i_steady;
        *problem.var_mut(PHI_M, cell) = problem.functions.p.phi_m_init;
        *problem.var_mut(PHI_C_U, cell) = 0.0;
        *problem.var_mut(PHI_C_B, cell) = 0.0;
    }

    problem.update_ghost_cells();
}

fn trace_header(mut trace_writer: impl Write) -> Result<()> {
    writeln!(
        &mut trace_writer,
        "t C_u^{{tot}} C_b^{{tot}} C_s^{{tot}} phi_i^{{tot}} phi_m^{{tot}} phi_{{C_u}}^{{tot}} phi_{{C_b}}^{{tot}} -F_{{phi_i}}(x=1) -F_{{phi_{{C_b}}}}(x=0) m j_{{phi_i}} state"
    )?;
    Ok(())
}

fn trace(problem: &Problem1D<Chemotaxis>, state: &State, mut trace_writer: impl Write) -> Result<()> {
    let c_u_total = problem.integrate_solution(C_U);
    let c_b_total = problem.integrate_solution(C_B);
    let c_s_total = problem.integrate_solution(C_S);
    let phi_i_total = problem.integrate_solution(PHI_I);
    let phi_m_total = problem.integrate_solution(PHI_M);
    let phi_c_u_total = problem.integrate_solution(PHI_C_U);
    let phi_c_b_total = problem.integrate_solution(PHI_C_B);
    writeln!(
        &mut trace_writer,
        "{:.6e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {}",
        problem.time,
        c_u_total,
        c_b_total,
        c_s_total,
        phi_i_total,
        phi_m_total,
        phi_c_u_total,
        phi_c_b_total,
        -problem.boundary_flux_right(PHI_I.into()),
        -problem.boundary_flux_left(PHI_C_B.into()),
        problem.functions.p.m,
        problem.functions.p.j_phi_i,
        state.to_f64()
    )?;
    Ok(())
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let chemotaxis_params = {
        let config_file = fs::File::open(&opt.config_path)?;
        let reader = BufReader::new(config_file);
        serde_json::from_reader(reader)?
    };

    println!("{:#?}", opt);
    println!("{:#?}", chemotaxis_params);

    let chemotaxis = Chemotaxis::without_forcing(chemotaxis_params);

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

    let mut state = State::HomeostasisInitial;

    set_initial_conditions(&mut problem);

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let mut ssd = SteadyStateDetector::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let trace_file = fs::File::create(dir_path.join("trace.csv"))?;
    let mut trace_writer = BufWriter::new(trace_file);

    trace_header(&mut trace_writer)?;
    trace(&problem, &state, &mut trace_writer)?;

    let ssd_threshold = opt.ssd_threshold;

    let mut i = 1;
    let mut outputs = 1;

    while problem.time < t_max {
        match state {
            State::HomeostasisInitial => {
                if ssd.is_steady_state(&problem, ssd_threshold) {
                    println!("Steady state reached at t = {} (within threshold {:e})", problem.time, ssd_threshold);
                    state = State::Inflammation(problem.time);
                    // update params that are changed immediately to inflammation values - just m
                    let p = &mut problem.functions.p;
                    p.m = p.m_i_factor * p.m_h;
                }
            }
            State::Inflammation(t) => {
                // update params that involve time lag to inflammation values - just j_phi_i
                if problem.time >= t + problem.functions.p.t_j_phi_i_lag {
                    let p = &mut problem.functions.p;
                    p.j_phi_i = p.j_phi_i_i_factor * p.j_phi_i_h;
                }

                if problem.time >= t + problem.functions.p.t_inflammation {
                    state = State::HomeostasisReturn(problem.time);
                    // update params that are changed immediately back to homeostasis values - just m
                    let p = &mut problem.functions.p;
                    p.m = p.m_h;
                }
            }
            State::HomeostasisReturn(t) => {
                // update params that involve time lag back to homeostasis values - just j_phi_i
                if problem.time >= t + problem.functions.p.t_j_phi_i_lag {
                    let p = &mut problem.functions.p;
                    p.j_phi_i = p.j_phi_i_h;
                }

                if ssd.is_steady_state(&problem, ssd_threshold) {
                    println!("Steady state reached at t = {} (within threshold {:e})", problem.time, ssd_threshold);
                    break;
                }
            }
        }

        let dt = problem.calculate_dt();
        ssp_rk33.step(&mut problem, dt);

        if problem.time >= outputs as f64 * output_time_interval {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", outputs)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;
            trace(&problem, &state, &mut trace_writer)?;
            outputs += 1;
        }
        i += 1;
    }

    Ok(())
}
