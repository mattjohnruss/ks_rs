use ks_rs::adr::one_dim::{
    DomainParams,
    Problem1D,
    ProblemFunctions,
    Variable,
    Cell,
    BoundaryCondition,
};
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    //SspRungeKutta33,
    RungeKutta44,
};
use std::fs;
use std::io::{Write, BufWriter, BufReader};
use std::path::Path;
use structopt::StructOpt;
use serde::{Serialize, Deserialize};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(StructOpt, Debug)]
#[structopt(name = "chemotaxis", rename_all = "verbatim")]
struct Opt {
    #[structopt(long, default_value = "101")]
    n_cell: usize,
    #[structopt(long, default_value = "100.0")]
    t_max: f64,
    #[structopt(long, default_value = "1.0e-6")]
    dt: f64,
    #[structopt(long, default_value = "1000")]
    output_interval: usize,
    #[structopt(long, default_value = "res")]
    dir: String,
    #[structopt(long = "config")]
    config_path: String,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(default)]
struct BindingFluctuatingBCs {
    d: f64,
    pe: f64,
    alpha_1: f64,
    beta_1: f64,
    t_p: f64,
    n_p: usize,
}

impl Default for BindingFluctuatingBCs {
    fn default() -> Self {
        Self {
            d: 1.0,
            pe: 1.0,
            alpha_1: 10.0,
            beta_1: 1.0,
            t_p: 1.0,
            n_p: 3,
        }
    }
}

#[allow(non_camel_case_types)]
#[repr(usize)]
#[derive(Clone, Copy, Debug)]
enum BindingFluctuatingBCsVariable {
    C_U = 0,
    C_B = 1,
}

impl BindingFluctuatingBCsVariable {
    const N_VARIABLE: usize = 2;
}

use BindingFluctuatingBCsVariable::*;

impl From<Variable> for BindingFluctuatingBCsVariable {
    fn from(var: Variable) -> Self {
        unsafe { std::mem::transmute(var.0) }
    }
}

impl From<BindingFluctuatingBCsVariable> for Variable {
    fn from(var: BindingFluctuatingBCsVariable) -> Self {
        Variable(var as usize)
    }
}

impl ProblemFunctions for BindingFluctuatingBCs {
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => self.d,
            C_B => 0.0,
        }
    }

    fn velocity_p_at_midpoint(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => self.pe,
            C_B => 0.0,
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => - self.alpha_1 * problem.var(C_U, cell) + self.beta_1 * problem.var(C_B, cell),
            C_B => self.alpha_1 * problem.var(C_U, cell) - self.beta_1 * problem.var(C_B, cell),
        }
    }

    fn left_bc(&self, problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        fn bcs(t: f64, t_p: f64, n_p: usize) -> f64 {
            if t < t_p * n_p as f64 {
                use std::f64::consts::PI;
                0.5 * (1.0 - (2.0 * PI * t / t_p).sin().signum())
            } else {
                0.0
            }
        }

        match var.into() {
            C_U => BoundaryCondition::Dirichlet(bcs(problem.time, self.t_p, self.n_p)),
            C_B => BoundaryCondition::Flux(0.0),
        }
    }

    fn right_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(0.0),
            C_B => BoundaryCondition::Flux(0.0),
        }
    }
}

fn set_initial_conditions(problem: &mut Problem1D<BindingFluctuatingBCs>) {
    for cell in problem.interior_cells() {
        *problem.var_mut(C_U, cell) = 0.0;
        *problem.var_mut(C_B, cell) = 0.0;
    }

    problem.update_ghost_cells();
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let binding_pulsating_bcs: BindingFluctuatingBCs = {
        let config_file = fs::File::open(&opt.config_path)?;
        let reader = BufReader::new(config_file);
        serde_json::from_reader(reader)?
    };

    println!("{:#?}", opt);
    println!("{:#?}", binding_pulsating_bcs);

    let n_cell = opt.n_cell;
    let t_max = opt.t_max;
    let dt = opt.dt;
    let output_interval = opt.output_interval;
    let dir = opt.dir;

    let domain = DomainParams { n_cell, width: 1.0 };

    let mut problem = Problem1D::new(BindingFluctuatingBCsVariable::N_VARIABLE, domain, binding_pulsating_bcs);
    problem.set_variable_names(&["$C_u$", "$C_b$"])?;

    set_initial_conditions(&mut problem);

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    //let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);
    let mut rk44 = RungeKutta44::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let trace_file = fs::File::create(dir_path.join("trace.csv"))?;
    let mut trace_writer = BufWriter::new(trace_file);
    writeln!(&mut trace_writer, "t C_u_total C_b_total")?;

    let mut i = 1;

    while problem.time < t_max {
        //ssp_rk33.step(&mut problem, dt);
        rk44.step(&mut problem, dt);

        if i % output_interval == 0 {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;

            let c_u_total = problem.integrate_solution(C_U);
            let c_b_total = problem.integrate_solution(C_B);

            if let BoundaryCondition::Dirichlet(c_u_0) = problem.functions.left_bc(&problem, C_U.into()) {
                writeln!(&mut trace_writer, "{} {} {} {}", problem.time, c_u_total, c_b_total, c_u_0)?;
            } else {
                eprintln!("No Dirichlet condition for C_U - shouldn't happen.");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    Ok(())
}
