use ks_rs::adr::one_dim::{
    BoundaryCondition, Cell, DomainParams, Problem1D, ProblemFunctions, Variable,
};
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    SspRungeKutta33,
    //RungeKutta44,
};
use serde::{Deserialize, Serialize};
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
    #[structopt(long, default_value = "1.0e-6")]
    dt: f64,
    #[structopt(long, default_value = "1000")]
    outputs_per_cycle: usize,
    #[structopt(long, default_value = "res")]
    dir: String,
    #[structopt(long = "config")]
    config_path: String,
    #[structopt(long, short = "s")]
    suppress_full_output: bool,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(default)]
struct BindingFluctuatingBCs {
    d: f64,
    pe: f64,
    alpha: f64,
    k: f64,
    t_p: f64,
    n_p: usize,
}

impl Default for BindingFluctuatingBCs {
    fn default() -> Self {
        Self {
            d: 1.0,
            pe: 1.0,
            alpha: 10.0,
            k: 10.0,
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
        match var.0 as usize {
            0 => C_U,
            1 => C_B,
            _ => panic!("invalid variable number"),
        }
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

    fn velocity_p_at_midpoint(
        &self,
        _problem: &Problem1D<Self>,
        var: Variable,
        _cell: Cell,
    ) -> f64 {
        match var.into() {
            C_U => self.pe,
            C_B => 0.0,
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        let beta = self.alpha / self.k;
        match var.into() {
            C_U => -self.alpha * problem.var(C_U, cell) + beta * problem.var(C_B, cell),
            C_B => self.alpha * problem.var(C_U, cell) - beta * problem.var(C_B, cell),
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

fn trace(problem: &Problem1D<BindingFluctuatingBCs>, mut trace_writer: impl Write) -> Result<()> {
    let c_u_total = problem.integrate_solution(C_U);
    let c_b_total = problem.integrate_solution(C_B);

    if let BoundaryCondition::Dirichlet(c_u_0) = problem.functions.left_bc(problem, C_U.into()) {
        writeln!(
            &mut trace_writer,
            "{:.6e} {:.8e} {:.8e} {:.6e}",
            problem.time, c_u_total, c_b_total, c_u_0
        )?;
    } else {
        return Err("No Dirichlet condition for C_U - shouldn't happen.".into());
    }

    Ok(())
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
    let outputs_per_cycle = opt.outputs_per_cycle;
    let dir = opt.dir;
    let suppress_full_output = opt.suppress_full_output;

    let domain = DomainParams { n_cell, width: 1.0 };

    let mut problem = Problem1D::new(
        BindingFluctuatingBCsVariable::N_VARIABLE,
        domain,
        binding_pulsating_bcs,
    );
    problem.set_variable_names(&["$C_u$", "$C_b$"])?;

    set_initial_conditions(&mut problem);

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);
    //let mut rk44 = RungeKutta44::new(problem.n_dof);

    if !suppress_full_output {
        let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
        let mut buf_writer = BufWriter::new(file);
        problem.output(&mut buf_writer)?;
        buf_writer.flush()?;
    }

    let trace_file = fs::File::create(dir_path.join("trace.csv"))?;
    let mut trace_writer = BufWriter::new(trace_file);
    writeln!(&mut trace_writer, "t C_u_total C_b_total C_u_0")?;

    trace(&problem, &mut trace_writer)?;

    let mut i = 1;
    let mut outputs = 1;

    let t_max_cycles_plus_five = problem.functions.t_p * (problem.functions.n_p + 5) as f64;

    let output_time_interval = problem.functions.t_p / outputs_per_cycle as f64;

    while problem.time < t_max_cycles_plus_five.min(t_max) {
        ssp_rk33.step(&mut problem, dt);
        //rk44.step(&mut problem, dt);

        if problem.time >= outputs as f64 * output_time_interval {
            println!("Outputting at time = {}, i = {}", problem.time, i);
            if !suppress_full_output {
                let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", outputs)))?;
                let mut buf_writer = BufWriter::new(file);
                problem.output(&mut buf_writer)?;
            }
            trace(&problem, &mut trace_writer)?;
            outputs += 1;
        }
        i += 1;
    }

    Ok(())
}
