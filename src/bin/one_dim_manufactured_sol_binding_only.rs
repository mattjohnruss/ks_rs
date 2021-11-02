use ks_rs::adr::one_dim::{
    BoundaryCondition, Cell, DomainParams, Problem1D, ProblemFunctions, Variable,
};
use ks_rs::timestepping::{ExplicitTimeStepper, SspRungeKutta33};

use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use structopt::StructOpt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BindingOnly {
    pub pe: f64,
    pub alpha_plus: f64,
    pub alpha_minus: f64,
}

#[allow(non_camel_case_types)]
#[repr(usize)]
#[derive(Clone, Copy, Debug)]
pub enum BindingOnlyVariable {
    C_U = 0,
    C_B = 1,
}

impl BindingOnlyVariable {
    pub const N_VARIABLE: usize = 2;
}

pub use BindingOnlyVariable::*;

impl From<Variable> for BindingOnlyVariable {
    fn from(var: Variable) -> Self {
        match var.0 as usize {
            0 => C_U,
            1 => C_B,
            _ => panic!("invalid variable number"),
        }
    }
}

impl From<BindingOnlyVariable> for Variable {
    fn from(var: BindingOnlyVariable) -> Self {
        Variable(var as usize)
    }
}

impl ProblemFunctions for BindingOnly {
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => 1.0,
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
        match var.into() {
            C_U => {
                -self.alpha_plus * problem.var(C_U, cell)
                    + self.alpha_minus * problem.var(C_B, cell)
            }
            C_B => {
                self.alpha_plus * problem.var(C_U, cell) - self.alpha_minus * problem.var(C_B, cell)
            }
        }
    }

    fn forcing(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => {
                -self.pe + self.alpha_plus - (self.alpha_plus + self.alpha_minus) * problem.x(cell)
            }
            C_B => -self.alpha_plus + (self.alpha_plus + self.alpha_minus) * problem.x(cell),
        }
    }

    fn left_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(1.0),
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

fn set_initial_conditions(problem: &mut Problem1D<BindingOnly>) {
    for cell in problem.interior_cells() {
        let x = problem.x(cell);

        *problem.var_mut(C_U, cell) = 1.0 - x;
        *problem.var_mut(C_B, cell) = x;
    }

    problem.update_ghost_cells();
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let model: BindingOnly = {
        let config_file = fs::File::open(&opt.config_path)?;
        let reader = BufReader::new(config_file);
        serde_json::from_reader(reader)?
    };

    println!("{:#?}", opt);
    println!("{:#?}", model);

    let n_cell = opt.n_cell;
    let t_max = opt.t_max;
    let output_time_interval = opt.output_time_interval;
    let dir = opt.dir;

    let domain = DomainParams { n_cell, width: 1.0 };

    let mut problem = Problem1D::new(BindingOnlyVariable::N_VARIABLE, domain, model);
    problem.set_variable_names(&["$C_u$", "$C_b$"])?;

    set_initial_conditions(&mut problem);

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let mut i = 1;
    let mut outputs = 1;

    while problem.time < t_max {
        let dt = problem.calculate_dt();
        ssp_rk33.step(&mut problem, dt);

        if problem.time >= outputs as f64 * output_time_interval {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", outputs)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;
            outputs += 1;
        }
        i += 1;
    }

    Ok(())
}
