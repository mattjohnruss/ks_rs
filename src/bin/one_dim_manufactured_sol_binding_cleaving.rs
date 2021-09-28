use ks_rs::adr::one_dim::{BoundaryCondition, Cell, DomainParams, Problem1D, ProblemFunctions, Variable};
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    SspRungeKutta33,
};

use std::fs;
use std::io::{Write, BufWriter, BufReader};
use std::path::Path;
use structopt::StructOpt;
use serde::{Serialize, Deserialize};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BindingCleaving {
    pub d_c_s: f64,
    pub d_phi_i: f64,
    pub pe: f64,
    pub alpha_plus: f64,
    pub alpha_minus: f64,
    pub gamma_ui: f64,
    pub gamma_bi: f64,
}

#[allow(non_camel_case_types)]
#[repr(usize)]
#[derive(Clone, Copy, Debug)]
pub enum BindingCleavingVariable {
    C_U = 0,
    C_B = 1,
    C_S = 2,
    PHI_I = 3,
}

impl BindingCleavingVariable {
    pub const N_VARIABLE: usize = 4;
}

pub use BindingCleavingVariable::*;

impl From<Variable> for BindingCleavingVariable {
    fn from(var: Variable) -> Self {
        match var.0 as usize {
            0 => C_U,
            1 => C_B,
            2 => C_S,
            3 => PHI_I,
            _ => panic!("invalid variable number"),
        }
    }
}

impl From<BindingCleavingVariable> for Variable {
    fn from(var: BindingCleavingVariable) -> Self {
        Variable(var as usize)
    }
}

impl ProblemFunctions for BindingCleaving {
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => 1.0,
            C_B => 0.0,
            C_S => self.d_c_s,
            PHI_I => self.d_phi_i,
        }
    }

    fn velocity_p_at_midpoint(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => self.pe,
            C_B => 0.0,
            C_S => self.pe,
            PHI_I => 0.0,
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => {
                - self.alpha_plus * problem.var(C_U, cell)
                    + self.alpha_minus * problem.var(C_B, cell)
                    - self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
            }
            C_B => {
                self.alpha_plus * problem.var(C_U, cell)
                    - self.alpha_minus * problem.var(C_B, cell)
                    - self.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
            }
            C_S => {
                    self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    + self.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
            }
            PHI_I => {
                0.0
            }
        }
    }

    fn forcing(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => {
                - self.pe + self.alpha_plus + self.gamma_ui - (self.alpha_plus + self.alpha_minus + self.gamma_ui) * problem.x(cell)
            }
            C_B => {
                - self.alpha_plus + (self.alpha_plus + self.alpha_minus + self.gamma_bi) * problem.x(cell)
            }
            C_S => {
                let x = problem.x(cell);
                2.0 * self.d_c_s + self.pe * (1.0 - 2.0 * x) - self.gamma_bi * x + (x - 1.0) * self.gamma_ui
            }
            PHI_I => {
                0.0
            }
        }
    }

    fn left_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(1.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => BoundaryCondition::Flux(0.0),
        }
    }

    fn right_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(0.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => BoundaryCondition::Flux(0.0),
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

fn set_initial_conditions(problem: &mut Problem1D<BindingCleaving>) {
    for cell in problem.interior_cells() {
        let x = problem.x(cell);

        *problem.var_mut(C_U, cell) = 1.0 - x;
        *problem.var_mut(C_B, cell) = x;
        *problem.var_mut(C_S, cell) = x * (1.0 - x);
        *problem.var_mut(PHI_I, cell) = 1.0;
    }

    problem.update_ghost_cells();
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let model: BindingCleaving = {
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

    let mut problem = Problem1D::new(BindingCleavingVariable::N_VARIABLE, domain, model);
    problem.set_variable_names(&["$C_u$", "$C_b$", "$C_s$", "$\\\\phi_i$"])?;

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
