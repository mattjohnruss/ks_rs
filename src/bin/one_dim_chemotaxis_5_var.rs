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
    SspRungeKutta33,
};
use std::fs;
use std::io::{Write, BufWriter};
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

struct Chemotaxis {
    pe: f64,
    alpha: f64,
    beta: f64,
    gamma_ui: f64,
    gamma_um: f64,
    gamma_bi: f64,
    gamma_bm: f64,
    q_u: f64,
    q_b: f64,
    q_s: f64,
    d_su: f64,
    d_iu: f64,
    d_mu: f64,
    nu_u: f64,
    nu_b: f64,
    nu_s: f64,
    chi_u: f64,
    chi_b: f64,
    chi_s: f64,
    r: f64,
    m: f64,
    phi_i_init: f64,
}

#[allow(non_camel_case_types)]
#[repr(usize)]
#[derive(Clone, Copy, Debug)]
enum ChemotaxisVariable {
    C_U = 0,
    C_B = 1,
    C_S = 2,
    PHI_I = 3,
    PHI_M = 4,
}

impl ChemotaxisVariable {
    const N_VARIABLE: usize = 5;
}

use ChemotaxisVariable::*;

impl From<Variable> for ChemotaxisVariable {
    fn from(var: Variable) -> Self {
        unsafe { std::mem::transmute(var.0) }
    }
}

impl From<ChemotaxisVariable> for Variable {
    fn from(var: ChemotaxisVariable) -> Self {
        Variable(var as usize)
    }
}

impl Default for Chemotaxis {
    fn default() -> Self {
        Chemotaxis {
            pe: 1.0,
            alpha: 10.0,
            beta: 5.0,
            gamma_ui: 1.0,
            gamma_um: 1.0,
            gamma_bi: 1.0,
            gamma_bm: 1.0,
            q_u: 1.0,
            q_b: 1.0,
            q_s: 1.0,
            d_su: 10.0,
            d_iu: 0.1,
            d_mu: 0.1,
            nu_u: 0.0,
            nu_b: 10.0,
            nu_s: 0.0,
            chi_u: 1.0,
            chi_b: 1.0,
            chi_s: 1.0,
            r: 10.0,
            m: 5.0,
            phi_i_init: 0.1,
        }
    }
}

impl ProblemFunctions for Chemotaxis {
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => 1.0,
            C_B => 0.0,
            C_S => self.d_su,
            PHI_I => self.d_iu,
            PHI_M => self.d_mu,
        }
    }

    fn velocity_p_at_midpoint(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => self.pe,
            C_B => 0.0,
            C_S => self.pe,
            PHI_I => 0.0,
            PHI_M => {
                self.nu_u * self.chi_u * problem.dvar_dx_p_at_midpoint(C_U.into(), cell)
                    + self.nu_b * self.chi_b * problem.dvar_dx_p_at_midpoint(C_B.into(), cell)
                    + self.nu_s * self.chi_s * problem.dvar_dx_p_at_midpoint(C_S.into(), cell)
            }
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => {
                - self.alpha * problem.var(C_U, cell)
                    + self.beta * problem.var(C_B, cell)
                    - self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    - self.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    - self.q_u * problem.var(PHI_I, cell) * problem.var(C_U, cell)
            }
            C_B => {
                self.alpha * problem.var(C_U, cell)
                    - self.beta * problem.var(C_B, cell)
                    - self.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    - self.gamma_bm * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.q_b * problem.var(PHI_I, cell) * problem.var(C_B, cell)
            }
            C_S => {
                self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    + self.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    + self.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    + self.gamma_bm * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.q_s * problem.var(PHI_I, cell) * problem.var(C_S, cell)
            }
            PHI_I => {
                self.r * problem.var(PHI_I, cell) * (1.0 - problem.var(PHI_I, cell))
                    - self.m * problem.var(PHI_I, cell)
            },
            PHI_M => {
                self.m * problem.var(PHI_I, cell)
            },
        }
    }

    fn left_bc(&self, problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        let flux = match var.into() {
            C_U => 1.0,
            C_B => 0.0,
            C_S => -problem.var(C_S, Cell(1)),
            PHI_I => 0.0,
            PHI_M => -problem.var(PHI_M, Cell(1)),
        };
        BoundaryCondition::Flux(flux)
    }

    fn right_bc(&self, problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        let flux = match var.into() {
            C_U => problem.var(C_U, Cell(problem.domain.n_cell)),
            C_B => 0.0,
            C_S => problem.var(C_S, Cell(problem.domain.n_cell)),
            PHI_I => 0.0,
            PHI_M => 0.0,
        };
        BoundaryCondition::Flux(flux)
    }
}

fn set_initial_conditions(problem: &mut Problem1D<Chemotaxis>)
{
    for cell in problem.interior_cells() {
        //let x = problem.x(cell);

        //*problem.var_mut(C_U, cell) = 1.0 - x;
        *problem.var_mut(C_U, cell) = 0.0;
        *problem.var_mut(C_B, cell) = 0.0;
        *problem.var_mut(C_S, cell) = 0.0;
        *problem.var_mut(PHI_I, cell) = problem.functions.phi_i_init;
        *problem.var_mut(PHI_M, cell) = 0.0;
    }

    problem.update_ghost_cells();
}

fn main() -> Result<()> {
    let domain = DomainParams { n_cell: 101, width: 1.0 };

    let chemotaxis = Chemotaxis::default();

    let mut problem = Problem1D::new(ChemotaxisVariable::N_VARIABLE, domain, chemotaxis);
    problem.set_variable_names(&["C_u", "C_b", "C_s", "phi_i", "phi_m"])?;

    set_initial_conditions(&mut problem);

    let dir = String::from("res_test");
    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let output_interval = 1000;
    let mut i = 1;
    let dt = 1e-6;
    let t_max = 100.0;

    while problem.time < t_max {
        ssp_rk33.step(&mut problem, dt);

        if i % output_interval == 0 {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let mut buf_writer = BufWriter::new(file);
            problem.output(&mut buf_writer)?;
        }
        i += 1;
    }

    Ok(())
}
