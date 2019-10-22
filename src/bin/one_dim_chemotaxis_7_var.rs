use ks_rs::adr::one_dim::{
    DomainParams,
    Problem1D,
    ProblemFunctions,
    Variable,
    Cell,
    Face,
    BoundaryCondition,
};
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
struct Chemotaxis {
    phi_i_max_over_c_0: f64,
    pe: f64,
    alpha_1: f64,
    beta_1: f64,
    alpha_2: f64,
    beta_2: f64,
    alpha_3: f64,
    beta_3: f64,
    alpha_4: f64,
    beta_4: f64,
    k_inhib: f64,
    n_inhib: f64,
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
    d_phi_c_u: f64,
    d_phi_c_b: f64,
    nu_u: f64,
    nu_b: f64,
    nu_s: f64,
    chi_u: f64,
    chi_b: f64,
    chi_s: f64,
    r: f64,
    m: f64,
    p: f64,
    s: f64,
    j_phi_c_b_left: f64,
    phi_i_init: f64,
}

impl Default for Chemotaxis {
    fn default() -> Self {
        Chemotaxis {
            phi_i_max_over_c_0: 0.01,
            pe: 1.0,
            alpha_1: 10.0,
            beta_1: 5.0,
            alpha_2: 10.0,
            beta_2: 5.0,
            alpha_3: 10.0,
            beta_3: 5.0,
            alpha_4: 10.0,
            beta_4: 5.0,
            k_inhib: 1.0,
            n_inhib: 0.0,
            gamma_ui: 0.0,
            gamma_um: 0.0,
            gamma_bi: 0.0,
            gamma_bm: 0.0,
            q_u: 0.0,
            q_b: 0.0,
            q_s: 0.0,
            d_su: 100.0,
            d_iu: 0.01,
            d_mu: 0.01,
            d_phi_c_u: 0.01,
            d_phi_c_b: 0.01,
            nu_u: 0.0,
            nu_b: 1.0,
            nu_s: 0.0,
            chi_u: 0.0,
            chi_b: 1.0,
            chi_s: 0.0,
            r: 10.0,
            m: 5.0,
            p: 10.0,
            s: 0.0,
            j_phi_c_b_left: 1.0,
            phi_i_init: 0.1,
        }
    }
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
    PHI_C_U = 5,
    PHI_C_B = 6,
}

impl ChemotaxisVariable {
    const N_VARIABLE: usize = 7;
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

impl ProblemFunctions for Chemotaxis {
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => 1.0,
            C_B => 0.0,
            C_S => self.d_su,
            PHI_I => self.d_iu,
            PHI_M => self.d_mu,
            PHI_C_U => self.d_phi_c_u,
            PHI_C_B => self.d_phi_c_b,
        }
    }

    fn velocity_p_at_midpoint(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => self.pe,
            C_B => 0.0,
            C_S => self.pe,
            PHI_I => 0.0,
            PHI_M => 0.0,
            PHI_C_U => 0.0,
            PHI_C_B => {
                self.nu_u * self.chi_u * problem.dvar_dx_p_at_midpoint(C_U.into(), cell)
                    + self.nu_b * self.chi_b * problem.dvar_dx_p_at_midpoint(C_B.into(), cell)
                    + self.nu_s * self.chi_s * problem.dvar_dx_p_at_midpoint(C_S.into(), cell)
            }
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        let c_0_over_phi_i_max = 1.0 / self.phi_i_max_over_c_0;

        let inhib = self.k_inhib.powf(self.n_inhib) / (self.k_inhib.powf(self.n_inhib) / problem.var(C_S, cell).powf(self.n_inhib));

        match var.into() {
            C_U => {
                - self.alpha_1 * problem.var(C_U, cell)
                    + self.beta_1 * problem.var(C_B, cell)
                    - self.alpha_2 * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    + self.phi_i_max_over_c_0 * self.beta_2 * problem.var(PHI_C_U, cell)
                    - inhib * self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    - inhib * self.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    - self.q_u * problem.var(PHI_I, cell) * problem.var(C_U, cell)
            }
            C_B => {
                self.alpha_1 * problem.var(C_U, cell)
                    - self.beta_1 * problem.var(C_B, cell)
                    - self.alpha_4 * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    + self.phi_i_max_over_c_0 * self.beta_4 * problem.var(PHI_C_B, cell)
                    - inhib * self.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    - inhib * self.gamma_bm * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.q_b * problem.var(PHI_I, cell) * problem.var(C_B, cell)
            }
            C_S => {
                inhib * self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    + inhib * self.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    + inhib * self.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    + inhib * self.gamma_bm * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.q_s * problem.var(PHI_I, cell) * problem.var(C_S, cell)
            }
            PHI_I => {
                self.r * problem.var(PHI_I, cell) * (1.0 - problem.var(PHI_I, cell))
                    - self.m * problem.var(PHI_I, cell)
            },
            PHI_M => {
                self.m * problem.var(PHI_I, cell)
                    - c_0_over_phi_i_max * self.alpha_2 * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    + self.beta_2 * problem.var(PHI_C_U, cell)
                    - c_0_over_phi_i_max * self.alpha_4 * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    + self.beta_4 * problem.var(PHI_C_B, cell)
            },
            PHI_C_U => {
                c_0_over_phi_i_max * self.alpha_2 * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    - self.beta_2 * problem.var(PHI_C_U, cell)
                    - self.alpha_3 * problem.var(PHI_C_U, cell)
                    + self.beta_3 * problem.var(PHI_C_B, cell)
            }
            PHI_C_B => {
                self.alpha_3 * problem.var(PHI_C_U, cell)
                    - self.beta_3 * problem.var(PHI_C_B, cell)
                    + c_0_over_phi_i_max * self.alpha_4 * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    - self.beta_4 * problem.var(PHI_C_B, cell)
            }
        }
    }

    fn left_bc(&self, problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(1.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => BoundaryCondition::Flux(0.0),
            PHI_M => BoundaryCondition::Flux(0.0),
            PHI_C_U => BoundaryCondition::Flux(0.0),
            PHI_C_B => {
                let flux = -self.j_phi_c_b_left * problem.var_point_value_at_face_for_dirichlet_bcs(PHI_C_B.into(), Cell(1), Face::West);
                BoundaryCondition::Flux(flux)
            }
        }
    }

    fn right_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(0.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => BoundaryCondition::Flux(0.0),
            PHI_M => BoundaryCondition::Flux(0.0),
            PHI_C_U => BoundaryCondition::Flux(0.0),
            PHI_C_B => BoundaryCondition::Flux(0.0),
        }
    }
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

        *problem.var_mut(C_U, cell) = cos_ramp(x, 5.0);
        *problem.var_mut(C_B, cell) = 0.0;
        *problem.var_mut(C_S, cell) = 0.0;
        *problem.var_mut(PHI_I, cell) = problem.functions.phi_i_init;
        *problem.var_mut(PHI_M, cell) = 0.0;
        *problem.var_mut(PHI_C_U, cell) = 0.0;
        *problem.var_mut(PHI_C_B, cell) = 0.0;
    }

    problem.update_ghost_cells();
}

#[allow(dead_code)]
fn update_params(problem: &mut Problem1D<Chemotaxis>) {
    problem.functions.m = if problem.time < 3.5 {
        2.0
    } else {
        7.0
    };
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let chemotaxis: Chemotaxis = {
        let config_file = fs::File::open(&opt.config_path)?;
        let reader = BufReader::new(config_file);
        serde_json::from_reader(reader)?
    };

    let n_cell = opt.n_cell;
    let t_max = opt.t_max;
    let dt = opt.dt;
    let output_interval = opt.output_interval;
    let dir = opt.dir;

    let domain = DomainParams { n_cell, width: 1.0 };

    let mut problem = Problem1D::new(ChemotaxisVariable::N_VARIABLE, domain, chemotaxis);
    problem.set_variable_names(&["C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b"])?;

    set_initial_conditions(&mut problem);

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let mut i = 1;

    while problem.time < t_max {
        //update_params(&mut problem);
        ssp_rk33.step(&mut problem, dt);

        if i % output_interval == 0 {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;
        }
        i += 1;
    }

    Ok(())
}
