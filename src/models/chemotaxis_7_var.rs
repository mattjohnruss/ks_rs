use crate::adr::one_dim::{
    Problem1D,
    ProblemFunctions,
    Variable,
    Cell,
    BoundaryCondition,
    Face,
};
use serde::{Serialize, Deserialize};

#[derive(Clone, Debug)]
pub struct Chemotaxis<F>
{
    pub p: ChemotaxisParameters,
    pub f: F,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ChemotaxisParameters {
    pub phi_bar_over_c_bar: f64,
    pub pe: f64,
    pub alpha_plus: f64,
    pub alpha_minus: f64,
    pub beta_plus: f64,
    pub beta_minus: f64,
    pub n_ccr7: f64,
    //pub k_inhib: f64,
    //pub n_inhib: f64,
    pub gamma_ui: f64,
    pub gamma_um: f64,
    pub gamma_bi: f64,
    pub gamma_bm: f64,
    pub q_u: f64,
    pub q_b: f64,
    pub q_s: f64,
    pub d_c_s: f64,
    pub d_phi_i: f64,
    pub d_phi_m: f64,
    pub d_phi_c_u: f64,
    pub d_phi_c_b: f64,
    pub chi_u: f64,
    pub chi_b: f64,
    pub chi_s: f64,
    pub r: f64,
    pub m_h: f64,
    pub m_i: f64,
    #[serde(skip)]
    pub m: f64,
    //pub p: f64,
    //pub s: f64,
    pub j_phi_c_b_bar_h: f64,
    pub j_phi_c_b_bar_i: f64,
    #[serde(skip)]
    pub j_phi_c_b_bar: f64,
    pub j_phi_i_bar_h: f64,
    pub j_phi_i_bar_i: f64,
    #[serde(skip)]
    pub j_phi_i_bar: f64,
    pub phi_i_init: f64,
    pub phi_m_init: f64,
    pub t_1: f64,
    pub t_2: f64,
}

#[allow(non_camel_case_types)]
#[repr(usize)]
#[derive(Clone, Copy, Debug)]
pub enum ChemotaxisVariable {
    C_U = 0,
    C_B = 1,
    C_S = 2,
    PHI_I = 3,
    PHI_M = 4,
    PHI_C_U = 5,
    PHI_C_B = 6,
}

impl ChemotaxisVariable {
    pub const N_VARIABLE: usize = 7;
}

pub use ChemotaxisVariable::*;

impl From<Variable> for ChemotaxisVariable {
    fn from(var: Variable) -> Self {
        match var.0 as usize {
            0 => C_U,
            1 => C_B,
            2 => C_S,
            3 => PHI_I,
            4 => PHI_M,
            5 => PHI_C_U,
            6 => PHI_C_B,
            _ => panic!("invalid variable number"),
        }
    }
}

impl From<ChemotaxisVariable> for Variable {
    fn from(var: ChemotaxisVariable) -> Self {
        Variable(var as usize)
    }
}

impl<F> ProblemFunctions for Chemotaxis<F>
where
    F: Fn(Variable, f64, f64) -> f64
{
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var.into() {
            C_U => 1.0,
            C_B => 0.0,
            C_S => self.p.d_c_s,
            PHI_I => self.p.d_phi_i,
            PHI_M => self.p.d_phi_m,
            PHI_C_U => self.p.d_phi_c_u,
            PHI_C_B => self.p.d_phi_c_b,
        }
    }

    fn velocity_p_at_midpoint(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var.into() {
            C_U => self.p.pe,
            C_B => 0.0,
            C_S => self.p.pe,
            PHI_I => 0.0,
            PHI_M => 0.0,
            PHI_C_U => 0.0,
            PHI_C_B => {
                self.p.chi_u * problem.dvar_dx_p_at_midpoint(C_U.into(), cell)
                    + self.p.chi_b * problem.dvar_dx_p_at_midpoint(C_B.into(), cell)
                    + self.p.chi_s * problem.dvar_dx_p_at_midpoint(C_S.into(), cell)
            }
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        let c_bar_over_phi_bar = 1.0 / self.p.phi_bar_over_c_bar;

        //let inhib = self.p.k_inhib.powf(self.p.n_inhib) / (self.p.k_inhib.powf(self.p.n_inhib) / problem.var(C_S, cell).powf(self.p.n_inhib));
        let inhib = 1.0;

        match var.into() {
            C_U => {
                - self.p.alpha_plus * problem.var(C_U, cell)
                    + self.p.alpha_minus * problem.var(C_B, cell)
                    - self.p.n_ccr7 * self.p.beta_plus * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    + self.p.n_ccr7 * self.p.phi_bar_over_c_bar * self.p.beta_minus * problem.var(PHI_C_U, cell)
                    - inhib * self.p.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    - inhib * self.p.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    - self.p.q_u * problem.var(PHI_I, cell) * problem.var(C_U, cell)
            }
            C_B => {
                self.p.alpha_plus * problem.var(C_U, cell)
                    - self.p.alpha_minus * problem.var(C_B, cell)
                    - self.p.n_ccr7 * self.p.beta_plus * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    + self.p.n_ccr7 * self.p.phi_bar_over_c_bar * self.p.beta_minus * problem.var(PHI_C_B, cell)
                    - inhib * self.p.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    - inhib * self.p.gamma_bm * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.p.q_b * problem.var(PHI_I, cell) * problem.var(C_B, cell)
            }
            C_S => {
                inhib * self.p.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    + inhib * self.p.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    + inhib * self.p.gamma_bi * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    + inhib * self.p.gamma_bm * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.p.q_s * problem.var(PHI_I, cell) * problem.var(C_S, cell)
            }
            PHI_I => {
                self.p.r * problem.var(PHI_I, cell) * (1.0 - problem.var(PHI_I, cell))
                    - self.p.m * problem.var(PHI_I, cell)
            },
            PHI_M => {
                self.p.m * problem.var(PHI_I, cell)
                    - c_bar_over_phi_bar * self.p.beta_plus * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    + self.p.beta_minus * problem.var(PHI_C_U, cell)
                    - c_bar_over_phi_bar * self.p.beta_plus * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    + self.p.beta_minus * problem.var(PHI_C_B, cell)
            },
            PHI_C_U => {
                - self.p.alpha_plus * problem.var(PHI_C_U, cell)
                    + self.p.alpha_minus * problem.var(PHI_C_B, cell)
                    + c_bar_over_phi_bar * self.p.beta_plus * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    - self.p.beta_minus * problem.var(PHI_C_U, cell)
            }
            PHI_C_B => {
                self.p.alpha_plus * problem.var(PHI_C_U, cell)
                    - self.p.alpha_minus * problem.var(PHI_C_B, cell)
                    + c_bar_over_phi_bar * self.p.beta_plus * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    - self.p.beta_minus * problem.var(PHI_C_B, cell)
            }
        }
    }

    fn forcing(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        (self.f)(var, problem.x(cell), problem.time)
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
                let flux = -self.p.j_phi_c_b_bar * problem.var_point_value_at_face_for_dirichlet_bcs(PHI_C_B.into(), Cell(1), Face::West);
                BoundaryCondition::Flux(flux)
            }
        }
    }

    fn right_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(0.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => BoundaryCondition::Flux(-self.p.j_phi_i_bar),
            PHI_M => BoundaryCondition::Flux(0.0),
            PHI_C_U => BoundaryCondition::Flux(0.0),
            PHI_C_B => BoundaryCondition::Flux(0.0),
        }
    }
}
