use crate::adr::one_dim::{
    Problem1D,
    ProblemFunctions,
    Variable,
    Cell,
    BoundaryCondition,
    Face,
};
use serde::{Serialize, Deserialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct Chemotaxis {
    pub phi_bar_over_c_0: f64,
    pub pe: f64,
    pub alpha_plus: f64,
    pub alpha_minus: f64,
    pub beta_plus: f64,
    pub beta_minus: f64,
    pub n_ccr7: f64,
    pub k_inhib: f64,
    pub n_inhib: f64,
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
    pub nu_u: f64,
    pub nu_b: f64,
    pub nu_s: f64,
    pub chi_u: f64,
    pub chi_b: f64,
    pub chi_s: f64,
    pub r: f64,
    pub m_h: f64,
    pub m_i: f64,
    #[serde(skip)]
    pub m: f64,
    pub t_1: f64,
    pub t_2: f64,
    //pub p: f64,
    //pub s: f64,
    pub j_phi_c_b_left_h: f64,
    pub j_phi_c_b_left_i: f64,
    #[serde(skip)]
    pub j_phi_c_b_left: f64,
    pub j_phi_i_right_h: f64,
    pub j_phi_i_right_i: f64,
    #[serde(skip)]
    pub j_phi_i_right: f64,
    pub phi_i_init: f64,
    pub phi_m_init: f64,
}

impl Default for Chemotaxis {
    fn default() -> Self {
        Chemotaxis {
            phi_bar_over_c_0: 0.01,
            pe: 1.0,
            alpha_plus: 10.0,
            alpha_minus: 5.0,
            beta_plus: 10.0,
            beta_minus: 5.0,
            n_ccr7: 30000.0,
            k_inhib: 1.0,
            n_inhib: 0.0,
            gamma_ui: 0.0,
            gamma_um: 0.0,
            gamma_bi: 0.0,
            gamma_bm: 0.0,
            q_u: 0.0,
            q_b: 0.0,
            q_s: 0.0,
            d_c_s: 100.0,
            d_phi_i: 0.01,
            d_phi_m: 0.01,
            d_phi_c_u: 0.01,
            d_phi_c_b: 0.01,
            nu_u: 0.0,
            nu_b: 1.0,
            nu_s: 0.0,
            chi_u: 0.0,
            chi_b: 1.0,
            chi_s: 0.0,
            r: 0.0,
            m_h: 2.0,
            m_i: 7.0,
            m: 5.0,
            t_1: 1.0,
            t_2: 2.0,
            //p: 10.0,
            //s: 0.0,
            j_phi_c_b_left_h: 1.0,
            j_phi_c_b_left_i: 1.0,
            j_phi_c_b_left: 1.0,
            j_phi_i_right_h: 0.0,
            j_phi_i_right_i: 0.0,
            j_phi_i_right: 0.0,
            phi_i_init: 0.1,
            phi_m_init: 0.1,
        }
    }
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
            C_S => self.d_c_s,
            PHI_I => self.d_phi_i,
            PHI_M => self.d_phi_m,
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
        let c_0_over_phi_bar = 1.0 / self.phi_bar_over_c_0;

        //let inhib = self.k_inhib.powf(self.n_inhib) / (self.k_inhib.powf(self.n_inhib) / problem.var(C_S, cell).powf(self.n_inhib));
        let inhib = 1.0;

        match var.into() {
            C_U => {
                - self.alpha_plus * problem.var(C_U, cell)
                    + self.alpha_minus * problem.var(C_B, cell)
                    - self.n_ccr7 * self.beta_plus * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    + self.n_ccr7 * self.phi_bar_over_c_0 * self.beta_minus * problem.var(PHI_C_U, cell)
                    - inhib * self.gamma_ui * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    - inhib * self.gamma_um * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    - self.q_u * problem.var(PHI_I, cell) * problem.var(C_U, cell)
            }
            C_B => {
                self.alpha_plus * problem.var(C_U, cell)
                    - self.alpha_minus * problem.var(C_B, cell)
                    - self.n_ccr7 * self.beta_plus * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    + self.n_ccr7 * self.phi_bar_over_c_0 * self.beta_minus * problem.var(PHI_C_B, cell)
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
                    - c_0_over_phi_bar * self.beta_plus * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    + self.beta_minus * problem.var(PHI_C_U, cell)
                    - c_0_over_phi_bar * self.beta_plus * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    + self.beta_minus * problem.var(PHI_C_B, cell)
            },
            PHI_C_U => {
                - self.alpha_plus * problem.var(PHI_C_U, cell)
                    + self.alpha_minus * problem.var(PHI_C_B, cell)
                    + c_0_over_phi_bar * self.beta_plus * problem.var(C_U, cell) * problem.var(PHI_M, cell)
                    - self.beta_minus * problem.var(PHI_C_U, cell)
            }
            PHI_C_B => {
                self.alpha_plus * problem.var(PHI_C_U, cell)
                    - self.alpha_minus * problem.var(PHI_C_B, cell)
                    + c_0_over_phi_bar * self.beta_plus * problem.var(C_B, cell) * problem.var(PHI_M, cell)
                    - self.beta_minus * problem.var(PHI_C_B, cell)
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
            PHI_I => BoundaryCondition::Flux(-self.j_phi_i_right),
            PHI_M => BoundaryCondition::Flux(0.0),
            PHI_C_U => BoundaryCondition::Flux(0.0),
            PHI_C_B => BoundaryCondition::Flux(0.0),
        }
    }
}
