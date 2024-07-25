use crate::adr::one_dim::{BoundaryCondition, Cell, Face, Problem1D, ProblemFunctions, Variable};
use serde::{Deserialize, Serialize};

type ForcingFn = fn(Variable, f64, f64, &ChemotaxisParameters) -> f64;

#[derive(Clone)]
pub struct Chemotaxis {
    pub p: ChemotaxisParameters,
    pub f: ForcingFn,
    pub state: State,
}

impl Chemotaxis {
    pub fn without_forcing(p: ChemotaxisParameters) -> Self {
        Self { p, f: |_var, _x, _t, _p| 0.0, state: State::HomeostasisInitial }
    }

    pub fn with_forcing(p: ChemotaxisParameters, f: ForcingFn) -> Self {
        Self { p, f, state: State::HomeostasisInitial }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ChemotaxisParameters {
    pub phi_bar_over_c_bar: f64,
    pub phi_bar_over_phi_max: f64,
    pub c_bar_over_e: f64,
    pub pe: f64,
    pub alpha_plus: f64,
    pub alpha_minus: f64,
    pub beta_plus: f64,
    pub beta_minus: f64,
    pub n_ccr7: f64,
    pub n: f64,
    pub a_bar: f64,
    pub gamma: f64,
    pub q_u: f64,
    pub q_b: f64,
    pub q_s: f64,
    pub d_f: f64,
    pub d_t: f64,
    pub d_c_s: f64,
    pub d_phi_i: f64,
    pub d_phi_m: f64,
    pub d_phi_c_u: f64,
    pub d_phi_c_b: f64,
    pub d_phi_c_s: f64,
    pub d_j: f64,
    pub chi_b: f64,
    pub chi_s: f64,
    pub mu_m: f64,
    pub m_h: f64,
    pub m_i_factor: f64,
    #[serde(skip)]
    pub m: f64,
    pub j_phi_i_h: f64,
    pub j_phi_i_i_factor: f64,
    #[serde(skip)]
    pub j_phi_i: f64,
    pub phi_i_init: f64,
    pub phi_m_init: f64,
    pub t_inflammation: f64,
    pub t_j_phi_i_lag: f64,
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
    PHI_C_S = 7,
    J = 8,
}

impl ChemotaxisVariable {
    pub const N_VARIABLE: usize = 9;
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
            7 => PHI_C_S,
            8 => J,
            _ => panic!("invalid variable number"),
        }
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
            C_S => self.p.d_c_s,
            PHI_I => self.p.d_phi_i,
            PHI_M => self.p.d_phi_m,
            PHI_C_U => self.p.d_phi_c_u,
            PHI_C_B => self.p.d_phi_c_b,
            PHI_C_S => self.p.d_phi_c_s,
            J => self.p.d_j,
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
            PHI_C_B => self.p.chi_b * problem.dvar_dx_p_at_midpoint(C_B.into(), cell),
            PHI_C_S => self.p.chi_s * problem.dvar_dx_p_at_midpoint(C_S.into(), cell),
            J => self.p.pe,
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        let c_bar_over_phi_bar = 1.0 / self.p.phi_bar_over_c_bar;

        let inhib = 1.0 / (1.0 + problem.var(J, cell).powf(self.p.n));

        let ecm_occupancy = self.p.c_bar_over_e
            * (problem.var(C_B, cell)
                + self.p.n_ccr7 * self.p.phi_bar_over_c_bar * problem.var(PHI_C_B, cell));

        match var.into() {
            C_U => {
                - self.p.alpha_plus * (1.0 - ecm_occupancy) * problem.var(C_U, cell)
                    + self.p.alpha_minus * problem.var(C_B, cell)
                    - self.p.n_ccr7
                        * self.p.beta_plus
                        * problem.var(C_U, cell)
                        * problem.var(PHI_M, cell)
                    + self.p.n_ccr7
                        * self.p.phi_bar_over_c_bar
                        * self.p.beta_minus
                        * problem.var(PHI_C_U, cell)
                    - inhib * self.p.gamma * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    - inhib * self.p.gamma * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                    - self.p.q_u * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                    - self.p.d_f * problem.var(C_U, cell)
            }
            C_B => {
                self.p.alpha_plus * (1.0 - ecm_occupancy) * problem.var(C_U, cell)
                    - self.p.alpha_minus * problem.var(C_B, cell)
                    - self.p.n_ccr7
                        * self.p.beta_plus
                        * problem.var(C_B, cell)
                        * problem.var(PHI_M, cell)
                    + self.p.n_ccr7
                        * self.p.phi_bar_over_c_bar
                        * self.p.beta_minus
                        * problem.var(PHI_C_B, cell)
                    - inhib * self.p.gamma * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    - inhib * self.p.gamma * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                    - self.p.q_b * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                    - self.p.d_f * problem.var(C_B, cell)
            }
            C_S => {
                - self.p.n_ccr7
                    * self.p.beta_plus
                    * problem.var(C_S, cell)
                    * problem.var(PHI_M, cell)
                + self.p.n_ccr7
                    * self.p.phi_bar_over_c_bar
                    * self.p.beta_minus
                    * problem.var(PHI_C_S, cell)
                + inhib * self.p.gamma * problem.var(PHI_I, cell) * problem.var(C_U, cell)
                + inhib * self.p.gamma * problem.var(PHI_M, cell) * problem.var(C_U, cell)
                + inhib * self.p.gamma * problem.var(PHI_I, cell) * problem.var(C_B, cell)
                + inhib * self.p.gamma * problem.var(PHI_M, cell) * problem.var(C_B, cell)
                - self.p.q_s * problem.var(PHI_I, cell) * problem.var(C_S, cell)
                - self.p.d_t * problem.var(C_S, cell)
            }
            PHI_I => {
                - self.p.m * problem.var(PHI_I, cell)
            }
            PHI_M => {
                self.p.m * problem.var(PHI_I, cell)
                    + self.p.mu_m * problem.var(PHI_M, cell)
                    - c_bar_over_phi_bar
                        * self.p.beta_plus
                        * problem.var(C_U, cell)
                        * problem.var(PHI_M, cell)
                    + self.p.beta_minus * problem.var(PHI_C_U, cell)
                    - c_bar_over_phi_bar
                        * self.p.beta_plus
                        * problem.var(C_B, cell)
                        * problem.var(PHI_M, cell)
                    + self.p.beta_minus * problem.var(PHI_C_B, cell)
                    - c_bar_over_phi_bar
                        * self.p.beta_plus
                        * problem.var(C_S, cell)
                        * problem.var(PHI_M, cell)
                    + self.p.beta_minus * problem.var(PHI_C_S, cell)
            }
            PHI_C_U => {
                - self.p.alpha_plus * (1.0 - ecm_occupancy) * problem.var(PHI_C_U, cell)
                    + self.p.alpha_minus * problem.var(PHI_C_B, cell)
                    + c_bar_over_phi_bar
                        * self.p.beta_plus
                        * problem.var(C_U, cell)
                        * problem.var(PHI_M, cell)
                    - self.p.beta_minus * problem.var(PHI_C_U, cell)
            }
            PHI_C_B => {
                self.p.alpha_plus * (1.0 - ecm_occupancy) * problem.var(PHI_C_U, cell)
                    - self.p.alpha_minus * problem.var(PHI_C_B, cell)
                    + c_bar_over_phi_bar
                        * self.p.beta_plus
                        * problem.var(C_B, cell)
                        * problem.var(PHI_M, cell)
                    - self.p.beta_minus * problem.var(PHI_C_B, cell)
            }
            PHI_C_S => {
                c_bar_over_phi_bar
                    * self.p.beta_plus
                    * problem.var(C_S, cell)
                    * problem.var(PHI_M, cell)
                    - self.p.beta_minus * problem.var(PHI_C_S, cell)
            }
            J => self.p.a_bar * problem.var(PHI_M, cell),
        }
    }

    fn forcing(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        (self.f)(var, problem.x(cell), problem.time, &self.p)
    }

    fn left_bc(&self, _problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(1.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => BoundaryCondition::Flux(0.0),
            PHI_M => BoundaryCondition::Flux(0.0),
            PHI_C_U => BoundaryCondition::Dirichlet(0.0),
            PHI_C_B => BoundaryCondition::Dirichlet(0.0),
            PHI_C_S => BoundaryCondition::Dirichlet(0.0),
            J => BoundaryCondition::Dirichlet(0.0),
        }
    }

    fn right_bc(&self, problem: &Problem1D<Self>, var: Variable) -> BoundaryCondition {
        match var.into() {
            C_U => BoundaryCondition::Dirichlet(0.0),
            C_B => BoundaryCondition::Flux(0.0),
            C_S => BoundaryCondition::Dirichlet(0.0),
            PHI_I => {
                let cell = Cell(problem.domain.n_cell);

                let phi_total = problem.var_point_value_at_face_for_dirichlet_bcs(PHI_I.into(), cell, Face::East) +
                    problem.var_point_value_at_face_for_dirichlet_bcs(PHI_M.into(), cell, Face::East) +
                    problem.var_point_value_at_face_for_dirichlet_bcs(PHI_C_U.into(), cell, Face::East) +
                    problem.var_point_value_at_face_for_dirichlet_bcs(PHI_C_B.into(), cell, Face::East) +
                    problem.var_point_value_at_face_for_dirichlet_bcs(PHI_C_S.into(), cell, Face::East);

                // factor representing the occupancy at the boundary
                let f = 1.0 - self.p.phi_bar_over_phi_max * phi_total;
                BoundaryCondition::Flux(-f * self.p.j_phi_i)
            }
            PHI_M => BoundaryCondition::Flux(0.0),
            PHI_C_U => BoundaryCondition::Flux(0.0),
            PHI_C_B => BoundaryCondition::Flux(0.0),
            PHI_C_S => BoundaryCondition::Flux(0.0),
            J => BoundaryCondition::Dirichlet(0.0),
        }
    }
}

/// Indicates whether the system is currently in homeostasis or inflammation. The `f64` value in
/// each variant is the time that state was entered.
#[derive(Clone, Debug)]
pub enum State {
    HomeostasisInitial,
    Inflammation(f64),
    HomeostasisReturn(f64),
}

impl State {
    pub fn to_f64(&self) -> f64 {
        match self {
            State::HomeostasisInitial | State::HomeostasisReturn(_) => 0.0,
            State::Inflammation(_) => 1.0,
        }
    }
}
