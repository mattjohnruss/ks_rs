use crate::stencil::{first_order, second_order};
use crate::stencil::apply;
use ndarray::prelude::*;
use std::io::prelude::*;
use crate::utilities::minmod;
use std::fmt;
use std::cell::RefCell;
use std::rc::Rc;

#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
#[repr(usize)]
enum Variable {
    RhoBar = 0,
    C = 1,
}

#[derive(Debug, Clone, Copy)]
enum Face {
    East,
    West,
}

#[allow(dead_code)]
#[derive(Debug)]
pub enum KellerSegelICs {
    /// (rho_bar_init, c_init)
    Constant(f64, f64),
    /// (rho_bar_init, c_init, pert_size)
    Perturbed(f64, f64, f64),
    /// (height, width)
    Gaussian(f64, f64),
    Exact,
}

// FIXME This is ridiculous. We are only using Rc<RefCell<Box<...>>> for the functions to make
// the borrow checker happy
pub struct KellerSegelExactSolution {
    rho_bar_solution: Rc<RefCell<Box<dyn Fn(f64, f64, &KellerSegelParameters) -> f64>>>,
    c_solution: Rc<RefCell<Box<dyn Fn(f64, f64, &KellerSegelParameters) -> f64>>>,
}

impl KellerSegelExactSolution
{
    pub fn new<F1, F2>(rho_bar_solution: F1, c_solution: F2) -> Self
    where F1: Fn(f64, f64, &KellerSegelParameters) -> f64 + 'static,
          F2: Fn(f64, f64, &KellerSegelParameters) -> f64 + 'static
    {
        KellerSegelExactSolution {
            rho_bar_solution: Rc::new(RefCell::new(Box::new(rho_bar_solution))),
            c_solution: Rc::new(RefCell::new(Box::new(c_solution))),
        }
    }
}

/// Minimal impl of Debug so we can derive Debug on the other types
/// Cannot derive Debug on Fn
impl fmt::Debug for KellerSegelExactSolution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "KellerSegelExactSolutionFn {{ rho_bar_solution, c_solution }}")
    }
}

pub struct KellerSegelForces {
    rho_bar_force: Box<dyn Fn(f64, f64, &KellerSegelParameters) -> f64>,
    c_force: Box<dyn Fn(f64, f64, &KellerSegelParameters) -> f64>,
}

impl KellerSegelForces
{
    pub fn new<F1, F2>(rho_bar_force: F1, c_force: F2) -> Self
    where F1: Fn(f64, f64, &KellerSegelParameters) -> f64 + 'static,
          F2: Fn(f64, f64, &KellerSegelParameters) -> f64 + 'static
    {
        KellerSegelForces {
            rho_bar_force: Box::new(rho_bar_force),
            c_force: Box::new(c_force),
        }
    }
}

/// Minimal impl of Debug so we can derive Debug on the other types
/// Cannot derive Debug on Fn
impl fmt::Debug for KellerSegelForces {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "KellerSegelForces {{ rho_bar_force, c_force }}")
    }
}

#[derive(Debug)]
pub struct KellerSegelParameters {
    pub ics: KellerSegelICs,
    pub diffusivity: f64,
    pub r: f64,
    pub chi: f64,
    pub gamma_rho: f64,
    pub gamma_c: f64,
    pub n_interior_cell_1d: usize,
    pub length: f64,
    pub dx: f64,
    pub exact_solution: Option<KellerSegelExactSolution>,
    pub forces: Option<KellerSegelForces>,
}

impl Default for KellerSegelParameters {
    fn default() -> Self {
        let n_interior_cell_1d = 3;
        let length = 1.0;

        KellerSegelParameters {
            ics: KellerSegelICs::Constant(1.0, 1.0),
            diffusivity: 1.0,
            r: 0.0,
            chi: 1.0,
            gamma_rho: 1.0,
            gamma_c: 1.0,
            n_interior_cell_1d,
            length,
            dx: length / n_interior_cell_1d as f64,
            exact_solution: None,
            forces: None,
        }
    }
}

#[derive(Debug)]
pub struct KellerSegelProblem1D {
    pub p: KellerSegelParameters,
    data: Array<f64, ndarray::Ix2>,
    rhs_buffer: Array<f64, ndarray::Ix2>,
    pub time: f64,
}

impl KellerSegelProblem1D {
    pub fn with_params(p: KellerSegelParameters) -> Self {
        let n_interior_cell_1d = p.n_interior_cell_1d;
        let n_cell = n_interior_cell_1d + 2;
        let n_variable = 2;

        KellerSegelProblem1D {
            p,
            data: Array::zeros((n_variable, n_cell)),
            rhs_buffer: Array::zeros((2, n_interior_cell_1d)),
            time: 0.0,
        }
    }

    pub fn output<W: Write>(&self, mut buffer: W) -> std::io::Result<()> {
        match self.p.exact_solution {
            Some(_) => buffer.write_all("t x rho\\\\_bar c rho\\\\_bar\\\\_exact c\\\\_exact\n".as_bytes())?,
            None => buffer.write_all("t x rho\\\\_bar c\n".as_bytes())?,
        }

        for cell in 1..=self.p.n_interior_cell_1d {
            let x = self.x(cell);
            let rho_bar = self.u(Variable::RhoBar, cell);
            let c = self.u(Variable::C, cell);
            let time = self.time;

            match &self.p.exact_solution {
                Some(exact_solution) => {
                    let rho_bar_exact = exact_solution.rho_bar_solution.borrow()(time, x, &self.p);
                    let c_exact = exact_solution.c_solution.borrow()(time, x, &self.p);
                    buffer.write_all(format!("{} {} {} {} {} {}\n", time, x, rho_bar, c, rho_bar_exact, c_exact).as_bytes())?;
                }
                None => buffer.write_all(format!("{} {} {} {}\n", time, x, rho_bar, c).as_bytes())?,
            }
        }
        Ok(())
    }

    /// x-coordinate of the centre of cell `cell`
    fn x(&self, cell: usize) -> f64 {
        (self.index(cell) as f64 - 0.5) * self.p.dx
    }

    /// Get the value of `var` in cell `cell` at time history `t`
    fn u(&self, var: Variable, cell: usize) -> f64 {
        let idx = self.index(cell);
        self.data[(var as usize, idx)]
    }

    /// Get the value of `var` in cell `cell` at time history `t` (mut)
    fn u_mut(&mut self, var: Variable, cell: usize) -> &mut f64 {
        let idx = self.index(cell);
        &mut self.data[(var as usize, idx)]
    }

    /// Translate from cell number to index of the internal data storage
    fn index(&self, cell: usize) -> usize {
        cell
    }

    /// Set the initial conditions according to the value of `p.ics` and associated
    /// parameters
    pub fn set_initial_conditions(&mut self) {
        self.data.fill(0.0);

        match self.p.ics {
            KellerSegelICs::Constant(rho_bar_init, c_init) => {
                for cell in 1..=self.p.n_interior_cell_1d {
                    *self.u_mut(Variable::RhoBar, cell) = rho_bar_init;
                    *self.u_mut(Variable::C, cell) = c_init;
                }
            }
            KellerSegelICs::Perturbed(rho_bar_init, c_init, pert_size) => {
                use rand::distributions::Uniform;
                use rand::{thread_rng, Rng};

                let uniform = Uniform::new_inclusive(-pert_size, pert_size);

                for cell in 1..=self.p.n_interior_cell_1d {
                    // Only the chemoattractant density is perturbed
                    *self.u_mut(Variable::RhoBar, cell) = rho_bar_init;
                    *self.u_mut(Variable::C, cell) = c_init + thread_rng().sample(uniform);
                }
            }
            KellerSegelICs::Gaussian(height, width) => {
                for cell in 1..=self.p.n_interior_cell_1d {
                    let x = self.x(cell) - 0.5 * self.p.length;

                    *self.u_mut(Variable::RhoBar, cell) = height * (-width * x * x).exp();
                    *self.u_mut(Variable::C, cell) = 0.5 * height * (-0.5 * width * x * x).exp();
                }
            }
            KellerSegelICs::Exact => {
                // FIXME This is ridiculous. The problem is that we need to borrow self mutably in
                // order to assign values, but we need to (immutably) borrow self to get to the
                // exact solution functions in self.p.exact_solution...
                if self.p.exact_solution.is_some() {
                    let (rho_bar_solution, c_solution) = {
                        let KellerSegelExactSolution { rho_bar_solution, c_solution } = self.p.exact_solution.as_ref().unwrap();
                        (rho_bar_solution.clone(), c_solution.clone())
                    };

                    for cell in 1..=self.p.n_interior_cell_1d {
                        let rho_bar = rho_bar_solution.borrow()(self.time, self.x(cell), &self.p);
                        let c = c_solution.borrow()(self.time, self.x(cell), &self.p);
                        *self.u_mut(Variable::RhoBar, cell) = rho_bar;
                        *self.u_mut(Variable::C, cell) = c;
                    }
                } else {
                    panic!("ICs from exact solution requested, but an exact solution was not provided.");
                }
            }
        }
    }

    /// Equation (2.3)(a)
    fn f_p(&self, cell: usize) -> f64 {
        let u_p = self.u_p_at_midpoint(cell);
        let drho_dx_p = self.drho_dx_p_at_midpoint(cell);
        self.p.chi * self.rho_point_value_i_p(cell) * u_p - self.p.diffusivity * drho_dx_p
    }

    /// Equation (2.3)(b)
    fn f_m(&self, cell: usize) -> f64 {
        let u_m = self.u_m_at_midpoint(cell);
        let drho_dx_m = self.drho_dx_m_at_midpoint(cell);
        self.p.chi * self.rho_point_value_i_m(cell) * u_m - self.p.diffusivity * drho_dx_m
    }

    /// Equation (2.4)(a)
    fn drho_dx_p_at_midpoint(&self, cell: usize) -> f64 {
        apply(&first_order::FORWARD_1, cell, |i| {
            self.u(Variable::RhoBar, i) / self.p.dx
        })
    }

    fn drho_dx_m_at_midpoint(&self, cell: usize) -> f64 {
        self.drho_dx_p_at_midpoint(cell - 1)
    }

    /// Equation (2.4)(c)
    fn u_p_at_midpoint(&self, cell: usize) -> f64 {
        apply(&first_order::FORWARD_1, cell, |i| {
            self.u(Variable::C, i) / self.p.dx
        })
    }

    fn u_m_at_midpoint(&self, cell: usize) -> f64 {
        self.u_p_at_midpoint(cell - 1)
    }

    /// Equation (2.5)(a)
    fn rho_point_value_i_p(&self, cell: usize) -> f64 {
        match self.u_p_at_midpoint(cell) {
            u_p if u_p > 0.0 => self.rho_point_value_at_face(cell, Face::East),
            _ => self.rho_point_value_at_face(cell + 1, Face::West),
        }
    }

    fn rho_point_value_i_m(&self, cell: usize) -> f64 {
        match self.u_m_at_midpoint(cell) {
            u_m if u_m > 0.0 => self.rho_point_value_at_face(cell - 1, Face::East),
            _ => self.rho_point_value_at_face(cell, Face::West),
        }
    }

    /// Equation (2.7)
    /// Indices have been shifted compared to the paper
    fn rho_point_value_at_face(&self, cell: usize, face: Face) -> f64 {
        let rho_bar = self.u(Variable::RhoBar, cell);
        match face {
            Face::East => rho_bar + 0.5 * self.p.dx * self.drho_dx(cell),
            Face::West => rho_bar - 0.5 * self.p.dx * self.drho_dx(cell),
        }
    }

    /// Equation (2.8)(a)
    fn drho_dx(&self, cell: usize) -> f64 {
        let drho_bar_central = apply(&second_order::CENTRAL_1, cell, |i| {
            self.u(Variable::RhoBar, i)
        });

        let test_p = self.u(Variable::RhoBar, cell) + 0.5 * drho_bar_central;
        let test_m = self.u(Variable::RhoBar, cell) - 0.5 * drho_bar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            drho_bar_central / self.p.dx
        } else {
            let drho_bar_forward = apply(&first_order::FORWARD_1, cell, |i| {
                self.u(Variable::RhoBar, i)
            });

            let drho_bar_backward = apply(&first_order::BACKWARD_1, cell, |i| {
                self.u(Variable::RhoBar, i)
            });

            minmod(&[2.0 * drho_bar_forward / self.p.dx,
                     drho_bar_central / self.p.dx,
                     2.0 * drho_bar_backward / self.p.dx ])
        }
    }

    fn fc_p(&self, cell: usize) -> f64 {
        - self.u_p_at_midpoint(cell)
    }

    fn fc_m(&self, cell:usize) -> f64 {
        - self.u_m_at_midpoint(cell)
    }

    fn rhs_rho_bar(&self, cell: usize) -> f64 {
        let mut result = 0.0;

        // flux
        let f_m = if cell == 1 {
            0.0
        } else {
            self.f_m(cell)
        };

        let f_p = if cell == self.p.n_interior_cell_1d {
            0.0
        } else {
            self.f_p(cell)
        };

        result += -(f_p - f_m) / self.p.dx;

        // logistic growth
        let u = self.u(Variable::RhoBar, cell);
        result += self.p.r * u * (1.0 - u);

        // forcing
        if let Some(forces) = &self.p.forces {
            result += forces.rho_bar_force.as_ref()(self.time, self.x(cell), &self.p);
        }

        result
    }

    fn rhs_c(&self, cell: usize) -> f64 {
        let mut result = 0.0;

        // laplacian (flux or what?)
        let fc_m = if cell == 1 {
            0.0
        } else {
            self.fc_m(cell)
        };

        let fc_p = if cell == self.p.n_interior_cell_1d {
            0.0
        } else {
            self.fc_p(cell)
        };

        result += -(fc_p - fc_m) / self.p.dx;

        // reaction terms
        result += -self.p.gamma_c * self.u(Variable::C, cell);
        result += self.p.gamma_rho * self.u(Variable::RhoBar, cell);

        // forcing
        if let Some(forces) = &self.p.forces {
            result += forces.c_force.as_ref()(self.time, self.x(cell), &self.p);
        }

        result
    }

    fn update_ghost_cells(&mut self) {
        // TODO currently zeroth-order extrapolation
        *self.u_mut(Variable::RhoBar, 0) = self.u(Variable::RhoBar, 1);
        *self.u_mut(Variable::C, 0) = self.u(Variable::C, 1);

        *self.u_mut(Variable::RhoBar, self.p.n_interior_cell_1d + 1) = self.u(Variable::RhoBar, self.p.n_interior_cell_1d);
        *self.u_mut(Variable::C, self.p.n_interior_cell_1d + 1) = self.u(Variable::C, self.p.n_interior_cell_1d);
    }

    fn step_euler_forward_helper(&mut self, dt: f64) {
        self.update_ghost_cells();

        for cell in 1..=self.p.n_interior_cell_1d {
            let idx = self.index(cell) - 1;
            self.rhs_buffer[(Variable::RhoBar as usize, idx)] = self.rhs_rho_bar(cell);
            self.rhs_buffer[(Variable::C as usize, idx)] = self.rhs_c(cell);
        }

        {
            let mut slice = self.data.slice_mut(s![.., 1..=self.p.n_interior_cell_1d]);
            self.rhs_buffer *= dt;
            slice += &self.rhs_buffer;
        }
    }

    pub fn step_euler_forward(&mut self, dt: f64) {
        self.time += dt;
        self.step_euler_forward_helper(dt);
    }

    pub fn step_ssp_rk3(&mut self, dt: f64) {
        // FIXME This is probably quite inefficient with all the copies etc
        self.time += dt;

        // save current solution
        let u_old = self.data.slice(s![.., 1..=self.p.n_interior_cell_1d]).to_owned();

        // an Euler step gives w_1
        self.step_euler_forward_helper(dt);

        // another Euler step gives the brackets in the w_2 eqn
        self.step_euler_forward_helper(dt);

        // calulate w_2
        let w_2 = {
            let temp = self.data.slice(s![.., 1..=self.p.n_interior_cell_1d]);
            0.75 * &u_old + 0.25 * &temp
        };

        // assign w_2 back into the main storage
        let mut slice = self.data.slice_mut(s![.., 1..=self.p.n_interior_cell_1d]);
        slice.assign(&w_2);

        // another Euler step gives the brackets in the w_3 eqn
        self.step_euler_forward_helper(dt);

        // calulate w_{n+1}
        let w_np1 = {
            let temp = self.data.slice(s![.., 1..=self.p.n_interior_cell_1d]);
            (1.0 / 3.0) * &u_old + (2.0 / 3.0) * &temp
        };

        // assign w_{n+1} back into the main storage
        let mut slice = self.data.slice_mut(s![.., 1..=self.p.n_interior_cell_1d]);
        slice.assign(&w_np1);
    }
}
