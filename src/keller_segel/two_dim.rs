use crate::stencil::{Stencil, first_order, second_order};
use crate::utilities::minmod;
use ndarray::prelude::*;
use std::fmt;
use std::io::prelude::*;

#[derive(Debug, Clone, Copy)]
#[repr(usize)]
enum Variable {
    RhoBar = 0,
    C = 1,
}

#[derive(Debug, Clone, Copy)]
enum Face {
    North,
    East,
    South,
    West,
}

#[derive(Debug)]
pub enum ICs {
    /// (rho_bar_init, c_init)
    Constant(f64, f64),
    /// (rho_bar_init, c_init, pert_size)
    Perturbed(f64, f64, f64),
    /// (height, width)
    Gaussian(f64, f64),
    Exact,
}

// NOTE The closures are boxed because a closure cannot accept itself as a parameter
pub struct ExactSolution {
    rho_bar_solution: Box<dyn Fn(f64, f64, f64, &Parameters) -> f64>,
    c_solution: Box<dyn Fn(f64, f64, f64, &Parameters) -> f64>,
}

impl ExactSolution {
    pub fn new<F1, F2>(rho_bar_solution: F1, c_solution: F2) -> Self
    where
        F1: Fn(f64, f64, f64, &Parameters) -> f64 + 'static,
        F2: Fn(f64, f64, f64, &Parameters) -> f64 + 'static,
    {
        ExactSolution {
            rho_bar_solution: Box::new(rho_bar_solution),
            c_solution: Box::new(c_solution),
        }
    }
}

/// Minimal impl of Debug so we can derive Debug on the other types
/// Cannot derive Debug on Fn
impl fmt::Debug for ExactSolution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ExactSolution {{ rho_bar_solution, c_solution }}"
        )
    }
}

pub struct Forces {
    rho_bar_force: Box<dyn Fn(f64, f64, f64, &Parameters) -> f64>,
    c_force: Box<dyn Fn(f64, f64, f64, &Parameters) -> f64>,
}

impl Forces {
    pub fn new<F1, F2>(rho_bar_force: F1, c_force: F2) -> Self
    where
        F1: Fn(f64, f64, f64, &Parameters) -> f64 + 'static,
        F2: Fn(f64, f64, f64, &Parameters) -> f64 + 'static,
    {
        Forces {
            rho_bar_force: Box::new(rho_bar_force),
            c_force: Box::new(c_force),
        }
    }
}

/// Minimal impl of Debug so we can derive Debug on the other types
/// Cannot derive Debug on Fn
impl fmt::Debug for Forces {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Forces {{ rho_bar_force, c_force }}")
    }
}

#[derive(Debug)]
pub struct Parameters {
    pub ics: ICs,
    pub diffusivity: f64,
    pub r: f64,
    pub chi: f64,
    pub gamma_rho: f64,
    pub gamma_c: f64,
    pub n_interior_cell_1d: usize,
    pub length: f64,
    pub dx: f64,
    pub dy: f64,
    pub exact_solution: Option<ExactSolution>,
    pub forces: Option<Forces>,
}

impl Default for Parameters {
    fn default() -> Self {
        let n_interior_cell_1d = 3;
        let length = 1.0;

        Parameters {
            ics: ICs::Constant(1.0, 1.0),
            diffusivity: 1.0,
            r: 0.0,
            chi: 1.0,
            gamma_rho: 1.0,
            gamma_c: 1.0,
            n_interior_cell_1d,
            length,
            dx: length / n_interior_cell_1d as f64,
            dy: length / n_interior_cell_1d as f64,
            exact_solution: None,
            forces: None,
        }
    }
}

#[derive(Debug)]
pub struct Problem2D {
    pub p: Parameters,
    data: Array<f64, ndarray::Ix3>,
    rhs_buffer: Array<f64, ndarray::Ix3>,
    pub time: f64,
}

impl Problem2D {
    pub fn with_params(p: Parameters) -> Self {
        let n_interior_cell_1d = p.n_interior_cell_1d;
        let n_cell_1d = n_interior_cell_1d + 2;
        let n_variable = 2;

        Problem2D {
            p,
            data: Array::zeros((n_variable, n_cell_1d, n_cell_1d)),
            rhs_buffer: Array::zeros((n_variable, n_interior_cell_1d, n_interior_cell_1d)),
            time: 0.0,
        }
    }

    pub fn output<W: Write>(&self, mut buffer: W) -> std::io::Result<()> {
        match self.p.exact_solution {
            Some(_) => buffer
                .write_all("t x y rho\\\\_bar c rho\\\\_bar\\\\_exact c\\\\_exact\n".as_bytes())?,
            None => buffer.write_all("t x y rho\\\\_bar c\n".as_bytes())?,
        }

        for cell_x in 1..=self.p.n_interior_cell_1d {
            for cell_y in 1..=self.p.n_interior_cell_1d {
                let x = self.x(cell_x);
                let y = self.y(cell_y);

                let rho_bar = self.u(Variable::RhoBar, cell_x, cell_y);
                let c = self.u(Variable::C, cell_x, cell_y);
                let time = self.time;

                match &self.p.exact_solution {
                    Some(exact_solution) => {
                        let rho_bar_exact = (exact_solution.rho_bar_solution)(time, x, y, &self.p);
                        let c_exact = (exact_solution.c_solution)(time, x, y, &self.p);
                        buffer.write_all(
                            format!(
                                "{:.6e} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e}\n",
                                time, x, y, rho_bar, c, rho_bar_exact, c_exact
                            )
                            .as_bytes(),
                            )?;
                    }
                    None => {
                        buffer.write_all(format!("{:.6e} {:.6e} {:.6e} {:.6e} {:.6e}\n", time, x, y, rho_bar, c).as_bytes())?
                    }
                }
            }
            buffer.write_all(b"\n")?;
        }
        Ok(())
    }

    /// x-coordinate of the centre of cell `cell_x`
    fn x(&self, cell_x: usize) -> f64 {
        (cell_x as f64 - 0.5) * self.p.dx
    }

    /// y-coordinate of the centre of cell `cell_y`
    fn y(&self, cell_y: usize) -> f64 {
        (cell_y as f64 - 0.5) * self.p.dy
    }

    /// Get the value of `var` in cell `(cell_x, cell_y)`
    fn u(&self, var: Variable, cell_x: usize, cell_y: usize) -> f64 {
        let (idx_x, idx_y) = self.index(cell_x, cell_y);
        self.data[(var as usize, idx_x, idx_y)]
    }

    /// Get the value of `var` in cell `(cell_x, cell_y)` (mut)
    fn u_mut(&mut self, var: Variable, cell_x: usize, cell_y: usize) -> &mut f64 {
        let (idx_x, idx_y) = self.index(cell_x, cell_y);
        &mut self.data[(var as usize, idx_x, idx_y)]
    }

    /// Translate from cell number to index of the internal data storage
    fn index(&self, cell_x: usize, cell_y: usize) -> (usize, usize) {
        (cell_x, cell_y)
    }

    /// Set the initial conditions according to the value of `p.ics` and associated
    /// parameters
    pub fn set_initial_conditions(&mut self) {
        self.data.fill(0.0);

        match self.p.ics {
            ICs::Constant(rho_bar_init, c_init) => {
                for cell_x in 1..=self.p.n_interior_cell_1d {
                    for cell_y in 1..=self.p.n_interior_cell_1d {
                        *self.u_mut(Variable::RhoBar, cell_x, cell_y) = rho_bar_init;
                        *self.u_mut(Variable::C, cell_x, cell_y) = c_init;
                    }
                }
            }
            ICs::Perturbed(rho_bar_init, c_init, pert_size) => {
                use rand::distributions::Uniform;
                use rand::{thread_rng, Rng};

                let uniform = Uniform::new_inclusive(-pert_size, pert_size);

                for cell_x in 1..=self.p.n_interior_cell_1d {
                    for cell_y in 1..=self.p.n_interior_cell_1d {
                        // Only the chemoattractant density is perturbed
                        *self.u_mut(Variable::RhoBar, cell_x, cell_y) = rho_bar_init;
                        *self.u_mut(Variable::C, cell_x, cell_y) = c_init + thread_rng().sample(uniform);
                    }
                }
            }
            ICs::Gaussian(height, width) => {
                for cell_x in 1..=self.p.n_interior_cell_1d {
                    for cell_y in 1..=self.p.n_interior_cell_1d {
                        let x = self.x(cell_x) - 0.5 * self.p.length;
                        let y = self.y(cell_y) - 0.5 * self.p.length;

                        *self.u_mut(Variable::RhoBar, cell_x, cell_y) = height * (-width * (x * x + y * y)).exp();
                        *self.u_mut(Variable::C, cell_x, cell_y) = 0.5 * height * (-0.5 * width * (x * x + y * y)).exp();
                    }
                }
            }
            ICs::Exact => {
                if self.p.exact_solution.is_some() {
                    for cell_x in 1..=self.p.n_interior_cell_1d {
                        for cell_y in 1..=self.p.n_interior_cell_1d {
                            // NOTE: We can't move the unwrap() outside the loop as it (immutably)
                            // borrows self, which disallows the mutable borrow inside the loop. If
                            // both are inside the loop, then (presumably due to NLL?) the compiler
                            // knows the immutable borrow isn't used again and allows the mutable
                            // borrow.
                            let x = self.x(cell_x);
                            let y = self.y(cell_y);

                            let rho_bar = (self.p.exact_solution.as_ref().unwrap().rho_bar_solution)(self.time, x, y, &self.p);
                            let c = (self.p.exact_solution.as_ref().unwrap().c_solution)(self.time, x, y, &self.p);

                            *self.u_mut(Variable::RhoBar, cell_x, cell_y) = rho_bar;
                            *self.u_mut(Variable::C, cell_x, cell_y) = c;
                        }
                    }
                } else {
                    panic!("ICs from exact solution requested, but an exact solution was not provided.");
                }
            }
        }
    }

    /// Equation (2.3)(a)
    fn flux_rho_x_p(&self, cell_x: usize, cell_y: usize) -> f64 {
        let u_p = self.u_p_at_midpoint(cell_x, cell_y);
        let drho_dx_p = self.drho_dx_p_at_midpoint(cell_x, cell_y);
        self.p.chi * self.rho_point_value_x_p(cell_x, cell_y) * u_p - self.p.diffusivity * drho_dx_p
    }

    /// Equation (2.3)(a)
    fn flux_rho_x_m(&self, cell_x: usize, cell_y: usize) -> f64 {
        let u_m = self.u_m_at_midpoint(cell_x, cell_y);
        let drho_dx_m = self.drho_dx_m_at_midpoint(cell_x, cell_y);
        self.p.chi * self.rho_point_value_x_m(cell_x, cell_y) * u_m - self.p.diffusivity * drho_dx_m
    }

    /// Equation (2.3)(b)
    fn flux_rho_y_p(&self, cell_x: usize, cell_y: usize) -> f64 {
        let v_p = self.v_p_at_midpoint(cell_x, cell_y);
        let drho_dy_p = self.drho_dy_p_at_midpoint(cell_x, cell_y);
        self.p.chi * self.rho_point_value_y_p(cell_x, cell_y) * v_p - self.p.diffusivity * drho_dy_p
    }

    /// Equation (2.3)(b)
    fn flux_rho_y_m(&self, cell_x: usize, cell_y: usize) -> f64 {
        let v_m = self.v_m_at_midpoint(cell_x, cell_y);
        let drho_dy_m = self.drho_dy_m_at_midpoint(cell_x, cell_y);
        self.p.chi * self.rho_point_value_y_m(cell_x, cell_y) * v_m - self.p.diffusivity * drho_dy_m
    }

    fn flux_c_x_p(&self, cell_x: usize, cell_y: usize) -> f64 {
        -self.u_p_at_midpoint(cell_x, cell_y)
    }

    fn flux_c_x_m(&self, cell_x: usize, cell_y: usize) -> f64 {
        -self.u_m_at_midpoint(cell_x, cell_y)
    }

    fn flux_c_y_p(&self, cell_x: usize, cell_y: usize) -> f64 {
        -self.v_p_at_midpoint(cell_x, cell_y)
    }

    fn flux_c_y_m(&self, cell_x: usize, cell_y: usize) -> f64 {
        -self.v_m_at_midpoint(cell_x, cell_y)
    }

    /// Equation (2.4)(a)
    fn drho_dx_p_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        first_order::Forward1::apply(cell_x, |i| {
            self.u(Variable::RhoBar, i, cell_y) / self.p.dx
        })
    }

    /// Equation (2.4)(a)
    fn drho_dx_m_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        self.drho_dx_p_at_midpoint(cell_x - 1, cell_y)
    }

    /// Equation (2.4)(b)
    fn drho_dy_p_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        first_order::Forward1::apply(cell_y, |i| {
            self.u(Variable::RhoBar, cell_x, i) / self.p.dy
        })
    }

    /// Equation (2.4)(b)
    fn drho_dy_m_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        self.drho_dy_p_at_midpoint(cell_x, cell_y - 1)
    }

    /// Equation (2.4)(c)
    fn u_p_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        first_order::Forward1::apply(cell_x, |i| {
            self.u(Variable::C, i, cell_y) / self.p.dx
        })
    }

    /// Equation (2.4)(c)
    fn u_m_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        self.u_p_at_midpoint(cell_x - 1, cell_y)
    }

    /// Equation (2.4)(d)
    fn v_p_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        first_order::Forward1::apply(cell_y, |i| {
            self.u(Variable::C, cell_x, i) / self.p.dy
        })
    }

    /// Equation (2.4)(d)
    fn v_m_at_midpoint(&self, cell_x: usize, cell_y: usize) -> f64 {
        self.v_p_at_midpoint(cell_x, cell_y - 1)
    }

    /// Equation (2.5)(a)
    fn rho_point_value_x_p(&self, cell_x: usize, cell_y: usize) -> f64 {
        match self.u_p_at_midpoint(cell_x, cell_y) {
            u_p if u_p > 0.0 => self.rho_point_value_at_face(cell_x, cell_y, Face::East),
            _ => self.rho_point_value_at_face(cell_x + 1, cell_y, Face::West),
        }
    }

    /// Equation (2.5)(a)
    fn rho_point_value_x_m(&self, cell_x: usize, cell_y: usize) -> f64 {
        match self.u_m_at_midpoint(cell_x, cell_y) {
            u_m if u_m > 0.0 => self.rho_point_value_at_face(cell_x - 1, cell_y, Face::East),
            _ => self.rho_point_value_at_face(cell_x, cell_y, Face::West),
        }
    }

    /// Equation (2.5)(b)
    fn rho_point_value_y_p(&self, cell_x: usize, cell_y: usize) -> f64 {
        match self.v_p_at_midpoint(cell_x, cell_y) {
            v_p if v_p > 0.0 => self.rho_point_value_at_face(cell_x, cell_y, Face::North),
            _ => self.rho_point_value_at_face(cell_x, cell_y + 1, Face::South),
        }
    }

    /// Equation (2.5)(b)
    fn rho_point_value_y_m(&self, cell_x: usize, cell_y: usize) -> f64 {
        match self.u_m_at_midpoint(cell_x, cell_y) {
            u_m if u_m > 0.0 => self.rho_point_value_at_face(cell_x, cell_y - 1, Face::North),
            _ => self.rho_point_value_at_face(cell_x, cell_y, Face::South),
        }
    }

    /// Equation (2.7)
    /// Indices have been shifted compared to the paper
    fn rho_point_value_at_face(&self, cell_x: usize, cell_y: usize, face: Face) -> f64 {
        let rho_bar = self.u(Variable::RhoBar, cell_x, cell_y);
        match face {
            Face::North => rho_bar + 0.5 * self.p.dy * self.drho_dy(cell_x, cell_y),
            Face::East => rho_bar + 0.5 * self.p.dx * self.drho_dx(cell_x, cell_y),
            Face::South => rho_bar - 0.5 * self.p.dy * self.drho_dy(cell_x, cell_y),
            Face::West => rho_bar - 0.5 * self.p.dx * self.drho_dx(cell_x, cell_y),
        }
    }

    /// Equation (2.8)(a)
    fn drho_dx(&self, cell_x: usize, cell_y: usize) -> f64 {
        let drho_bar_central = second_order::Central1::apply(cell_x, |i| {
            self.u(Variable::RhoBar, i, cell_y)
        });

        let test_p = self.u(Variable::RhoBar, cell_x, cell_y) + 0.5 * drho_bar_central;
        let test_m = self.u(Variable::RhoBar, cell_x, cell_y) - 0.5 * drho_bar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            drho_bar_central / self.p.dx
        } else {
            let drho_bar_forward = first_order::Forward1::apply(cell_x, |i| {
                self.u(Variable::RhoBar, i, cell_y)
            });

            let drho_bar_backward = first_order::Backward1::apply(cell_x, |i| {
                self.u(Variable::RhoBar, i, cell_y)
            });

            minmod(&[
                2.0 * drho_bar_forward / self.p.dx,
                drho_bar_central / self.p.dx,
                2.0 * drho_bar_backward / self.p.dx,
            ])
        }
    }

    /// Equation (2.8)(b)
    fn drho_dy(&self, cell_x: usize, cell_y: usize) -> f64 {
        let drho_bar_central = second_order::Central1::apply(cell_y, |i| {
            self.u(Variable::RhoBar, cell_x, i)
        });

        let test_p = self.u(Variable::RhoBar, cell_x, cell_y) + 0.5 * drho_bar_central;
        let test_m = self.u(Variable::RhoBar, cell_x, cell_y) - 0.5 * drho_bar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            drho_bar_central / self.p.dy
        } else {
            let drho_bar_forward = first_order::Forward1::apply(cell_y, |i| {
                self.u(Variable::RhoBar, cell_x, i)
            });

            let drho_bar_backward = first_order::Backward1::apply(cell_y, |i| {
                self.u(Variable::RhoBar, cell_x, i)
            });

            minmod(&[
                2.0 * drho_bar_forward / self.p.dy,
                drho_bar_central / self.p.dy,
                2.0 * drho_bar_backward / self.p.dy,
            ])
        }
    }

    fn rhs_rho_bar(&self, cell_x: usize, cell_y: usize) -> f64 {
        let mut result = 0.0;

        // x-flux
        result += {
            let flux_rho_x_m = if cell_x == 1 {
                0.0
            } else {
                self.flux_rho_x_m(cell_x, cell_y)
            };
            let flux_rho_x_p = if cell_x == self.p.n_interior_cell_1d {
                0.0
            } else {
                self.flux_rho_x_p(cell_x, cell_y)
            };

            -(flux_rho_x_p - flux_rho_x_m) / self.p.dx
        };

        // y-flux
        result += {
            let flux_rho_y_m = if cell_y == 1 {
                0.0
            } else {
                self.flux_rho_y_m(cell_x, cell_y)
            };
            let flux_rho_y_p = if cell_y == self.p.n_interior_cell_1d {
                0.0
            } else {
                self.flux_rho_y_p(cell_x, cell_y)
            };

            -(flux_rho_y_p - flux_rho_y_m) / self.p.dy
        };

        // logistic growth
        result += {
            let u = self.u(Variable::RhoBar, cell_x, cell_y);
            self.p.r * u * (1.0 - u)
        };

        // forcing
        if let Some(forces) = &self.p.forces {
            let x = self.x(cell_x);
            let y = self.y(cell_y);
            result += forces.rho_bar_force.as_ref()(self.time, x, y, &self.p);
        }

        result
    }

    fn rhs_c(&self, cell_x: usize, cell_y: usize) -> f64 {
        let mut result = 0.0;

        // x-flux
        let flux_c_x_m = if cell_x == 1 {
            0.0
        } else {
            self.flux_c_x_m(cell_x, cell_y)
        };

        let flux_c_x_p = if cell_x == self.p.n_interior_cell_1d {
            0.0
        } else {
            self.flux_c_x_p(cell_x, cell_y)
        };

        result += -(flux_c_x_p - flux_c_x_m) / self.p.dx;

        // y-flux
        let flux_c_y_m = if cell_y == 1 {
            0.0
        } else {
            self.flux_c_y_m(cell_x, cell_y)
        };

        let flux_c_y_p = if cell_y == self.p.n_interior_cell_1d {
            0.0
        } else {
            self.flux_c_y_p(cell_x, cell_y)
        };

        result += -(flux_c_y_p - flux_c_y_m) / self.p.dy;

        // reaction terms
        result += -self.p.gamma_c * self.u(Variable::C, cell_x, cell_y);
        result += self.p.gamma_rho * self.u(Variable::RhoBar, cell_x, cell_y);

        // forcing
        if let Some(forces) = &self.p.forces {
            let x = self.x(cell_x);
            let y = self.y(cell_y);
            result += forces.c_force.as_ref()(self.time, x, y, &self.p);
        }

        result
    }

    fn update_ghost_cells(&mut self) {
        // TODO currently zeroth-order extrapolation
        for cell in 1..=self.p.n_interior_cell_1d {
            // bottom
            *self.u_mut(Variable::RhoBar, cell, 0) = self.u(Variable::RhoBar, cell, 1);
            *self.u_mut(Variable::C, cell, 0) = self.u(Variable::C, cell, 1);

            // top
            *self.u_mut(Variable::RhoBar, cell, self.p.n_interior_cell_1d + 1) =
                self.u(Variable::RhoBar, cell, self.p.n_interior_cell_1d);
            *self.u_mut(Variable::C, cell, self.p.n_interior_cell_1d + 1) =
                self.u(Variable::C, cell, self.p.n_interior_cell_1d);

            // left
            *self.u_mut(Variable::RhoBar, 0, cell) = self.u(Variable::RhoBar, 1, cell);
            *self.u_mut(Variable::C, 0, cell) = self.u(Variable::C, 1, cell);

            // right
            *self.u_mut(Variable::RhoBar, self.p.n_interior_cell_1d + 1, cell) =
                self.u(Variable::RhoBar, self.p.n_interior_cell_1d, cell);
            *self.u_mut(Variable::C, self.p.n_interior_cell_1d + 1, cell) =
                self.u(Variable::C, self.p.n_interior_cell_1d, cell);
        }
    }

    fn step_euler_forward_helper(&mut self, dt: f64) {
        self.update_ghost_cells();

        for cell_x in 1..=self.p.n_interior_cell_1d {
            for cell_y in 1..=self.p.n_interior_cell_1d {
                // Subtract one because self.index() currently returns 1..=self.p.n_interior_cell_1d
                let (idx_x, idx_y) = self.index(cell_x, cell_y);
                let (idx_x, idx_y) = (idx_x - 1, idx_y - 1);
                self.rhs_buffer[(Variable::RhoBar as usize, idx_x, idx_y)] = self.rhs_rho_bar(cell_x, cell_y);
                self.rhs_buffer[(Variable::C as usize, idx_x, idx_y)] = self.rhs_c(cell_x, cell_y);
            }
        }

        {
            let n = self.p.n_interior_cell_1d;
            let mut slice = self.data.slice_mut(s![.., 1..=n, 1..=n]);
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

        let n = self.p.n_interior_cell_1d;

        // save current solution
        let u_old = self
            .data
            .slice(s![.., 1..=n, 1..=n])
            .to_owned();

        // an Euler step gives w_1
        self.step_euler_forward_helper(dt);

        // another Euler step gives the brackets in the w_2 eqn
        self.step_euler_forward_helper(dt);

        // calulate w_2
        let w_2 = {
            let temp = self.data.slice(s![.., 1..=n, 1..=n]);
            0.75 * &u_old + 0.25 * &temp
        };

        // assign w_2 back into the main storage
        let mut slice = self.data.slice_mut(s![.., 1..=n, 1..=n]);
        slice.assign(&w_2);

        // another Euler step gives the brackets in the w_3 eqn
        self.step_euler_forward_helper(dt);

        // calulate w_{n+1}
        let w_np1 = {
            let temp = self.data.slice(s![.., 1..=n, 1..=n]);
            (1.0 / 3.0) * &u_old + (2.0 / 3.0) * &temp
        };

        // assign w_{n+1} back into the main storage
        let mut slice = self.data.slice_mut(s![.., 1..=n, 1..=n]);
        slice.assign(&w_np1);
    }
}
