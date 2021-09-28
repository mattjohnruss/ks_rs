use crate::stencil::{Stencil, first_order, second_order};
use crate::utilities::minmod;
use crate::timestepping::ExplicitTimeSteppable;
use ndarray::prelude::*;
use std::io::prelude::*;
use crate::utilities::IterMinMax;

/// Wrapper type for variable numbers
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Variable(pub usize);

/// Wrapper type for cell indices
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cell(pub usize);

impl Cell {
    /// Get the cell to the left of `self`, i.e. with index decremented by 1
    fn left(&self) -> Cell {
        Cell(self.0 - 1)
    }

    /// Get the cell to the right of `self`, i.e. with index incremented by 1
    fn right(&self) -> Cell {
        Cell(self.0 + 1)
    }
}

/// The faces of 1D cells
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Face {
    East,
    West,
}

/// Parameters for the physical domain
#[derive(Debug, Clone, Copy)]
pub struct DomainParams {
    /// The number of cells in the interior of the domain (excludes ghost cells)
    pub n_cell: usize,
    /// The width of the domain
    pub width: f64,
}

pub enum BoundaryCondition {
    Dirichlet(f64),
    Flux(f64),
}

impl BoundaryCondition {
    fn is_dirichlet(&self) -> bool {
        matches!(self, BoundaryCondition::Dirichlet(_))
    }
}

/// Allows for used-specified functions for diffusivity, advection velocity, reaction terms and
/// forcing. All functions have access to the current variable, cell and public methods of the
/// underlying `Problem`.
pub trait ProblemFunctions: Sized {
    // TODO for a position-dependent diffusivity to be correct, we would need _p/_m versions like
    // for the velocity, since this one function is currently called in both flux_p and flux_m.
    // Works fine if it's a constant (in each cell).
    /// Diffusivity. Defaults to one.
    fn diffusivity(&self, _problem: &Problem1D<Self>, _var: Variable, _cell: Cell) -> f64 {
        1.0
    }

    /// Advection velocity at the `East` face of the given cell. Defaults to zero.
    fn velocity_p_at_midpoint(&self, _problem: &Problem1D<Self>, _var: Variable, _cell: Cell) -> f64 {
        0.0
    }

    /// Advection velocity at the `West` face of the given cell. Defaults to the velocity at the
    /// `East` face of the neighbouring cell.
    fn velocity_m_at_midpoint(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        self.velocity_p_at_midpoint(problem, var, cell.left())
    }

    /// Reaction terms, sampled at the centre of the cell. Defaults to zero.
    fn reactions(&self, _problem: &Problem1D<Self>, _var: Variable, _cell: Cell) -> f64 {
        0.0
    }

    /// Additional forcing terms, sampled at the centre of the cell. Useful for testing
    /// manufactured solutions. Defaults to zero.
    fn forcing(&self, _problem: &Problem1D<Self>, _var: Variable, _cell: Cell) -> f64 {
        0.0
    }

    /// Value of the flux at the left-hand boundary
    fn left_bc(&self, _problem: &Problem1D<Self>, _var: Variable) -> BoundaryCondition {
        BoundaryCondition::Flux(0.0)
    }

    /// Value of the flux at the right-hand boundary
    fn right_bc(&self, _problem: &Problem1D<Self>, _var: Variable) -> BoundaryCondition {
        BoundaryCondition::Flux(0.0)
    }
}

/// Trivial implementation of `ProblemFunctions` that leaves all functions as the default
/// zero, corresponding to a pure diffusion problem with <math>D = 1</math> for all components
#[derive(Copy, Clone)]
pub struct DiffusionZeroFluxBcs;

impl ProblemFunctions for DiffusionZeroFluxBcs { }

/// The main problem type
#[derive(Clone)]
pub struct Problem1D<F> {
    data: Array<f64, Ix2>,
    ghost_data: Array<f64, Ix2>,
    n_variable: usize,
    pub time: f64,
    pub variable_names: Vec<String>,
    pub domain: DomainParams,
    pub n_dof: usize,
    dx: f64,
    pub functions: F,
}

impl<F> Problem1D<F>
    where F: ProblemFunctions 
{
    /// Create a new `Problem1D` with the given number of variables and cells
    pub fn new(n_variable: usize, domain: DomainParams, functions: F) -> Self {
        let variable_names = (0..n_variable)
            .map(|i| format!("var_{}", i))
            .collect();

        let n_dof = n_variable * domain.n_cell;

        Problem1D {
            data: Array::zeros((n_variable, domain.n_cell)),
            ghost_data: Array::zeros((n_variable, 2)),
            n_variable,
            time: 0.0,
            variable_names,
            domain,
            n_dof,
            dx: domain.width / domain.n_cell as f64,
            functions,
        }
    }

    /// Output the header, consisting of space-separated variable names, to the given
    /// writer
    fn output_header(&self, buffer: &mut impl Write) -> std::io::Result<()> {
        buffer.write_all("t x".as_bytes())?;
        for name in &self.variable_names {
            buffer.write_all(format!(" {}", name).as_bytes())?;
        }
        buffer.write_all("\n".as_bytes())?;
        Ok(())
    }

    /// Generic function that outputs the data produced by the closure `f` to the given writer
    fn output_data(&self, buffer: &mut impl Write, f: impl Fn(Variable, f64, Cell, Face) -> f64) -> std::io::Result<()> {
        let time = self.time;

        for cell in self.interior_cells() {
            // Values at West face
            {
                let x_m = self.x(cell) - 0.5 * self.dx;
                buffer.write_all(format!("{:.6e} {:.6e}", time, x_m).as_bytes())?;
                for var in 0..self.n_variable {
                    let value = f(Variable(var), x_m, cell, Face::West);
                    buffer.write_all(format!(" {:.6e}", value).as_bytes())?;
                }
            }

            buffer.write_all(b"\n")?;

            // Values at East face
            {
                let x_p = self.x(cell) + 0.5 * self.dx;
                buffer.write_all(format!("{:.6e} {:.6e}", time, x_p).as_bytes())?;
                for var in 0..self.n_variable {
                    let value = f(Variable(var), x_p, cell, Face::East);
                    buffer.write_all(format!(" {:.6e}", value).as_bytes())?;
                }
            }

            buffer.write_all(b"\n\n")?;
        }
        Ok(())
    }

    /// Output the values at face boundaries to the given writer
    fn output_face_data(&self, buffer: &mut impl Write) -> std::io::Result<()> {
        self.output_data(buffer, |var, _x, cell, face| self.var_point_value_at_face(var, cell, face))
    }

    /// Output the cell averages to the given writer
    fn output_cell_average_data(&self, buffer: &mut impl Write) -> std::io::Result<()> {
        self.output_data(buffer, |var, _x, cell, _face| self.var(var, cell))
    }

    /// Output the result of a function to the given writer
    fn output_fn_data(&self, buffer: &mut impl Write, f: impl Fn(Variable, f64, f64) -> f64) -> std::io::Result<()> {
        self.output_data(buffer, |var, x, _cell, _face| f(var, x, self.time))
    }

    /// Output the current state of the problem to the given writer
    pub fn output(&self, buffer: &mut impl Write) -> std::io::Result<()> {
        self.output_header(buffer)?;
        self.output_face_data(buffer)?;
        Ok(())
    }

    /// Output the current state of the problem to the given writer
    pub fn output_cell_averages(&self, buffer: &mut impl Write) -> std::io::Result<()> {
        self.output_header(buffer)?;
        self.output_cell_average_data(buffer)?;
        Ok(())
    }

    /// Output the result of a function to the given writer
    pub fn output_fn(&self, buffer: &mut impl Write, f: impl Fn(Variable, f64, f64) -> f64) -> std::io::Result<()> {
        self.output_header(buffer)?;
        self.output_fn_data(buffer, f)?;
        Ok(())
    }

    /// Calculate the integral of `var` over `cell` using a simple rectangle rule. There's
    /// no advantage to using e.g. the trapezium rule since the point value at the centre
    /// of each cell already represents the cell average of the variable and the point
    /// values at the cell edges are computed so as to preserve the cell average anyway.
    pub fn integrate_cell(&self, var: impl Into<Variable>, cell: Cell) -> f64 {
        self.dx * self.var(var.into(), cell)
    }

    /// Calculate the integral of `var` over the whole domain using the rectangle rule on
    /// each cell
    pub fn integrate_solution(&self, var: impl Into<Variable>) -> f64 {
        let var = var.into();
        self.interior_cells()
            .fold(0.0, |total, cell| total + self.integrate_cell(var, cell))
    }

    /// Set the names of the variables. Used in the headers of output files etc. Returns `Err(...)`
    /// if the number of variable names given doesn't match the number of variables
    pub fn set_variable_names(&mut self, names: &[&str]) -> Result<(), String> {
        let n = names.len();
        if names.len() == self.n_variable {
            self.variable_names = names.iter().map(|&s| String::from(s)).collect();
            Ok(())
        } else {
            Err(format!("{} variable names given, but there are {} variables", n, self.n_variable))
        }
    }
    
    /// Get the x-coordinate of the centre of the given cell
    pub fn x(&self, cell: Cell) -> f64 {
        (cell.0 as f64 - 0.5) * self.dx
    }

    /// Get the width of the given cell
    pub fn dx(&self, _cell: Cell) -> f64 {
        // Currently, only uniform meshes are supported
        self.dx
    }

    /// Get the value of the given variable at the centre of the given cell
    #[inline]
    pub fn var(&self, var: impl Into<Variable>, cell: Cell) -> f64 {
        let idx = self.index(cell);
        let var = var.into();
        match idx {
            idx if idx == 0 => {
                unsafe { *self.ghost_data.uget((var.0, 0)) }
            },
            idx if idx == self.domain.n_cell + 1 => {
                unsafe { *self.ghost_data.uget((var.0, 1)) }
            },
            _ => {
                unsafe { *self.data.uget((var.0, idx - 1)) }
            },
        }
    }

    /// Get a mutable reference to the value of the given variable at the centre of the
    /// given cell
    #[inline]
    pub fn var_mut(&mut self, var: impl Into<Variable>, cell: Cell) -> &mut f64 {
        let idx = self.index(cell);
        let var = var.into();
        match idx {
            idx if idx == 0 => unsafe { self.ghost_data.uget_mut((var.0, 0)) },
            idx if idx == self.domain.n_cell + 1 => unsafe { self.ghost_data.uget_mut((var.0, 1)) },
            _ => unsafe { self.data.uget_mut((var.0, idx - 1)) },
        }
    }

    /// Convert from a `Cell` to a `usize` index
    fn index(&self, cell: Cell) -> usize {
        cell.0
    }

    /// Calculate all terms on the right-hand side of the equation for the given variable
    /// in the given cell
    #[inline(always)]
    fn rhs(&self, var: Variable, cell: Cell) -> f64 {
        let mut result = 0.0;

        // flux (includes diffusion and advection)
        result += {
            // This logic also works for equations with zero flux (i.e. ODEs parametrised by
            // space; e.g. C_b in the chemokine model) if the diffusivity, advection velocity and
            // flux BCs are set to zero
            let flux_m = if cell == Cell(1) {
                // At the left-hand boundary
                match self.functions.left_bc(self, var) {
                    BoundaryCondition::Flux(flux) => flux,
                    BoundaryCondition::Dirichlet(_) => self.flux_m_for_dirichlet_bc(var, cell),
                }
            } else {
                self.flux_m(var, cell)
            };

            let flux_p = if cell == Cell(self.domain.n_cell) {
                // At the right-hand boundary
                match self.functions.right_bc(self, var) {
                    BoundaryCondition::Flux(flux) => flux,
                    BoundaryCondition::Dirichlet(_) => self.flux_p_for_dirichlet_bc(var, cell),
                }
            } else {
                self.flux_p(var, cell)
            };

            -(flux_p - flux_m) / self.dx
        };

        // reactions
        result += self.functions.reactions(self, var, cell);

        // forcing
        result += self.functions.forcing(self, var, cell);

        result
    }

    fn flux_p(&self, var: Variable, cell: Cell) -> f64 {
        let velocity_p = self.functions.velocity_p_at_midpoint(self, var, cell);
        let dvar_dx_p = self.dvar_dx_p_at_midpoint(var, cell);
        self.var_point_value_x_p(var, cell) * velocity_p - self.functions.diffusivity(self, var, cell) * dvar_dx_p
    }

    fn flux_m(&self, var: Variable, cell: Cell) -> f64 {
        let velocity_m = self.functions.velocity_m_at_midpoint(self, var, cell);
        let dvar_dx_m = self.dvar_dx_m_at_midpoint(var, cell);
        self.var_point_value_x_m(var, cell) * velocity_m - self.functions.diffusivity(self, var, cell) * dvar_dx_m
    }

    // Flux functions for Dirichlet boundary conditions. These don't upwind the point values,
    // instead getting the appropriate point value directly, since upwinding could cause an
    // out-of-bounds array access when used at a boundary.
    fn flux_p_for_dirichlet_bc(&self, var: Variable, cell: Cell) -> f64 {
        let velocity_p = self.functions.velocity_p_at_midpoint(self, var, cell);
        let dvar_dx_p = self.dvar_dx_p_at_midpoint(var, cell);
        let var_point_value_at_face = self.var_point_value_at_face(var, cell, Face::East);
        var_point_value_at_face * velocity_p - self.functions.diffusivity(self, var, cell) * dvar_dx_p
    }

    fn flux_m_for_dirichlet_bc(&self, var: Variable, cell: Cell) -> f64 {
        let velocity_m = self.functions.velocity_m_at_midpoint(self, var, cell);
        let dvar_dx_m = self.dvar_dx_m_at_midpoint(var, cell);
        let var_point_value_at_face = self.var_point_value_at_face(var, cell, Face::West);
        var_point_value_at_face * velocity_m - self.functions.diffusivity(self, var, cell) * dvar_dx_m
    }

    #[inline]
    pub fn dvar_dx_p_at_midpoint(&self, var: Variable, cell: Cell) -> f64 {
        first_order::Forward1::apply(cell.0, |i| {
            self.var(var, Cell(i)) / self.dx
        })
    }

    #[inline]
    pub fn dvar_dx_m_at_midpoint(&self, var: Variable, cell: Cell) -> f64 {
        self.dvar_dx_p_at_midpoint(var, cell.left())
    }

    #[inline]
    fn var_point_value_x_p(&self, var: Variable, cell: Cell) -> f64 {
        match self.functions.velocity_p_at_midpoint(self, var, cell) {
            v if v > 0.0 => self.var_point_value_at_face(var, cell, Face::East),
            _ => self.var_point_value_at_face(var, cell.right(), Face::West),
        }
    }

    #[inline]
    fn var_point_value_x_m(&self, var: Variable, cell: Cell) -> f64 {
        match self.functions.velocity_m_at_midpoint(self, var, cell) {
            v if v > 0.0 => self.var_point_value_at_face(var, cell.left(), Face::East),
            _ => self.var_point_value_at_face(var, cell, Face::West),
        }
    }

    // We omit the dx factor here since it cancels with the 1/dx in dvar_limited* (also omitted).
    // Chooses the appropriate form of the gradient depending on the boundary conditions.
    #[inline]
    pub fn var_point_value_at_face(&self, var: Variable, cell: Cell, face: Face) -> f64 {
        let value = self.var(var, cell);
        let sign = match face {
            Face::East => 1.0,
            Face::West => -1.0,
        };

        if (cell == Cell(1) && self.functions.left_bc(self, var).is_dirichlet())
            || (cell == Cell(self.domain.n_cell) && self.functions.right_bc(self, var).is_dirichlet()) {
            value + sign * 0.5 * self.dvar_limited_for_dirichlet_bcs(var, cell)
        } else {
            value + sign * 0.5 * self.dvar_limited(var, cell)
        }
    }

    // We omit the dx factor here since it cancels with the 1/dx in dvar_limited_for_dirichlet_bcs
    // (also omitted). This version assumes that `cell` is on a boundary and that Dirichlet BCs
    // are applied at that boundary. This separate function is needed so that Dirichlet BCs
    // themselves can depend on values at the boundary without causing infinite function recursion.
    #[inline]
    pub fn var_point_value_at_face_for_dirichlet_bcs(&self, var: Variable, cell: Cell, face: Face) -> f64 {
        debug_assert!(cell == Cell(1) || cell == Cell(self.domain.n_cell), "Should only be called on the first or last cell");

        let value = self.var(var, cell);
        let sign = match face {
            Face::East => 1.0,
            Face::West => -1.0,
        };
        value + sign * 0.5 * self.dvar_limited_for_dirichlet_bcs(var, cell)
    }

    // We omit the 1/dx factor on the returned value here since it cancels the dx in
    // var_point_value_at_face (also omitted there, the only place this function is called from).
    fn dvar_limited(&self, var: Variable, cell: Cell) -> f64 {
        let dvar_central = second_order::Central1::apply(cell.0, |i| {
            self.var(var, Cell(i))
        });

        let test_p = self.var(var, cell) + 0.5 * dvar_central;
        let test_m = self.var(var, cell) - 0.5 * dvar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            dvar_central
        } else {
            let dvar_forward = first_order::Forward1::apply(cell.0, |i| {
                self.var(var, Cell(i))
            });

            let dvar_backward = first_order::Backward1::apply(cell.0, |i| {
                self.var(var, Cell(i))
            });

            minmod(&[
                   2.0 * dvar_forward,
                   dvar_central,
                   2.0 * dvar_backward,
            ])
        }
    }

    // We omit the 1/dx factor on the returned value here since it cancels the dx in
    // var_point_value_at_face* (also omitted there, the only places this function is called from).
    fn dvar_limited_for_dirichlet_bcs(&self, var: Variable, cell: Cell) -> f64 {
        let dvar_central = second_order::Central1::apply(cell.0, |i| {
            self.var(var, Cell(i))
        });

        let test_p = self.var(var, cell) + 0.5 * dvar_central;
        let test_m = self.var(var, cell) - 0.5 * dvar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            dvar_central
        } else {
            match cell {
                cell if cell == Cell(1) => {
                    let dvar_backward = first_order::Backward1::apply(cell.0, |i| {
                        self.var(var, Cell(i))
                    });

                    if dvar_central > 0.0 && 2.0 * dvar_backward > 0.0 {
                        (dvar_central).min(2.0 * dvar_backward)
                    } else if dvar_central < 0.0 && 2.0 * dvar_backward < 0.0 {
                        (dvar_central).max(2.0 * dvar_backward)
                    } else {
                        0.0
                    }
                }
                cell if cell == Cell(self.domain.n_cell) => {
                    let dvar_forward = first_order::Forward1::apply(cell.0, |i| {
                        self.var(var, Cell(i))
                    });

                    if 2.0 * dvar_forward > 0.0 && dvar_central > 0.0 {
                        (2.0 * dvar_forward).min(dvar_central)
                    } else if 2.0 * dvar_forward < 0.0 && dvar_central < 0.0 {
                        (2.0 * dvar_forward).max(dvar_central)
                    } else {
                        0.0
                    }
                }
                _ => panic!("Should only be called on the first or last cell")
            }
        }
    }

    // Calculates the appropriate value for the ghost cell so that when the flux is
    // calculated, the desired Dirichlet boundary condition is effectively applied.
    fn ghost_cell_value_for_dirichlet_bc(
        &self,
        var: Variable,
        cell: Cell,
        dirichlet_value: f64
    ) -> Option<f64> {
        let dvar_central = 0.5 * (self.var(var, cell.right()) - self.var(var, cell.left()));

        let test_p = self.var(var, cell) + 0.5 * dvar_central;
        let test_m = self.var(var, cell) - 0.5 * dvar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            Some(match cell {
                cell if cell == Cell(1) => {
                    4.0 * (dirichlet_value - self.var(var, cell)) + self.var(var, cell.right())
                }
                cell if cell == Cell(self.domain.n_cell) => {
                    4.0 * (dirichlet_value - self.var(var, cell)) + self.var(var, cell.left())
                }
                _ => panic!("Should only be called on the first or last cell")
            })
        } else {
            match cell {
                cell if cell == Cell(1) => {
                    let dvar_backward = first_order::Backward1::apply(cell.0, |i| {
                        self.var(var, Cell(i))
                    });

                    if dvar_central > 0.0 && 2.0 * dvar_backward > 0.0 {
                        // if both are positive, use the smallest
                        Some(if dvar_central < 2.0 * dvar_backward {
                            // use dvar_central
                            4.0 * (dirichlet_value - self.var(var, cell)) + self.var(var, cell.left())
                        } else {
                            // use 2.0 * dvar_backward
                            dirichlet_value
                        })
                    } else if dvar_central < 0.0 && 2.0 * dvar_backward < 0.0 {
                        // if both are negative, use the largest (smallest in magnitude)
                        Some(if dvar_central > 2.0 * dvar_backward {
                            // use dvar_central
                            4.0 * (dirichlet_value - self.var(var, cell)) + self.var(var, cell.left())
                        } else {
                            // use 2.0 * dvar_backward
                            dirichlet_value
                        })
                    } else {
                        // else they have different signs so use gradient 0.0
                        None
                    }
                }
                cell if cell == Cell(self.domain.n_cell) => {
                    let dvar_forward = first_order::Forward1::apply(cell.0, |i| {
                        self.var(var, Cell(i))
                    });

                    if 2.0 * dvar_forward > 0.0 && dvar_central > 0.0 {
                        // if both are positive, use the smallest
                        Some(if 2.0 * dvar_forward < dvar_central {
                            // use 2.0 * dvar_forward
                            dirichlet_value
                        } else {
                            // use dvar_central
                            4.0 * (dirichlet_value - self.var(var, cell)) + self.var(var, cell.left())
                        })
                    } else if 2.0 * dvar_forward < 0.0 && dvar_central < 0.0 {
                        // else if both are negative, use the largest (smallest in magnitude)
                        Some(if 2.0 * dvar_forward > dvar_central {
                            // use 2.0 * dvar_forward
                            dirichlet_value
                        } else {
                            // use dvar_central
                            4.0 * (dirichlet_value - self.var(var, cell)) + self.var(var, cell.left())
                        })
                    } else {
                        // else they have different signs so use gradient 0.0
                        None
                    }
                }
                _ => panic!("Should only be called on the first or last cell")
            }
        }
    }

    /// Update the values in the ghost cells. Automatically called before/after timestepping as
    /// apprappropriate, but must be called manually after modifying any values (e.g. setting
    /// initial conditions) for the piecewise linear reconstruction of the solution to be correct.
    pub fn update_ghost_cells(&mut self) {
        for var in 0..self.n_variable {
            let var = Variable(var);

            match self.functions.left_bc(self, var) {
                BoundaryCondition::Flux(_) => {
                    *self.var_mut(var, Cell(0)) = self.var(var, Cell(1));
                }
                BoundaryCondition::Dirichlet(dirichlet_value) => {
                    let ghost_cell_value = self.ghost_cell_value_for_dirichlet_bc(var, Cell(1), dirichlet_value);
                    match ghost_cell_value {
                        Some(ghost_cell_value) => *self.var_mut(var, Cell(0)) = ghost_cell_value,
                        None => *self.var_mut(var, Cell(1)) = dirichlet_value,
                    }
                }
            }

            match self.functions.right_bc(self, var) {
                BoundaryCondition::Flux(_) => {
                    *self.var_mut(var, Cell(self.domain.n_cell + 1)) =
                        self.var(var, Cell(self.domain.n_cell));
                }
                BoundaryCondition::Dirichlet(dirichlet_value) => {
                    let ghost_cell_value = self.ghost_cell_value_for_dirichlet_bc(var, Cell(self.domain.n_cell), dirichlet_value);
                    match ghost_cell_value {
                        Some(ghost_cell_value) => *self.var_mut(var, Cell(self.domain.n_cell + 1)) = ghost_cell_value,
                        None => *self.var_mut(var, Cell(self.domain.n_cell)) = dirichlet_value,
                    }
                }
            }
        }
    }

    pub fn interior_cells(&self) -> impl Iterator<Item=Cell> {
        (0..self.domain.n_cell).map(|i| Cell(i + 1))
    }

    /// Calculates the largest timestep such that the solution remains positive after a forward
    /// Euler step. The method used is a generalization of that in Chertock et al. (2018) and
    /// related papers allowing for arbitrary velocities and reactions.
    ///
    /// It currently does not implement the extra factors that arise when an arbitrary function
    /// g(rho) of each variable is allowed in the advective terms rather than the usual rho.
    pub fn calculate_dt(&self) -> f64 {
        // NOTE: We could avoid repeated allocations by putting these Vecs in the Problem struct
        // (using interior mutability). Initial testing suggests it's actually slower than the
        // naive approach, presumably because these allocations are optimised away. Similarly,
        // avoiding accumulating all a, k and d entries by keeping track of just their running
        // minimum sounds like it should be faster than this, but it appears also not to be.

        // NOTE: Entries of the vectors a, k, d will be infinity if their denominator is zero, e.g.
        // when a velocity is zero for a particular variable. This is the expected behaviour when a
        // particular term doesn't contribute to the timestep restriction. This isn't a problem
        // provided that at least one of the equations has, say, a nonzero diffusion coefficient.
        // In the end the threshold on dt depends on the smallest of these quantities, and clearly
        // min(x, inf) = x, provided x != inf.

        // Advective terms for each variable
        let mut a = Vec::with_capacity(self.n_variable);

        // Reaction terms for each variable
        let mut k = Vec::with_capacity(self.n_variable);

        // Diffusion terms for each variable
        let mut d = Vec::with_capacity(self.n_variable);

        for var in 0..self.n_variable {
            let var = Variable(var);

            // TODO: These 3 loops over the cells could be combined into one.

            let a_i_denom = self
                .interior_cells()
                .map(|cell| self.functions.velocity_p_at_midpoint(self, var, cell).abs())
                .max_all();

            let k_i_denom = self
                .interior_cells()
                .map(|cell| {
                    2.0 * self.functions.diffusivity(self, var, cell) / self.dx.powi(2)
                        - self.functions.reactions(self, var, cell) / self.var(var, cell)
                })
                .max_all();

            let diffusivity_max = self
                .interior_cells()
                .map(|cell| self.functions.diffusivity(self, var, cell))
                .max_all();

            a.push(self.dx / (8.0 * a_i_denom));
            k.push(1.0 / k_i_denom);
            d.push(0.25 * self.dx.powi(2) / diffusivity_max);
        }

        let a_min = a.iter().min_all();
        let k_min = k.iter().min_all();
        let d_min = d.iter().min_all();

        [a_min, k_min, d_min].iter().min_all()
    }
}

impl<F> ExplicitTimeSteppable for Problem1D<F>
    where F: ProblemFunctions
{
    fn time(&self) -> f64 {
        self.time
    }

    fn time_mut(&mut self) -> &mut f64 {
        &mut self.time
    }

    fn dofs(&self) -> ArrayView1<f64> {
        self.data.view().into_shape(self.data.len()).unwrap()
    }

    fn dofs_mut(&mut self) -> ArrayViewMut1<f64> {
        let len = self.data.len();
        self.data.view_mut().into_shape(len).unwrap()
    }

    fn rhs(&self, buffer: ArrayViewMut1<f64>) {
        let mut buffer = buffer.into_shape((self.n_variable, self.domain.n_cell)).unwrap();

        (0..self.n_variable).for_each(|var| {
            self.interior_cells().for_each(|cell| {
                let idx = cell.0 - 1;
                buffer[(var, idx)] = self.rhs(Variable(var), cell);
            });
        });
    }

    fn actions_before_explicit_stage(&mut self) {
        self.update_ghost_cells();
    }

    fn actions_after_explicit_timestep(&mut self) {
        self.update_ghost_cells();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn headers() {
        use std::io::BufWriter;

        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let functions = DiffusionZeroFluxBcs { };

        // Zero variable edge case
        let problem = Problem1D::new(0, domain, functions);
        let mut header = Vec::new();
        {
            let mut header_writer = BufWriter::new(&mut header);
            problem.output_header(&mut header_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&header).unwrap(), "t x\n");

        // A single variable
        let problem = Problem1D::new(1, domain, functions);
        let mut header = Vec::new();
        {
            let mut header_writer = BufWriter::new(&mut header);
            problem.output_header(&mut header_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&header).unwrap(), "t x var_0\n");

        // Multiple variables
        let problem = Problem1D::new(3, domain, functions);
        let mut header = Vec::new();
        {
            let mut header_writer = BufWriter::new(&mut header);
            problem.output_header(&mut header_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&header).unwrap(), "t x var_0 var_1 var_2\n");
    }

    #[test]
    fn variable_names() {
        use std::io::BufWriter;

        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let functions = DiffusionZeroFluxBcs { };

        // Zero variable edge case
        let mut problem = Problem1D::new(0, domain, functions);
        let result = problem.set_variable_names(&[]);
        assert!(result.is_ok());

        let mut header = Vec::new();
        {
            let mut header_writer = BufWriter::new(&mut header);
            problem.output_header(&mut header_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&header).unwrap(), "t x\n");

        // 3 variables, 3 names
        let mut problem = Problem1D::new(3, domain, functions);
        let result = problem.set_variable_names(&["a", "b", "c"]);
        assert!(result.is_ok());

        let mut header = Vec::new();
        {
            let mut header_writer = BufWriter::new(&mut header);
            problem.output_header(&mut header_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&header).unwrap(), "t x a b c\n");
    }

    #[test]
    fn variable_names_fail() {
        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let functions = DiffusionZeroFluxBcs { };

        // 3 variables, 2 names
        let mut problem = Problem1D::new(3, domain, functions);
        let result = problem.set_variable_names(&["a", "b"]);
        assert_eq!(result, Err("2 variable names given, but there are 3 variables".to_string()));

        // 3 variables, 4 names
        let mut problem = Problem1D::new(3, domain, functions);
        let result = problem.set_variable_names(&["a", "b", "c", "d"]);
        assert_eq!(result, Err("4 variable names given, but there are 3 variables".to_string()));
    }

    #[test]
    fn output() {
        use std::io::BufWriter;

        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let functions = DiffusionZeroFluxBcs { };

        let problem = Problem1D::new(3, domain, functions);
        let mut output = Vec::new();
        {
            let mut output_writer = BufWriter::new(&mut output);
            problem.output_face_data(&mut output_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&output).unwrap(), "0.000000e0 0.000000e0 0.000000e0 0.000000e0 0.000000e0\n0.000000e0 1.000000e0 0.000000e0 0.000000e0 0.000000e0\n\n");
    }

    #[test]
    fn output_cell_averages() {
        use std::io::BufWriter;

        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let functions = DiffusionZeroFluxBcs { };

        let problem = Problem1D::new(3, domain, functions);
        let mut output = Vec::new();
        {
            let mut output_writer = BufWriter::new(&mut output);
            problem.output_cell_average_data(&mut output_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&output).unwrap(), "0.000000e0 0.000000e0 0.000000e0 0.000000e0 0.000000e0\n0.000000e0 1.000000e0 0.000000e0 0.000000e0 0.000000e0\n\n");
    }

    #[test]
    fn x_coords() {
        let functions = DiffusionZeroFluxBcs { };

        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let problem = Problem1D::new(1, domain, functions);
        assert_abs_diff_eq!(problem.x(Cell(1)), 0.5);

        let domain = DomainParams { n_cell: 100, width: 8.5 };
        let problem = Problem1D::new(1, domain, functions);
        assert_abs_diff_eq!(problem.x(Cell(12)), 0.9775);
    }

    #[test]
    fn cells() {
        let cell = Cell(13);
        assert_eq!(cell.left(), Cell(12));
        assert_eq!(cell.right(), Cell(14));
    }

    #[test]
    fn cells_iter() {
        let functions = DiffusionZeroFluxBcs { };
        let domain = DomainParams { n_cell: 10, width: 1.0 };
        let problem = Problem1D::new(1, domain, functions);
        for (idx, cell) in (1..=10).zip(problem.interior_cells()) {
            assert_eq!(cell.0, idx);
        }
    }
}
