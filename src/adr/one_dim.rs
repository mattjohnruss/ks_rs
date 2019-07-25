use crate::stencil::{Stencil, first_order, second_order};
use crate::utilities::minmod;
use crate::timestepping::ExplicitTimeSteppable;
use ndarray::prelude::*;
use std::io::prelude::*;

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
enum Face {
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
    fn velocity_m_at_midpoint(&self, _problem: &Problem1D<Self>, _var: Variable, _cell: Cell) -> f64 {
        self.velocity_p_at_midpoint(_problem, _var, _cell.left())
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
}

/// Trivial implementation of `ProblemFunctions` that leaves all functions as the default
/// zero, corresponding to a pure diffusion problem with <math>D = 1</math> for all components
#[derive(Copy, Clone)]
pub struct DiffusionOnly;

impl ProblemFunctions for DiffusionOnly { }

/// The main problem type
pub struct Problem1D<F> {
    data: Array<f64, Ix2>,
    ghost_data: Array<f64, Ix2>,
    // TODO: rename to n_var?
    pub n_variable: usize,
    pub time: f64,
    // TODO: rename to var_names?
    // TODO: change to HashMap<Variable, String>?
    pub variable_names: Vec<String>,
    pub domain: DomainParams,
    dx: f64,
    functions: F,
}

impl<F> Problem1D<F>
    where F: ProblemFunctions 
{
    /// Create a new `Problem1D` with the given number of variables and cells
    pub fn new(n_variable: usize, domain: DomainParams, functions: F) -> Self {
        let variable_names = (0..n_variable)
            .map(|i| format!("var_{}", i))
            .collect();

        Problem1D {
            data: Array::zeros((n_variable, domain.n_cell)),
            ghost_data: Array::zeros((n_variable, 2)),
            n_variable,
            time: 0.0,
            variable_names,
            domain,
            dx: domain.width / domain.n_cell as f64,
            functions,
        }
    }

    /// Output the header, consisting of space-separated variable names, to the given
    /// writer
    fn output_header<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        buffer.write_all("t x".as_bytes())?;
        for name in &self.variable_names {
            buffer.write_all(format!(" {}", name).as_bytes())?;
        }
        buffer.write_all("\n".as_bytes())?;
        Ok(())
    }

    /// Output the data to the given writer
    fn output_data<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        let time = self.time;
        
        for cell in 1..=self.domain.n_cell {
            let cell = Cell(cell);
            let x = self.x(cell);
            buffer.write_all(format!("{:.6e} {:.6e}", time, x).as_bytes())?;

            for var in self.all_vars(cell) {
                buffer.write_all(format!(" {:.6e}", var).as_bytes())?;

            //for var in 0..self.n_variable {
                //let value = self.var(Variable(var), cell);
                //buffer.write_all(format!(" {:.6e}", value).as_bytes())?;
            }

            buffer.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Output the current state of the problem to the given writer
    pub fn output<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        self.output_header(buffer)?;
        self.output_data(buffer)?;
        Ok(())
    }
    
    /// Get the x-coordinate of the centre of the given cell
    pub fn x(&self, cell: Cell) -> f64 {
        (cell.0 as f64 - 0.5) * self.dx
    }

    /// Get the value of the given variable at the centre of the given cell
    // TODO: rename, but what to? u() is possibly confusing, since in KS u = dC_dx; var()
    // conflicts with the var parameters used everywhere; maybe value(), value_at_centre()?
    pub fn var(&self, var: Variable, cell: Cell) -> f64 {
        let idx = self.index(cell);
        match idx {
            idx if idx == 0 => {
                self.ghost_data[(var.0, 0)]
            },
            idx if idx == self.domain.n_cell + 1 => {
                self.ghost_data[(var.0, 1)]
            },
            _ => {
                self.data[(var.0, idx - 1)]
            },
        }
    }

    /// Get a mutable reference to the value of the given variable at the centre of the
    /// given cell
    pub fn var_mut(&mut self, var: Variable, cell: Cell) -> &mut f64 {
        let idx = self.index(cell);
        match idx {
            idx if idx == 0 => &mut self.ghost_data[(var.0, 0)],
            idx if idx == self.domain.n_cell + 1 => &mut self.ghost_data[(var.0, 1)],
            _ => &mut self.data[(var.0, idx - 1)],
        }
    }

    /// Get an ArrayView for all variables in the given cell
    fn all_vars(&self, cell: Cell) -> ArrayView1<f64> {
        let idx = self.index(cell) - 1;
        self.data.slice(s![.., idx])
    }

    /// Convert from a `Cell` to a `usize` index
    fn index(&self, cell: Cell) -> usize {
        cell.0
    }

    /// Calculate all terms on the right-hand side of the equation for the given variable
    /// in the given cell
    fn rhs(&self, var: Variable, cell: Cell) -> f64 {
        let mut result = 0.0;

        // flux (includes diffusion and advection)
        result += {
            let flux_m = if cell == Cell(1) {
                0.0
            } else {
                self.flux_m(var, cell)
            };

            let flux_p = if cell == Cell(self.domain.n_cell) {
                0.0
            } else {
                self.flux_p(var, cell)
            };

            -(flux_p - flux_m) / self.dx
        };

        // reactions
        result += self.functions.reactions(&self, var, cell);

        // forcing
        result += self.functions.forcing(&self, var, cell);

        result
    }

    fn flux_p(&self, var: Variable, cell: Cell) -> f64 {
        let velocity_p = self.functions.velocity_p_at_midpoint(&self, var, cell);
        let dvar_dx_p = self.dvar_dx_p_at_midpoint(var, cell);
        self.var_point_value_x_p(var, cell) * velocity_p - self.functions.diffusivity(&self, var, cell) * dvar_dx_p
    }

    fn flux_m(&self, var: Variable, cell: Cell) -> f64 {
        let velocity_m = self.velocity_m_at_midpoint(var, cell);
        let dvar_dx_m = self.dvar_dx_m_at_midpoint(var, cell);
        self.var_point_value_x_m(var, cell) * velocity_m - self.functions.diffusivity(&self, var, cell) * dvar_dx_m
    }

    fn velocity_m_at_midpoint(&self, var: Variable, cell: Cell) -> f64 {
        self.functions.velocity_p_at_midpoint(&self, var, cell.left())
    }

    pub fn dvar_dx_p_at_midpoint(&self, var: Variable, cell: Cell) -> f64 {
        first_order::Forward1::apply(cell.0, |i| {
            self.var(var, Cell(i)) / self.dx
        })
    }

    pub fn dvar_dx_m_at_midpoint(&self, var: Variable, cell: Cell) -> f64 {
        self.dvar_dx_p_at_midpoint(var, cell.left())
    }

    fn var_point_value_x_p(&self, var: Variable, cell: Cell) -> f64 {
        match self.functions.velocity_p_at_midpoint(&self, var, cell) {
            v if v > 0.0 => self.var_point_value_at_face(var, cell, Face::East),
            _ => self.var_point_value_at_face(var, cell.right(), Face::West),
        }
    }

    fn var_point_value_x_m(&self, var: Variable, cell: Cell) -> f64 {
        match self.velocity_m_at_midpoint(var, cell) {
            v if v > 0.0 => self.var_point_value_at_face(var, cell.left(), Face::East),
            _ => self.var_point_value_at_face(var, cell, Face::West),
        }
    }

    fn var_point_value_at_face(&self, var: Variable, cell: Cell, face: Face) -> f64 {
        let value = self.var(var, cell);
        match face {
            Face::East => value + 0.5 * self.dx * self.dvar_dx(var, cell),
            Face::West => value - 0.5 * self.dx * self.dvar_dx(var, cell),
        }
    }

    fn dvar_dx(&self, var: Variable, cell: Cell) -> f64 {
        let dvar_central = second_order::Central1::apply(cell.0, |i| {
            self.var(var, Cell(i))
        });

        let test_p = self.var(var, cell) + 0.5 * dvar_central;
        let test_m = self.var(var, cell) - 0.5 * dvar_central;

        if test_p >= 0.0 && test_m >= 0.0 {
            dvar_central / self.dx
        } else {
            let dvar_forward = first_order::Forward1::apply(cell.0, |i| {
                self.var(var, Cell(i))
            });

            let dvar_backward = first_order::Backward1::apply(cell.0, |i| {
                self.var(var, Cell(i))
            });

            minmod(&[
                   2.0 * dvar_forward / self.dx,
                   dvar_central / self.dx,
                   2.0 * dvar_backward / self.dx,
            ])
        }
    }

    fn update_ghost_cells(&mut self) {
        for var in 0..self.n_variable {
            let var = Variable(var);
            *self.var_mut(var, Cell(0)) = self.var(var, Cell(1));
            *self.var_mut(var, Cell(self.domain.n_cell + 1)) =
                self.var(var, Cell(self.domain.n_cell));
        }
    }

    pub fn interior_cells(&self) -> InteriorCellsIter {
        InteriorCellsIter {
            cell: Cell(1),
            n_cell: self.domain.n_cell,
        }
    }
}

impl<F> ExplicitTimeSteppable for Problem1D<F>
    where F: ProblemFunctions
{
    fn time(&self) -> &f64 {
        &self.time
    }

    fn time_mut(&mut self) -> &mut f64 {
        &mut self.time
    }

    fn dofs(&self) -> ArrayView1<f64> {
        self.data.view().into_shape(self.data.len()).unwrap()
    }

    fn set_dofs(&mut self, dofs: ArrayView1<f64>) {
        let len = self.data.len();
        let mut dofs_mut = self.data.view_mut().into_shape(len).unwrap();
        for (dof, dof_new) in dofs_mut.iter_mut().zip(dofs.iter()) {
            *dof = *dof_new
        }
    }

    fn increment_and_multiply_dofs(&mut self, increment: ArrayView1<f64>, factor: f64) {
        let len = self.data.len();
        let mut dofs_mut = self.data.view_mut().into_shape(len).unwrap();
        for (dof, inc) in dofs_mut.iter_mut().zip(increment.iter()) {
            *dof += factor * *inc;
        }
    }

    fn scale_dofs(&mut self, factor: f64) {
        let len = self.data.len();
        let mut dofs_mut = self.data.view_mut().into_shape(len).unwrap();
        for dof in dofs_mut.iter_mut() {
            *dof *= factor;
        }
    }

    fn rhs(&self) -> Array1<f64> {
        let mut rhs = Array::zeros((self.n_variable, self.domain.n_cell));

        for var in 0..self.n_variable {
            for cell in 1..=self.domain.n_cell {
                rhs[(var, cell - 1)] = self.rhs(Variable(var), Cell(cell));
            }
        }

        let len = self.data.len();
        rhs.into_shape(len).unwrap()
    }

    fn actions_before_explicit_stage(&mut self) {
        self.update_ghost_cells();
    }
}

/// Iterates over the interior cells of a problem
pub struct InteriorCellsIter {
    /// The current cell
    cell: Cell,
    /// The number of interior cells
    n_cell: usize,
}

impl Iterator for InteriorCellsIter {
    type Item = Cell;
    fn next(&mut self) -> Option<Cell> {
        let result = match self.cell {
            Cell(i) if i <= self.n_cell => Some(self.cell),
            _ => None,
        };
        self.cell.0 += 1;
        result
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
        let functions = DiffusionOnly { };

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
    fn output() {
        use std::io::BufWriter;

        let domain = DomainParams { n_cell: 1, width: 1.0 };
        let functions = DiffusionOnly { };

        let problem = Problem1D::new(3, domain, functions);
        let mut output = Vec::new();
        {
            let mut output_writer = BufWriter::new(&mut output);
            problem.output_data(&mut output_writer).unwrap();
        }
        assert_eq!(std::str::from_utf8(&output).unwrap(), "0.000000e0 5.000000e-1 0.000000e0 0.000000e0 0.000000e0\n");
    }

    #[test]
    fn x_coords() {
        let functions = DiffusionOnly { };

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
        let functions = DiffusionOnly { };
        let domain = DomainParams { n_cell: 10, width: 1.0 };
        let problem = Problem1D::new(1, domain, functions);
        for (idx, cell) in (1..10).zip(problem.interior_cells()) {
            assert_eq!(cell.0, idx);
        }
    }
}
