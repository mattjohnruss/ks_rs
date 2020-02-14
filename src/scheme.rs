use crate::mesh::Cell;

pub trait Scheme {
    const N_DOFS_PER_CELL: usize;

    fn rhs(&self, cell: &Cell, local_dof: usize) -> f64;
}
