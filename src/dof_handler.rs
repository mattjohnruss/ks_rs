use crate::mesh::*;
use crate::scheme::*;
//use std::collections::HashMap;

#[derive(Debug)]
pub struct DofHandler {
    dof_map: Vec<usize>,
}

impl DofHandler {
    pub fn new<S: Scheme>(mesh: &Mesh1D, _scheme: &S) -> Self {
        let n_dofs = mesh.n_cells() * S::N_DOFS_PER_CELL;
        let mut dof_map = Vec::with_capacity(n_dofs);

        for i in 0..S::N_DOFS_PER_CELL {
            for (j, _) in mesh.cells().enumerate() {
                dof_map.push(i * mesh.n_cells() + j);
            }
        }

        DofHandler { dof_map }
    }
}
