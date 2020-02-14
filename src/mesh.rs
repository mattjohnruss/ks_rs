#![allow(dead_code)]

#[derive(Clone, Debug)]
pub struct Cell<'a> {
    /// The index of the cell. Its vertices are then `idx` and `idx + 1` in the mesh's storage.
    idx: usize,

    /// The mesh that owns this cell
    mesh: &'a Mesh1D,
}

impl<'a> Cell<'a> {
    fn vertices(&self) -> [f64; 2] {
        [self.mesh.vertices[self.idx], self.mesh.vertices[self.idx + 1]]
    }

    fn centre(&self) -> f64 {
        let vertices = self.vertices();
        0.5 * (vertices[0] + vertices[1])
    }

    pub fn neighbour_left(&self) -> Option<Cell> {
        if self.idx >= 1 {
            Some(Cell {
                idx: self.idx - 1,
                mesh: self.mesh,
            })
        } else {
            None
        }
    }

    pub fn neighbour_right(&self) -> Option<Cell> {
        if self.idx <= self.mesh.n_cells - 2 {
            Some(Cell {
                idx: self.idx + 1,
                mesh: self.mesh,
            })
        } else {
            None
        }
    }

    pub fn dx_left(&self) -> f64 {
        let left = match self.neighbour_left() {
            Some(cell) => cell,
            None => panic!("No left neighbour"),
        };
        self.centre() - left.centre()
    }

    pub fn dx_right(&self) -> f64 {
        let left = match self.neighbour_left() {
            Some(cell) => cell,
            None => panic!("No left neighbour"),
        };
        self.centre() - left.centre()
    }
}

#[derive(Debug)]
pub struct CellIter<'a> {
    cell: Cell<'a>,
}

impl<'a> Iterator for CellIter<'a> {
    type Item = Cell<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        let cell = &mut self.cell;
        if cell.idx < cell.mesh.n_cells {
            let result = Some(cell.clone());
            cell.idx += 1;
            result
        }
        else {
            None
        }
    }
}

#[derive(Debug)]
pub struct Mesh1D {
    vertices: Vec<f64>,
    n_cells: usize,
}

impl Mesh1D {
    /// Create an empty mesh
    pub fn new_empty() -> Self {
        Mesh1D {
            vertices: vec![],
            n_cells: 0,
        }
    }

    /// Create a mesh with uniformly spaced vertices and cells
    pub fn new_uniform(start: f64, end: f64, n_cells: usize) -> Self {
        let mut vertices = Vec::with_capacity(n_cells + 1);
        let dx = (end - start) / n_cells as f64;

        for i in 0..=n_cells {
            vertices.push(i as f64 * dx);
        }

        Mesh1D {
            vertices,
            n_cells,
        }
    }

    /// Returns the number of cells in the mesh
    pub fn n_cells(&self) -> usize {
        self.n_cells
    }

    /// Returns an iterator over the cells
    pub fn cells(&self) -> CellIter {
        let cell = Cell { idx: 0, mesh: self };
        CellIter { cell }
    }

    /// Divides the cell with index `idx` in two by adding a vertex at its centre point
    fn refine(&mut self, idx: usize) {
        let v_1 = self.vertices[idx];
        let v_2 = self.vertices[idx + 1];
        let v_new = 0.5 * (v_1 + v_2);
        self.vertices.insert(idx + 1, v_new);
        self.n_cells += 1;
    }

    /// Perform uniform refinements of the mesh the requested number of times
    pub fn refine_uniformly(&mut self, n_refinements: usize) {
        for _ in 0..n_refinements {
            // We must loop by a step of 2, and up to twice the current number of cells
            for idx in (0..2 * self.n_cells).step_by(2) {
                self.refine(idx);
            }
        }
    }

    /// Merges cell with index `idx` with its neighbour to the right by removing their
    /// shared vertex
    fn coarsen(&mut self, idx: usize) {
        // Make sure the vertex we're removing isn't at one of the boundaries
        assert!(idx > 0);
        assert!(idx < self.vertices.len() - 1);

        self.vertices.remove(idx + 1);
        self.n_cells -= 1;
    }

    pub fn cell<'a>(&'a self, idx: usize) -> Cell<'a> {
        Cell { idx, mesh: self }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn empty_mesh() {
        let mesh = Mesh1D::new_empty();
        assert_eq!(mesh.vertices, Vec::<f64>::new());
    }

    #[test]
    fn uniform_mesh_vertices() {
        let mesh = Mesh1D::new_uniform(0.0, 1.0, 2);
        assert_abs_diff_eq!(mesh.vertices[0], 0.0);
        assert_abs_diff_eq!(mesh.vertices[1], 0.5);
        assert_abs_diff_eq!(mesh.vertices[2], 1.0);
    }

    #[test]
    fn uniform_mesh_cells_iter() {
        let mesh = Mesh1D::new_uniform(0.0, 1.0, 2);
        for (i, cell) in mesh.cells().enumerate() {
            let vertices = cell.vertices();
            assert_abs_diff_eq!(mesh.vertices[i], vertices[0]);
            assert_abs_diff_eq!(mesh.vertices[i + 1], vertices[1]);
        }
    }

    #[test]
    fn manual_refinement() {
        let mut mesh = Mesh1D::new_uniform(0.0, 1.0, 2);
        mesh.refine(0);
        assert_abs_diff_eq!(mesh.vertices[0], 0.0);
        assert_abs_diff_eq!(mesh.vertices[1], 0.25);
        assert_abs_diff_eq!(mesh.vertices[2], 0.5);
        assert_abs_diff_eq!(mesh.vertices[3], 1.0);

        mesh.coarsen(1);
        assert_abs_diff_eq!(mesh.vertices[0], 0.0);
        assert_abs_diff_eq!(mesh.vertices[1], 0.25);
        assert_abs_diff_eq!(mesh.vertices[2], 1.0);

        mesh.refine(1);
        assert_abs_diff_eq!(mesh.vertices[0], 0.0);
        assert_abs_diff_eq!(mesh.vertices[1], 0.25);
        assert_abs_diff_eq!(mesh.vertices[2], 0.625);
        assert_abs_diff_eq!(mesh.vertices[3], 1.0);
    }

    #[test]
    fn uniform_refinement() {
        let mut mesh = Mesh1D::new_uniform(0.0, 2.0, 4);
        mesh.refine_uniformly(1);
        assert_abs_diff_eq!(mesh.vertices[0], 0.0);
        assert_abs_diff_eq!(mesh.vertices[1], 0.25);
        assert_abs_diff_eq!(mesh.vertices[2], 0.5);
        assert_abs_diff_eq!(mesh.vertices[3], 0.75);
        assert_abs_diff_eq!(mesh.vertices[4], 1.0);
        assert_abs_diff_eq!(mesh.vertices[5], 1.25);
        assert_abs_diff_eq!(mesh.vertices[6], 1.5);
        assert_abs_diff_eq!(mesh.vertices[7], 1.75);
        assert_abs_diff_eq!(mesh.vertices[8], 2.0);
    }
}
