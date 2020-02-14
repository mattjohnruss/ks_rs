use ks_rs::mesh::*;
use ks_rs::dof_handler::*;
use ks_rs::scheme::*;

struct DummyScheme {
}

impl Scheme for DummyScheme {
    const N_DOFS_PER_CELL: usize = 2;

    fn rhs(&self, cell: &Cell, local_dof: usize) -> f64 {
        0.0
    }
}

fn main() {
    let mesh = Mesh1D::new_uniform(0.0, 1.0, 10);
    let scheme = DummyScheme {};

    let dof_handler = DofHandler::new(&mesh, &scheme);

    dbg!(&dof_handler);

    //for cell in mesh.cells() {
        //dbg!(cell);
    //}

    //mesh.refine_uniformly(1);
}
