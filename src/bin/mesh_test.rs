use ks_rs::mesh::*;
use std::rc::Rc;

struct DofHandler {
    mesh: Rc<Mesh1D>,
}

impl DofHandler {
    fn new(mesh: Rc<Mesh1D>) -> Self {
        DofHandler { mesh }
    }
}

fn main() {
    let mesh = Rc::new(Mesh1D::new_uniform(0.0, 1.0, 2));
    //let mesh = Rc::new(mesh);
    let dof_handler = DofHandler { mesh: mesh.clone() };

    for cell in mesh.cells() {
        dbg!(cell);
    }
}
