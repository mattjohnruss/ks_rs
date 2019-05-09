#![allow(dead_code)]

pub mod first_order {
    pub static FORWARD_1: [(isize, f64); 2] = [(0, -1.0), (1, 1.0)];
    pub static BACKWARD_1: [(isize, f64); 2] = [(0, 1.0), (-1, -1.0)];
}

pub mod second_order {
    pub static CENTRAL_1: [(isize, f64); 2] = [(-1, -0.5), (1, 0.5)];
    pub static FORWARD_1: [(isize, f64); 3] = [(0, -1.5), (1, 2.0), (2, -0.5)];
    pub static BACKWARD_1: [(isize, f64); 3] = [(0, 1.5), (-1, -2.0), (-2, 0.5)];
    pub static CENTRAL_2: [(isize, f64); 3] = [(-1, 1.0), (0, -2.0), (1, 1.0)];
}

pub fn apply<F>(stencil: &[(isize, f64)], i: usize, f: F) -> f64
where
    F: Fn(usize) -> f64,
{
    stencil.iter().fold(0.0, |acc, (k, w)| {
        let cell_offset = (i as isize + k) as usize;
        acc + w * f(cell_offset)
    })
}
