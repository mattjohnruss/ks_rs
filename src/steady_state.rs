use crate::timestepping::ExplicitTimeSteppable;
use ndarray::prelude::*;

#[derive(Clone)]
pub struct SteadyStateDetector {
    buffer: Array1<f64>,
}

impl SteadyStateDetector {
    pub fn new(n_dof: usize) -> Self {
        Self {
            buffer: Array::zeros(n_dof),
        }
    }

    pub fn is_steady_state(&mut self, obj: &impl ExplicitTimeSteppable, threshold: f64) -> bool {
        obj.rhs(self.buffer.view_mut());
        self.buffer.iter().all(|&value| value.abs() < threshold)
    }
}
