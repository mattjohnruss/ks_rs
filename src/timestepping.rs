use ndarray::prelude::*;

/// Represent an object whose data can be timestepped with an explicit scheme
/// TODO clean up this trait's interface. E.g., increment_and_multiply_dofs() could be a
/// provided function if say set_dofs() is replaced with dofs_mut(); time() could return by
/// value...
pub trait ExplicitTimeSteppable {
    /// Get the current time
    fn time(&self) -> &f64;

    /// Get a mut reference to the current time
    fn time_mut(&mut self) -> &mut f64;

    /// Get an `ArrayView1` into the current state of the degrees of freedom
    fn dofs(&self) -> ArrayView1<f64>;

    /// Set the object's degrees of freedom from `dofs`
    fn set_dofs(&mut self, dofs: ArrayView1<f64>);

    /// Perform `dofs += factor * increment`
    fn increment_and_multiply_dofs(&mut self, increment: ArrayView1<f64>, factor: f64);

    /// Perform `dofs *= factor`
    fn scale_dofs(&mut self, factor: f64);

    /// Gets an owned array of the current right-hand side values
    fn rhs(&self) -> Array1<f64>;

    /// Actions performed before each explicit timestep
    fn actions_before_explicit_timestep(&mut self) { }

    /// Actions performed before each stage of the scheme
    fn actions_before_explicit_stage(&mut self) { }

    /// Actions performed after each stage of the scheme
    fn actions_after_explicit_stage(&mut self) { }

    /// Actions performed after each explicit timestep
    fn actions_after_explicit_timestep(&mut self) { }
}

/// Trait for explicit timesteppers
pub trait ExplicitTimeStepper {
    /// Take a timestep of length `dt` on the timesteppable object `obj`
    fn step<T: ExplicitTimeSteppable>(&self, obj: &mut T, dt: f64);
}

/// Euler forward timestepper
pub struct EulerForward;

impl EulerForward {
    /// Create a new `EulerForward` timestepper
    pub fn new() -> Self {
        EulerForward {}
    }
}

impl ExplicitTimeStepper for EulerForward {
    fn step<T: ExplicitTimeSteppable>(&self, obj: &mut T, dt: f64) {
        obj.actions_before_explicit_timestep();
        obj.actions_before_explicit_stage();

        obj.increment_and_multiply_dofs(obj.rhs().view(), dt);
        *obj.time_mut() += dt;

        obj.actions_after_explicit_stage();
        obj.actions_after_explicit_timestep();
    }
}

/// Fourth-order, four stage Runge-Kutta timestepper
pub struct RungeKutta44;

impl RungeKutta44 {
    /// Create a new Runge-Kutta timestepper
    pub fn new() -> Self {
        RungeKutta44 {}
    }
}

impl ExplicitTimeStepper for RungeKutta44 {
    fn step<T: ExplicitTimeSteppable>(&self, obj: &mut T, dt: f64) {
        obj.actions_before_explicit_timestep();

        obj.actions_before_explicit_stage();

        let dofs = obj.dofs().to_vec();

        // Stage 1
        let k_1 = obj.rhs().to_vec();
        obj.increment_and_multiply_dofs(aview1(&k_1), 0.5 * dt);
        *obj.time_mut() += 0.5 * dt;
        obj.actions_after_explicit_stage();

        obj.actions_before_explicit_stage();

        // Stage 2
        let k_2 = obj.rhs().to_vec();
        obj.set_dofs(aview1(&dofs));
        obj.increment_and_multiply_dofs(aview1(&k_2), 0.5 * dt);
        obj.actions_after_explicit_stage();

        obj.actions_before_explicit_stage();

        // Stage 3
        let k_3 = obj.rhs().to_vec();
        obj.set_dofs(aview1(&dofs));
        obj.increment_and_multiply_dofs(aview1(&k_3), dt);
        *obj.time_mut() += 0.5 * dt;
        obj.actions_after_explicit_stage();

        obj.actions_before_explicit_stage();

        // Stage 4
        let k_4 = obj.rhs().to_vec();
        obj.set_dofs(aview1(&dofs));

        obj.increment_and_multiply_dofs(aview1(&k_1), dt / 6.0);
        obj.increment_and_multiply_dofs(aview1(&k_2), dt / 3.0);
        obj.increment_and_multiply_dofs(aview1(&k_3), dt / 3.0);
        obj.increment_and_multiply_dofs(aview1(&k_4), dt / 6.0);

        obj.actions_after_explicit_stage();

        obj.actions_after_explicit_timestep();
    }
}

/// Third-order, three stage strong stability-preserving Runge-Kutta timestepper
pub struct SspRungeKutta33;

impl SspRungeKutta33 {
    /// Create a new SSP-RK timestepper
    pub fn new() -> Self {
        SspRungeKutta33 {}
    }
}

impl ExplicitTimeStepper for SspRungeKutta33 {
    // NOTE from Hesthaven & Warburton (2008), p.158
    // NOTE double check the placement of actions_before/after_explicit_stage()
    fn step<T: ExplicitTimeSteppable>(&self, obj: &mut T, dt: f64) {
        obj.actions_before_explicit_timestep();

        // get the current dofs
        let v_0 = obj.dofs().to_vec();

        // work out f(v_0, t)
        obj.actions_before_explicit_stage();
        let rhs_1 = obj.rhs().to_vec();
        // multipliy by dt and add to current dofs (an Euler step)
        obj.increment_and_multiply_dofs(aview1(&rhs_1), dt);
        // don't actually need to save a copy of v_1
        //let v_1 = obj.dofs().to_vec();

        // increment time by dt in prep for the next stage
        *obj.time_mut() += dt;
        obj.actions_after_explicit_stage();

        // work out f(v_1, t + dt)
        obj.actions_before_explicit_stage();
        let rhs_2 = obj.rhs().to_vec();
        // pre: dofs == v_1; post: dofs == v_1 + dt*rhs_2
        obj.increment_and_multiply_dofs(aview1(&rhs_2), dt);
        // pre: dofs == v_1 + dt*rhs_2; post: dofs == 3.0*v_0 + v_1 + dt*rhs_2
        obj.increment_and_multiply_dofs(aview1(&v_0), 3.0);
        // pre: dofs == 3.0*v_0 + v_1 + dt*rhs_2; post: dofs == 0.25*(3.0*v_0 + v_1 + dt*rhs_2)
        obj.scale_dofs(0.25);
        let v_2 = obj.dofs().to_vec();

        // this is strange, but apparently the next rhs is evaluated back half a step...
        *obj.time_mut() -= 0.5 * dt;
        obj.actions_after_explicit_stage();

        // work out f(v_2, t + 0.5*dt)
        obj.actions_before_explicit_stage();
        let rhs_3 = obj.rhs().to_vec();
        // pre: dofs == v_2; post: dofs == 2.0*v_2
        obj.increment_and_multiply_dofs(aview1(&v_2), 1.0);
        // pre: dofs == 2.0*v_2; post: dofs == v_0 + 2.0*v_2
        obj.increment_and_multiply_dofs(aview1(&v_0), 1.0);
        // pre: dofs == v_0 + 2.0*v_2; post: dofs == v_0 + 2.0*v_2 + 2.0*dt*rhs_3
        obj.increment_and_multiply_dofs(aview1(&rhs_3), 2.0 * dt);
        // pre: dofs == v_0 + 2.0*v_2 2.0*dt*rhs_3; post: dofs == 1/3(v_0 + 2.0*v_2 + 2.0*dt*rhs_3)
        obj.scale_dofs(1.0/3.0);

        // set the time to t + dt by going forwards half a step
        *obj.time_mut() += 0.5 * dt;
        obj.actions_after_explicit_stage();

        obj.actions_after_explicit_timestep();
    }
}

// NOTE these tests don't actually have any assertions etc. They just do a test calculation
// and dump the results to file, which I check manually against known solutions.
#[cfg(test)]
mod test {
    use super::*;
    use std::fs;
    use std::io::{BufWriter, Write};

    struct DummyProblem {
        time: f64,
        data: Array1<f64>,
        n_dof: usize,
    }

    impl DummyProblem {
        fn new(n_dof: usize) -> Self {
            DummyProblem {
                time: 0.0,
                data: Array::zeros(n_dof),
                n_dof,
            }
        }

        // dx/dt = x
        fn calculate_rhs(&self, i: usize) -> f64 {
            self.data[i]
        }

        fn exact_solution(&self) -> f64 {
            self.time.exp()
        }
    }

    impl ExplicitTimeSteppable for DummyProblem {
        fn time(&self) -> &f64 {
            &self.time
        }

        fn time_mut(&mut self) -> &mut f64 {
            &mut self.time
        }

        fn dofs(&self) -> ArrayView1<f64> {
            self.data.view()
        }

        fn set_dofs(&mut self, dofs: ArrayView1<f64>) {
            self.data.assign(&dofs);
        }

        fn increment_and_multiply_dofs(&mut self, increment: ArrayView1<f64>, factor: f64) {
            self.data += &(factor * &increment);
        }

        fn scale_dofs(&mut self, factor: f64) {
            self.data *= factor;
        }

        fn rhs(&self) -> Array1<f64> {
            let mut rhs = Array1::zeros(self.n_dof);
            for i in 0..self.n_dof {
                rhs[i] = self.calculate_rhs(i);
            }
            rhs
        }
    }

    const N_STEP: usize = 100;
    const DT: f64 = 1.0 / N_STEP as f64;

    #[test]
    fn euler_forward() {
        let file = fs::File::create("euler_forward.csv").unwrap();
        let mut buf_writer = BufWriter::new(file);

        let mut problem = DummyProblem::new(1);
        let timestepper = EulerForward {};

        problem.data[0] = 1.0;
        writeln!(buf_writer, "{} {} {}", problem.time, problem.data[0], problem.exact_solution()).unwrap();

        for _ in 0..N_STEP {
            timestepper.step(&mut problem, DT);
            writeln!(buf_writer, "{} {} {}", problem.time, problem.data[0], problem.exact_solution()).unwrap();
        }
    }

    #[test]
    #[test]
    fn runge_kutta_4_4() {
        let file = fs::File::create("runge_kutta_4_4.csv").unwrap();
        let mut buf_writer = BufWriter::new(file);

        let mut problem = DummyProblem::new(1);
        let timestepper = RungeKutta44 {};

        problem.data[0] = 1.0;
        writeln!(buf_writer, "{} {} {}", problem.time, problem.data[0], problem.exact_solution()).unwrap();

        for _ in 0..N_STEP {
            timestepper.step(&mut problem, DT);
            writeln!(buf_writer, "{} {} {}", problem.time, problem.data[0], problem.exact_solution()).unwrap();
        }
    }

    #[test]
    fn ssp_runge_kutta_3_3() {
        let file = fs::File::create("ssp_runge_kutta_3_3.csv").unwrap();
        let mut buf_writer = BufWriter::new(file);

        let mut problem = DummyProblem::new(1);
        let timestepper = SspRungeKutta33 {};

        problem.data[0] = 1.0;
        writeln!(buf_writer, "{} {} {}", problem.time, problem.data[0], problem.exact_solution()).unwrap();

        for _ in 0..N_STEP {
            timestepper.step(&mut problem, DT);
            writeln!(buf_writer, "{} {} {}", problem.time, problem.data[0], problem.exact_solution()).unwrap();
        }
    }
}
