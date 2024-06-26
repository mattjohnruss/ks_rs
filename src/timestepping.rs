use ndarray::prelude::*;
use std::ops::MulAssign;

/// Represent an object whose data can be timestepped with an explicit scheme
pub trait ExplicitTimeSteppable {
    /// Get the current time
    fn time(&self) -> f64;

    /// Get a mut reference to the current time
    fn time_mut(&mut self) -> &mut f64;

    /// Get an `ArrayView1` into the current state of the degrees of freedom
    fn dofs(&self) -> ArrayView1<f64>;

    /// Get an `ArrayViewMut1` into the current state of the degrees of freedom
    fn dofs_mut(&mut self) -> ArrayViewMut1<f64>;

    /// Set the object's degrees of freedom from `dofs`
    fn set_dofs(&mut self, dofs: ArrayView1<f64>) {
        self.dofs_mut().assign(&dofs);
    }

    /// Perform `dofs += factor * increment`
    fn scaled_add_dofs(&mut self, factor: f64, increment: ArrayView1<f64>) {
        self.dofs_mut().scaled_add(factor, &increment);
    }

    /// Perform `dofs *= factor`
    fn scale_dofs(&mut self, factor: f64) {
        self.dofs_mut().mul_assign(factor);
    }

    /// Fills the given buffer with the current right-hand side values
    fn rhs(&self, buffer: ArrayViewMut1<f64>);

    /// Actions performed before each explicit timestep
    fn actions_before_explicit_timestep(&mut self) {}

    /// Actions performed before each stage of the scheme
    fn actions_before_explicit_stage(&mut self) {}

    /// Actions performed after each stage of the scheme
    fn actions_after_explicit_stage(&mut self) {}

    /// Actions performed after each explicit timestep
    fn actions_after_explicit_timestep(&mut self) {}
}

/// Trait for explicit timesteppers
pub trait ExplicitTimeStepper {
    /// Take a timestep of length `dt` on the timesteppable object `obj`
    fn step<T: ExplicitTimeSteppable>(&mut self, obj: &mut T, dt: f64);
}

/// Euler forward timestepper
#[derive(Clone)]
pub struct EulerForward {
    rhs_buffer: Array1<f64>,
}

impl EulerForward {
    /// Create a new forward Euler timestepper
    pub fn new(n_dof: usize) -> Self {
        EulerForward {
            rhs_buffer: Array::zeros(n_dof),
        }
    }
}

impl ExplicitTimeStepper for EulerForward {
    fn step<T: ExplicitTimeSteppable>(&mut self, obj: &mut T, dt: f64) {
        obj.actions_before_explicit_timestep();
        obj.actions_before_explicit_stage();

        obj.rhs(self.rhs_buffer.view_mut());
        obj.scaled_add_dofs(dt, self.rhs_buffer.view());
        *obj.time_mut() += dt;

        obj.actions_after_explicit_stage();
        obj.actions_after_explicit_timestep();
    }
}

/// Fourth-order, four stage Runge-Kutta timestepper
#[derive(Clone)]
pub struct RungeKutta44 {
    rhs_buffer: Array1<f64>,
}

impl RungeKutta44 {
    /// Create a new Runge-Kutta timestepper
    pub fn new(n_dof: usize) -> Self {
        RungeKutta44 {
            rhs_buffer: Array::zeros(n_dof),
        }
    }
}

impl ExplicitTimeStepper for RungeKutta44 {
    fn step<T: ExplicitTimeSteppable>(&mut self, obj: &mut T, dt: f64) {
        obj.actions_before_explicit_timestep();

        obj.actions_before_explicit_stage();

        let dofs = obj.dofs().to_owned();

        // Stage 1
        obj.rhs(self.rhs_buffer.view_mut());
        let k_1 = self.rhs_buffer.to_owned();
        obj.scaled_add_dofs(0.5 * dt, k_1.view());
        *obj.time_mut() += 0.5 * dt;
        obj.actions_after_explicit_stage();

        obj.actions_before_explicit_stage();

        // Stage 2
        obj.rhs(self.rhs_buffer.view_mut());
        let k_2 = self.rhs_buffer.to_owned();
        obj.set_dofs(dofs.view());
        obj.scaled_add_dofs(0.5 * dt, k_2.view());
        obj.actions_after_explicit_stage();

        obj.actions_before_explicit_stage();

        // Stage 3
        obj.rhs(self.rhs_buffer.view_mut());
        let k_3 = self.rhs_buffer.to_owned();
        obj.set_dofs(dofs.view());
        obj.scaled_add_dofs(dt, k_3.view());
        *obj.time_mut() += 0.5 * dt;
        obj.actions_after_explicit_stage();

        obj.actions_before_explicit_stage();

        // Stage 4
        obj.rhs(self.rhs_buffer.view_mut());
        let k_4 = self.rhs_buffer.to_owned();
        obj.set_dofs(dofs.view());

        obj.scaled_add_dofs(dt / 6.0, k_1.view());
        obj.scaled_add_dofs(dt / 3.0, k_2.view());
        obj.scaled_add_dofs(dt / 3.0, k_3.view());
        obj.scaled_add_dofs(dt / 6.0, k_4.view());

        obj.actions_after_explicit_stage();

        obj.actions_after_explicit_timestep();
    }
}

/// Third-order, three stage strong stability-preserving Runge-Kutta timestepper
#[derive(Clone)]
pub struct SspRungeKutta33 {
    rhs_buffer: Array1<f64>,
    v_0: Array1<f64>,
    v_2: Array1<f64>,
}

impl SspRungeKutta33 {
    /// Create a new SSP-RK timestepper
    pub fn new(n_dof: usize) -> Self {
        SspRungeKutta33 {
            rhs_buffer: Array::zeros(n_dof),
            v_0: Array::zeros(n_dof),
            v_2: Array::zeros(n_dof),
        }
    }
}

impl ExplicitTimeStepper for SspRungeKutta33 {
    // NOTE from Hesthaven & Warburton (2008), p.158
    // NOTE double check the placement of actions_before/after_explicit_stage()
    fn step<T: ExplicitTimeSteppable>(&mut self, obj: &mut T, dt: f64) {
        obj.actions_before_explicit_timestep();

        // get the current dofs
        self.v_0.assign(&obj.dofs());

        // work out f(v_0, t)
        obj.actions_before_explicit_stage();
        obj.rhs(self.rhs_buffer.view_mut());
        // multiply by dt and add to current dofs (an Euler step)
        obj.scaled_add_dofs(dt, self.rhs_buffer.view());
        // don't actually need to save a copy of v_1
        //let v_1 = obj.dofs().to_vec();

        // increment time by dt in prep for the next stage
        *obj.time_mut() += dt;
        obj.actions_after_explicit_stage();

        // work out f(v_1, t + dt)
        obj.actions_before_explicit_stage();
        obj.rhs(self.rhs_buffer.view_mut());
        // pre: dofs == v_1; post: dofs == v_1 + dt*rhs_2
        obj.scaled_add_dofs(dt, self.rhs_buffer.view());
        // pre: dofs == v_1 + dt*rhs_2; post: dofs == 3.0*v_0 + v_1 + dt*rhs_2
        obj.scaled_add_dofs(3.0, self.v_0.view());
        // pre: dofs == 3.0*v_0 + v_1 + dt*rhs_2; post: dofs == 0.25*(3.0*v_0 + v_1 + dt*rhs_2)
        obj.scale_dofs(0.25);
        self.v_2.assign(&obj.dofs());

        // this is strange, but apparently the next rhs is evaluated back half a step...
        *obj.time_mut() -= 0.5 * dt;
        obj.actions_after_explicit_stage();

        // work out f(v_2, t + 0.5*dt)
        obj.actions_before_explicit_stage();
        obj.rhs(self.rhs_buffer.view_mut());
        // pre: dofs == v_2; post: dofs == 2.0*v_2
        obj.scaled_add_dofs(1.0, self.v_2.view());
        // pre: dofs == 2.0*v_2; post: dofs == v_0 + 2.0*v_2
        obj.scaled_add_dofs(1.0, self.v_0.view());
        // pre: dofs == v_0 + 2.0*v_2; post: dofs == v_0 + 2.0*v_2 + 2.0*dt*rhs_3
        obj.scaled_add_dofs(2.0 * dt, self.rhs_buffer.view());
        // pre: dofs == v_0 + 2.0*v_2 2.0*dt*rhs_3; post: dofs == 1/3(v_0 + 2.0*v_2 + 2.0*dt*rhs_3)
        obj.scale_dofs(1.0 / 3.0);

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
        fn time(&self) -> f64 {
            self.time
        }

        fn time_mut(&mut self) -> &mut f64 {
            &mut self.time
        }

        fn dofs(&self) -> ArrayView1<f64> {
            self.data.view()
        }

        fn dofs_mut(&mut self) -> ArrayViewMut1<f64> {
            self.data.view_mut()
        }

        fn rhs(&self, mut buffer: ArrayViewMut1<f64>) {
            for i in 0..self.n_dof {
                buffer[i] = self.calculate_rhs(i);
            }
        }
    }

    const N_STEP: usize = 100;
    const DT: f64 = 1.0 / N_STEP as f64;

    #[test]
    #[ignore]
    fn euler_forward() {
        let file = fs::File::create("euler_forward.csv").unwrap();
        let mut buf_writer = BufWriter::new(file);

        let mut problem = DummyProblem::new(1);
        let mut timestepper = EulerForward::new(1);

        problem.data[0] = 1.0;
        writeln!(
            buf_writer,
            "{} {} {}",
            problem.time,
            problem.data[0],
            problem.exact_solution()
        )
        .unwrap();

        for _ in 0..N_STEP {
            timestepper.step(&mut problem, DT);
            writeln!(
                buf_writer,
                "{} {} {}",
                problem.time,
                problem.data[0],
                problem.exact_solution()
            )
            .unwrap();
        }
    }

    #[test]
    #[ignore]
    fn runge_kutta_4_4() {
        let file = fs::File::create("runge_kutta_4_4.csv").unwrap();
        let mut buf_writer = BufWriter::new(file);

        let mut problem = DummyProblem::new(1);
        let mut timestepper = RungeKutta44::new(1);

        problem.data[0] = 1.0;
        writeln!(
            buf_writer,
            "{} {} {}",
            problem.time,
            problem.data[0],
            problem.exact_solution()
        )
        .unwrap();

        for _ in 0..N_STEP {
            timestepper.step(&mut problem, DT);
            writeln!(
                buf_writer,
                "{} {} {}",
                problem.time,
                problem.data[0],
                problem.exact_solution()
            )
            .unwrap();
        }
    }

    #[test]
    #[ignore]
    fn ssp_runge_kutta_3_3() {
        let file = fs::File::create("ssp_runge_kutta_3_3.csv").unwrap();
        let mut buf_writer = BufWriter::new(file);

        let mut problem = DummyProblem::new(1);
        let mut timestepper = SspRungeKutta33::new(1);

        problem.data[0] = 1.0;
        writeln!(
            buf_writer,
            "{} {} {}",
            problem.time,
            problem.data[0],
            problem.exact_solution()
        )
        .unwrap();

        for _ in 0..N_STEP {
            timestepper.step(&mut problem, DT);
            writeln!(
                buf_writer,
                "{} {} {}",
                problem.time,
                problem.data[0],
                problem.exact_solution()
            )
            .unwrap();
        }
    }
}
