use ks_rs::adr::one_dim::{
    DomainParams,
    Problem1D,
    ProblemFunctions,
    Variable,
    Cell,
};
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    SspRungeKutta33,
};
use std::fs;
use std::io::{Write, BufWriter};
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

struct KellerSegel {
    d_rho: f64,
    d_c: f64,
    chi: f64,
    gamma_rho: f64,
    gamma_c: f64,
    r: f64,
}

impl KellerSegel {
    const RHO: Variable = Variable(0);
    const C: Variable = Variable(1);
    const N_VARIABLE: usize = 2;
}

impl ProblemFunctions for KellerSegel {
    fn diffusivity(&self, _problem: &Problem1D<Self>, var: Variable, _cell: Cell) -> f64 {
        match var {
            Self::RHO => self.d_rho,
            Self::C => self.d_c,
            _ => panic!("Invalid variable number"),
        }
    }

    fn velocity_p_at_midpoint(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var {
            Self::RHO => self.chi * problem.dvar_dx_p_at_midpoint(Self::C, cell),
            Self::C => 0.0,
            _ => panic!("Invalid variable number"),
        }
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        match var {
            Self::RHO => self.r * problem.var(Self::RHO, cell) * (1.0 - problem.var(Self::RHO, cell)),
            Self::C => self.gamma_rho * problem.var(Self::RHO, cell) - self.gamma_c * problem.var(Self::C, cell),
            _ => panic!("Invalid variable number"),
        }
    }
}

fn set_initial_conditions<F>(problem: &mut Problem1D<F>)
    where F: ProblemFunctions
{
    // ICs: uniform for rho, perturbed for c
    use rand::distributions::Uniform;
    use rand::{thread_rng, Rng};

    let pert_size = 0.01;
    let uniform = Uniform::new_inclusive(-pert_size, pert_size);

    for cell in problem.interior_cells() {
        *problem.var_mut(KellerSegel::RHO, cell) = 1.0;
        *problem.var_mut(KellerSegel::C, cell) = 1.0 + thread_rng().sample(uniform);
    }

    problem.update_ghost_cells();
}

fn main() -> Result<()> {
    let domain = DomainParams { n_cell: 501, width: 25.0 };

    // Use the parameters from the Painter and Hillen paper
    let ks = KellerSegel {
        d_rho: 0.1,
        d_c: 1.0,
        chi: 5.0,
        gamma_rho: 1.0,
        gamma_c: 1.0,
        r: 1.0,
    };

    let mut problem = Problem1D::new(KellerSegel::N_VARIABLE, domain, ks);
    problem.variable_names = vec![String::from("rho"), String::from("c")];

    set_initial_conditions(&mut problem);

    let dir = String::from("res_test");
    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;
    buf_writer.flush()?;

    let output_interval = 10;
    let mut i = 1;
    let dt = 1e-4;
    let t_max = 100.0;

    while problem.time < t_max {
        ssp_rk33.step(&mut problem, dt);

        if i % output_interval == 0 {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let mut buf_writer = BufWriter::new(file);
            problem.output(&mut buf_writer)?;
        }
        i += 1;
    }

    Ok(())
}
