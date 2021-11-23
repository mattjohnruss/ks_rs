use ks_rs::adr::one_dim::{
    BoundaryCondition, Cell, DomainParams, Face, Problem1D, ProblemFunctions, Variable,
};
use ks_rs::steady_state::SteadyStateDetector;
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    //EulerForward,
    SspRungeKutta33,
};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

struct DiffusionDirichlet;

impl DiffusionDirichlet {
    const C: Variable = Variable(0);
    const N_VARIABLE: usize = 1;
}

impl ProblemFunctions for DiffusionDirichlet {
    fn velocity_p_at_midpoint(
        &self,
        _problem: &Problem1D<Self>,
        _var: Variable,
        _cell: Cell,
    ) -> f64 {
        1.0
    }

    fn left_bc(&self, _problem: &Problem1D<Self>, _var: Variable) -> BoundaryCondition {
        BoundaryCondition::Dirichlet(0.0)
    }

    fn right_bc(&self, _problem: &Problem1D<Self>, _var: Variable) -> BoundaryCondition {
        BoundaryCondition::Dirichlet(1.0)
    }

    fn diffusivity(&self, _problem: &Problem1D<Self>, _var: Variable, _cell: Cell) -> f64 {
        2.0
    }

    fn reactions(&self, problem: &Problem1D<Self>, var: Variable, cell: Cell) -> f64 {
        100.0 - 100.0 * problem.var(var, cell)
    }
}

fn set_initial_conditions<F>(problem: &mut Problem1D<F>)
where
    F: ProblemFunctions,
{
    //*problem.var_mut(DiffusionDirichlet::C, Cell(1)) = 1.0;

    for cell in problem.interior_cells() {
        *problem.var_mut(DiffusionDirichlet::C, cell) = 1.0 - (1.0 - problem.x(cell)).powi(2);
    }

    problem.update_ghost_cells();
}

fn trace_header<W: Write>(buffer: &mut W) -> std::io::Result<()> {
    buffer.write_all(b"t c_left c_right\n")
}

fn trace<F, W>(problem: &Problem1D<F>, buffer: &mut W) -> std::io::Result<()>
where
    F: ProblemFunctions,
    W: Write,
{
    let left_val = problem.var_point_value_at_face(DiffusionDirichlet::C, Cell(1), Face::West);
    let right_val = problem.var_point_value_at_face(
        DiffusionDirichlet::C,
        Cell(problem.domain.n_cell),
        Face::East,
    );
    buffer.write_all(format!("{} {} {}\n", problem.time, left_val, right_val).as_bytes())
}

fn main() -> Result<()> {
    let domain = DomainParams {
        n_cell: 101,
        width: 1.0,
    };

    let diffusion_dirichlet = DiffusionDirichlet {};

    let mut problem = Problem1D::new(DiffusionDirichlet::N_VARIABLE, domain, diffusion_dirichlet);
    problem.variable_names = vec![String::from("c")];

    set_initial_conditions(&mut problem);

    let dir = ".";
    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    //let mut euler_forward = EulerForward::new(problem.n_dof);
    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let mut ssd = SteadyStateDetector::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let cell_averages_file =
        fs::File::create(dir_path.join(format!("output_averages_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    let mut cell_averages_buf_writer = BufWriter::new(cell_averages_file);
    problem.output(&mut buf_writer)?;
    problem.output_cell_averages(&mut cell_averages_buf_writer)?;
    buf_writer.flush()?;
    cell_averages_buf_writer.flush()?;

    let trace_file = fs::File::create(dir_path.join("trace.csv"))?;
    let mut trace_writer = BufWriter::new(trace_file);
    trace_header(&mut trace_writer)?;
    trace(&problem, &mut trace_writer)?;
    trace_writer.flush()?;

    let output_interval = 1000;
    let mut i = 1;
    let dt = 1.0e-6;
    let t_max = 1.0;
    let ssd_threshold = 1.0e-6;
    let mut reached_steady_state = false;

    while problem.time < t_max {
        if !reached_steady_state && ssd.is_steady_state(&problem, ssd_threshold) {
            println!("Steady state reached at t = {} (within threshold {:e})", problem.time, ssd_threshold);
            reached_steady_state = true;
        }

        //euler_forward.step(&mut problem, dt);
        ssp_rk33.step(&mut problem, dt);

        if i % output_interval == 0 {
            let file =
                fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let cell_averages_file = fs::File::create(
                dir_path.join(format!("output_averages_{:05}.csv", i / output_interval)),
            )?;
            let mut buf_writer = BufWriter::new(file);
            let mut cell_averages_buf_writer = BufWriter::new(cell_averages_file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;
            problem.output_cell_averages(&mut cell_averages_buf_writer)?;
            trace(&problem, &mut trace_writer)?;
        }
        i += 1;
    }

    Ok(())
}
