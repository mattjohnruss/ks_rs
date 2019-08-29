use ks_rs::adr::one_dim::{
    DomainParams,
    Problem1D,
    ProblemFunctions,
    Variable,
    Cell,
    BoundaryCondition,
    Face,
};
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    //EulerForward,
    SspRungeKutta33,
};
use std::fs;
use std::io::{Write, BufWriter};
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

struct DiffusionDirichlet;

impl DiffusionDirichlet {
    const C: Variable = Variable(0);
    const N_VARIABLE: usize = 1;
}

impl ProblemFunctions for DiffusionDirichlet {
    fn left_bc(&self, _problem: &Problem1D<Self>, _var: Variable) -> BoundaryCondition {
        BoundaryCondition::Dirichlet(1.0)
    }

    fn right_bc(&self, _problem: &Problem1D<Self>, _var: Variable) -> BoundaryCondition {
        BoundaryCondition::Dirichlet(0.0)
    }
}

fn set_initial_conditions<F>(problem: &mut Problem1D<F>)
    where F: ProblemFunctions
{
    //*problem.var_mut(DiffusionDirichlet::C, Cell(1)) = 1.0;

    for cell in problem.interior_cells() {
        *problem.var_mut(DiffusionDirichlet::C, cell) = (1.0 - problem.x(cell)).powi(2);
    }
}

fn trace_header<W: Write>(buffer: &mut W) -> std::io::Result<()> {
    buffer.write_all(b"t c_left c_right\n")
}

fn trace<F, W>(problem: &Problem1D<F>, buffer: &mut W) -> std::io::Result<()>
    where F: ProblemFunctions,
          W: Write
{
    let left_val = problem.var_point_value_at_face(DiffusionDirichlet::C, Cell(1), Face::West);
    let right_val = problem.var_point_value_at_face(DiffusionDirichlet::C, Cell(problem.domain.n_cell), Face::East);
    buffer.write_all(format!("{} {} {}\n", problem.time, left_val, right_val).as_bytes())
}

fn main() -> Result<()> {
    let domain = DomainParams { n_cell: 11, width: 1.0 };

    // Use the parameters from the Painter and Hillen paper
    let diffusion_dirichlet = DiffusionDirichlet {};

    let mut problem = Problem1D::new(DiffusionDirichlet::N_VARIABLE, domain, diffusion_dirichlet);
    problem.variable_names = vec![String::from("c")];

    set_initial_conditions(&mut problem);

    let dir = String::from("res_test");
    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    //let mut euler_forward = EulerForward::new(problem.n_dof);
    let mut ssp_rk33 = SspRungeKutta33::new(problem.n_dof);

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.update_ghost_cells();
    problem.output(&mut buf_writer)?;
    //problem.output_cell_averages(&mut buf_writer)?;
    buf_writer.flush()?;

    let trace_file = fs::File::create(dir_path.join("trace.csv"))?;
    let mut trace_writer = BufWriter::new(trace_file);
    trace_header(&mut trace_writer)?;
    trace_writer.flush()?;

    let output_interval = 100;
    let mut i = 1;
    let dt = 1e-5;
    let t_max = 1.0;

    while problem.time < t_max {
        //euler_forward.step(&mut problem, dt);
        ssp_rk33.step(&mut problem, dt);

        if i % output_interval == 0 {
            let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let mut buf_writer = BufWriter::new(file);
            println!("Outputting at time = {}, i = {}", problem.time, i);
            problem.output(&mut buf_writer)?;
            //problem.output_cell_averages(&mut buf_writer)?;
            trace(&problem, &mut trace_writer)?;
        }
        i += 1;
    }

    Ok(())
}
