use ks_rs::adr::one_dim::{
    DomainParams,
    Problem1D,
    DiffusionOnly,
    ProblemFunctions,
};
use ks_rs::timestepping::{
    ExplicitTimeStepper,
    EulerForward,
    //SspRungeKutta3,
    //RungeKutta4
};
use std::fs;
use std::io::BufWriter;
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn set_initial_conditions<F>(problem: &mut Problem1D<F>)
    where F: ProblemFunctions
{
    use ks_rs::adr::one_dim::{Variable, Cell};

    for var in 0..problem.n_variable {
        for cell in 1..=problem.domain.n_cell {
            let cell = Cell(cell);
            let x = problem.x(cell);
            *problem.var_mut(Variable(var), cell) = x / problem.domain.width;
        }
    }

    //let central_cell = (problem.domain.n_cell - 1) / 2 + 1;
    //for var in 0..problem.n_variable {
        //for cell in 1..=problem.domain.n_cell {
            //*problem.var_mut(Variable(var), Cell(cell)) = if cell == central_cell {
                //1.0
            //} else {
                //0.0
            //};
        //}
    //}
}

fn main() -> Result<()> {
    let domain = DomainParams { n_cell: 5, width: 1.0 };
    let mut problem = Problem1D::new(1, domain, DiffusionOnly);

    set_initial_conditions(&mut problem);

    let dir = String::from("res_test");
    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let euler_forward = EulerForward::new();
    //let ssp_rk3 = SspRungeKutta3::new();

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let mut buf_writer = BufWriter::new(file);
    problem.output(&mut buf_writer)?;

    let mut i = 1;
    let dt = 0.1;

    while problem.time < 1.0 {
        euler_forward.step(&mut problem, dt);
        //ssp_rk3.step(&mut problem, dt);
        let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", i)))?;
        let mut buf_writer = BufWriter::new(file);
        problem.output(&mut buf_writer)?;
        i += 1;
    }

    Ok(())
}
