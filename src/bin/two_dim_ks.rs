use ks_rs::keller_segel::two_dim::{ExactSolution, Forces, ICs, Parameters, Problem2D};
use std::convert::{TryFrom, TryInto};
use std::fmt;
use std::fs;
use std::io::BufWriter;
use std::path::Path;
use structopt::StructOpt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(StructOpt, Debug)]
#[structopt(name = "ks_rs", rename_all = "verbatim")]
struct Opt {
    /// Initial conditions
    #[structopt(
        long,
        default_value = "constant",
        possible_values = &["constant", "perturbed", "gaussian", "exact"]
    )]
    ics: String,

    /// Diffusivity
    #[structopt(short, long, default_value = "1.0")]
    diffusivity: f64,

    /// Logistic growth rate
    #[structopt(short, long, default_value = "0.0")]
    r: f64,

    /// Chemotaxis strength
    #[structopt(short, long, default_value = "1.0")]
    chi: f64,

    /// Something
    #[structopt(long, default_value = "1.0")]
    gamma_rho: f64,

    /// Something
    #[structopt(long, default_value = "1.0")]
    gamma_c: f64,

    /// Cells in the interior of the domain
    #[structopt(short, long, default_value = "10")]
    n_interior_cell_1d: usize,

    /// Domain length
    #[structopt(short, long, default_value = "1.0")]
    length: f64,

    /// Initial rho_bar value
    #[structopt(long, default_value = "1.0")]
    rho_bar_init: f64,

    /// Initial c value
    #[structopt(long, default_value = "1.0")]
    c_init: f64,

    /// Size of initial chemoattractant perturbation
    #[structopt(long, default_value = "0.01")]
    pert_size: f64,

    /// Height of Gaussian ICs
    #[structopt(long, default_value = "1000.0")]
    gaussian_height: f64,

    /// Width of Gaussian ICs
    #[structopt(long, default_value = "100.0")]
    gaussian_width: f64,

    /// Maximum time
    #[structopt(long)]
    t_max: f64,

    /// Timestep
    #[structopt(long)]
    dt: f64,

    /// Output interval (in timesteps)
    #[structopt(long, default_value = "1")]
    output_interval: usize,

    /// Output directory
    #[structopt(long, default_value = "res")]
    dir: String,

    /// Whether to solve for an exact solution
    #[structopt(long)]
    exact_solve: bool,
}

enum TryFromOptError {
    ICsError(String),
}

impl fmt::Debug for TryFromOptError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            TryFromOptError::ICsError(ics) => write!(f, "Specified ICs `{}` invalid", ics),
        }
    }
}

impl fmt::Display for TryFromOptError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

impl std::error::Error for TryFromOptError {}

impl TryFrom<Opt> for Parameters {
    type Error = TryFromOptError;

    fn try_from(opt: Opt) -> std::result::Result<Self, Self::Error> {
        let ics = match opt.ics.as_ref() {
            "constant" => ICs::Constant(opt.rho_bar_init, opt.c_init),
            "perturbed" => ICs::Perturbed(opt.rho_bar_init, opt.c_init, opt.pert_size),
            "gaussian" => ICs::Gaussian(opt.gaussian_height, opt.gaussian_width),
            "exact" => ICs::Exact,
            s => return Err(TryFromOptError::ICsError(s.into())),
        };

        let params = Parameters {
            ics,
            diffusivity: opt.diffusivity,
            r: opt.r,
            chi: opt.chi,
            gamma_rho: opt.gamma_rho,
            gamma_c: opt.gamma_c,
            n_interior_cell_1d: opt.n_interior_cell_1d,
            length: opt.length,
            dx: opt.length / opt.n_interior_cell_1d as f64,
            dy: opt.length / opt.n_interior_cell_1d as f64,
            exact_solution: None,
            forces: None,
        };

        Ok(params)
    }
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    let t_max = opt.t_max;
    let dt = opt.dt;
    let output_interval = opt.output_interval;
    let dir = opt.dir.clone();
    let exact_solve = opt.exact_solve;

    let mut problem = Problem2D::with_params(opt.try_into()?);

    if exact_solve {
        use std::f64::consts::PI;

        // Exact solution
        problem.p.exact_solution = Some(ExactSolution::new(
            |t, x, y, _p| {
                (x * (1.0 - x)).powi(2) * (y * (1.0 - y)).powi(2) * (PI * t).sin().powi(2)
            },
            |t, x, y, _p| {
                (x * (1.0 - x)).powi(2) * (y * (1.0 - y)).powi(2) * (PI * t).sin().powi(2)
            },
            //|t, _x, _y, _p| (PI * t).sin().powi(2),
            //|t, _x, _y, _p| (PI * t).sin().powi(2),
        ));

        // The forces that are required to find the above solution
        problem.p.forces = Some(Forces::new(
            |t, x, y, p| {
                let rho_t =
                    PI * (x * (1.0 - x)).powi(2) * (y * (1.0 - y)).powi(2) * (2.0 * PI * t).sin();
                let rho_xx = 2.0
                    * (1.0 + 6.0 * x * (x - 1.0))
                    * (y * (1.0 - y)).powi(2)
                    * (PI * t).sin().powi(2);
                let rho_yy = 2.0
                    * (1.0 + 6.0 * y * (y - 1.0))
                    * (x * (1.0 - x)).powi(2)
                    * (PI * t).sin().powi(2);
                let chemotaxis_x = 2.0
                    * p.chi
                    * (x * (1.0 - x)).powi(2)
                    * (3.0 + 14.0 * x * (x - 1.0))
                    * (y * (y - 1.0)).powi(4)
                    * (PI * t).sin().powi(4);
                let chemotaxis_y = 2.0
                    * p.chi
                    * (y * (y - 1.0)).powi(2)
                    * (3.0 + 14.0 * y * (y - 1.0))
                    * (x * (x - 1.0)).powi(4)
                    * (PI * t).sin().powi(4);

                rho_t + chemotaxis_x + chemotaxis_y - rho_xx - rho_yy
            },
            |t, x, y, p| {
                let rho =
                    (x * (1.0 - x)).powi(2) * (y * (1.0 - y)).powi(2) * (PI * t).sin().powi(2);
                let c = (x * (1.0 - x)).powi(2) * (y * (1.0 - y)).powi(2) * (PI * t).sin().powi(2);
                let c_t =
                    PI * (x * (1.0 - x)).powi(2) * (y * (1.0 - y)).powi(2) * (2.0 * PI * t).sin();
                let c_xx = 2.0
                    * (1.0 + 6.0 * x * (x - 1.0))
                    * (y * (1.0 - y)).powi(2)
                    * (PI * t).sin().powi(2);
                let c_yy = 2.0
                    * (1.0 + 6.0 * y * (y - 1.0))
                    * (x * (1.0 - x)).powi(2)
                    * (PI * t).sin().powi(2);

                c_t - c_xx - c_yy + p.gamma_c * c - p.gamma_rho * rho
            },
            //|t, _x, _y, _p| {
            //let rho_t = PI * (2.0 * PI * t).sin();
            //rho_t
            //},
            //|t, _x, _y, _p| {
            //let c_t = PI * (2.0 * PI * t).sin();
            //c_t
            //},
        ));
    }

    problem.set_initial_conditions();

    let dir_path = Path::new(&dir);
    fs::create_dir_all(dir_path)?;

    let file = fs::File::create(dir_path.join(format!("output_{:05}.csv", 0)))?;
    let buf_writer = BufWriter::new(file);
    problem.output(buf_writer)?;

    let mut i = 1;

    while problem.time < t_max {
        //problem.step_euler_forward(dt);
        problem.step_ssp_rk3(dt);

        if i % output_interval == 0 {
            println!(
                "Output {} at timestep {}, t = {}",
                i / output_interval,
                i,
                problem.time
            );
            let file =
                fs::File::create(dir_path.join(format!("output_{:05}.csv", i / output_interval)))?;
            let buf_writer = BufWriter::new(file);
            problem.output(buf_writer)?;
        }

        i += 1;
    }

    Ok(())
}
