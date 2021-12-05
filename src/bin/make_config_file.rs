use ks_rs::models::chemotaxis_7_var::ChemotaxisParameters;

use std::fs;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "make_config_file", rename_all = "verbatim")]
struct Opt {
    #[structopt(long = "base_config")]
    base_config_path: String,
    #[structopt(long = "full_config")]
    full_config_path: String,
    #[structopt(long)]
    j_phi_i_i_factor: f64,
    #[structopt(long)]
    m_i_factor: f64,
    #[structopt(long)]
    t_j_phi_i_lag: f64,
    #[structopt(long)]
    gamma: f64,
}

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result<()> {
    let opt = Opt::from_args();

    // Read the base parameters from file
    let mut chemotaxis_params: ChemotaxisParameters = {
        let base_config_file = fs::File::open(&opt.base_config_path)?;
        let reader = std::io::BufReader::new(base_config_file);
        serde_json::from_reader(reader)?
    };

    // Update using values from cmdline
    chemotaxis_params.j_phi_i_i_factor = opt.j_phi_i_i_factor;
    chemotaxis_params.m_i_factor = opt.m_i_factor;
    chemotaxis_params.t_j_phi_i_lag = opt.t_j_phi_i_lag;

    // Set all cleaving parameters to gamma
    chemotaxis_params.gamma_ui = opt.gamma;
    chemotaxis_params.gamma_um = opt.gamma;
    chemotaxis_params.gamma_bi = opt.gamma;
    chemotaxis_params.gamma_bm = opt.gamma;

    let full_config_file = std::io::BufWriter::new(fs::File::create(opt.full_config_path)?);
    serde_json::to_writer_pretty(full_config_file, &chemotaxis_params)?;

    Ok(())
}
