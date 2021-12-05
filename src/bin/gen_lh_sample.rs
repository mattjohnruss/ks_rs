use ks_rs::utilities::lhsu;

use std::fs;
use std::io::prelude::*;
use structopt::StructOpt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(StructOpt, Debug)]
#[structopt(name = "gen_lh_sample", rename_all = "verbatim")]
struct Opt {
    #[structopt(long)]
    min: Vec<f64>,
    #[structopt(long)]
    max: Vec<f64>,
    #[structopt(long)]
    n_parameter_sample: usize,
    #[structopt(long)]
    output_file: String,
}

fn main() -> Result<()> {
    let opt = Opt::from_args();
    let n_param = opt.min.len();
    let samples = lhsu(&opt.min, &opt.max, opt.n_parameter_sample);

    let file = fs::File::create(&opt.output_file)?;
    let mut writer = std::io::BufWriter::new(file);

    for s in samples.rows() {
        for (i, p) in s.iter().enumerate() {
            write!(&mut writer, "{}", p)?;
            if i != n_param - 1 {
                write!(&mut writer, " ")?;
            }
        }
        writeln!(&mut writer)?;
    }

    Ok(())
}
