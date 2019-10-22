use serde::Serialize;
use std::fs;
use std::io::BufWriter;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

// Starting with NaN works because f64::max ignores the NaN and picks the other value
fn max(v: &[f64]) -> f64 {
    v.iter().copied().fold(std::f64::NAN, f64::max)
}

fn min(v: &[f64]) -> f64 {
    v.iter().copied().fold(std::f64::NAN, f64::min)
}

fn all_positive(v: &[f64]) -> bool {
    v.iter().all(|&x| x >= 0.0)
}

fn all_negative(v: &[f64]) -> bool {
    v.iter().all(|&x| x <= 0.0)
}

pub fn minmod(v: &[f64]) -> f64 {
    if all_positive(v) {
        min(v)
    } else if all_negative(v) {
        max(v)
    } else {
        0.0
    }
}

pub fn dump_default_to_json_file<T>(filename: &str) -> Result<()>
    where T: Default + Serialize
{
    let file = fs::File::create(filename)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &T::default())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_test() {
        let v = vec![0.0, 0.1, 100.4, -100.4];
        assert_eq!(max(&v), 100.4);

        let v = vec![std::f64::NAN, 0.0, 0.1, 100.4, -100.4];
        assert_eq!(max(&v), 100.4);
    }

    #[test]
    fn min_test() {
        let v = vec![0.0, 0.1, 100.4, -100.4];
        assert_eq!(min(&v), -100.4);

        let v = vec![std::f64::NAN, 0.0, 0.1, 100.4, -100.4];
        assert_eq!(min(&v), -100.4);
    }

    #[test]
    fn all_positive_test() {
        let v = vec![0.0, 0.1, 100.4, -100.4];
        assert_eq!(all_positive(&v), false);

        let v = vec![0.0, 0.1, 100.4, 2100.4];
        assert_eq!(all_positive(&v), true);
    }

    #[test]
    fn all_negative_test() {
        let v = vec![0.0, 0.1, 100.4, -100.4];
        assert_eq!(all_negative(&v), false);

        let v = vec![-0.01, -0.1, -100.4, -2100.4];
        assert_eq!(all_negative(&v), true);
    }
}
