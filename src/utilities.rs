use serde::Serialize;
use std::fs;
use std::io::BufWriter;
use std::borrow::Borrow;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// Extension trait providing functions to find the minimum and maximum values from an iterator
/// over f64
pub(crate) trait IterMinMax
{
    fn min_all(self) -> f64;
    fn max_all(self) -> f64;
}

/// Implement `IterMinMax` for any iterator over types that can be borrowed as f64.
impl<F, T> IterMinMax for T
where
    F: Borrow<f64>,
    T: Iterator<Item=F>
{
    /// The minimum of the items in the iterator. Base case is `INFINITY` and this is returned if
    /// no item is less than this.
    fn min_all(self) -> f64 {
        self.fold(f64::INFINITY, |x_1, x_2| x_1.min(*x_2.borrow()))
    }

    /// The maximum of the items in the iterator. Base case is 0 and this is returned if no item is
    /// greater than this.
    fn max_all(self) -> f64 {
        self.fold(0.0, |x_1, x_2| x_1.max(*x_2.borrow()))
    }
}

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
