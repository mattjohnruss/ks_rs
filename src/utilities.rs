use serde::Serialize;
use std::fs;
use std::io::BufWriter;
use std::borrow::Borrow;

use ndarray::prelude::*;
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::rand::{seq::SliceRandom, thread_rng};

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

// Translation of the Python https://gist.github.com/jgomezdans/4739643, which itself is from the
// MATLAB https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fuk.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F4352%2Fversions%2F1%2Fcontents%2Flhsu.m
pub fn lhsu(x_min: &[f64], x_max: &[f64], n_sample: usize) -> Array2<f64> {
    assert_eq!(x_min.len(), x_max.len());
    let n_var = x_min.len();

    let uniform_samples = Array::random((n_sample, n_var), Uniform::new(0.0_f64, 1.0));
    let mut samples = Array::zeros((n_sample, n_var));

    let mut rng = thread_rng();

    let mut idx: Vec<_> = (1..=n_sample).map(|i| i as f64).collect();

    for j in 0..n_var {
        idx.shuffle(&mut rng);
        let idx = Array::from_vec(idx.clone());

        let p = (idx - uniform_samples.index_axis(Axis(1), j)) / n_sample as f64;
        let mut sample_j = samples.index_axis_mut(Axis(1), j);
        sample_j.assign(&(x_min[j] + (x_max[j] - x_min[j]) * p));
    }
    samples
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

    #[test]
    fn lhsu_test() {
        for _ in 0..1000 {
            let x_min = &[0.3, 1.01, -3.123, 103.987];
            let x_max = &[8.09123, 1.02, 0.0, 6123.123];
            let n_sample = 1000;

            let samples = lhsu(x_min, x_max, n_sample);

            for sample in samples.axis_iter(Axis(0)) {
                assert_eq!(sample.len(), 4);

                for (s, (&s_min, &s_max)) in sample.iter().zip(x_min.iter().zip(x_max.iter())) {
                    assert!((s_min..=s_max).contains(s));
                }
            }
        }
    }
}
