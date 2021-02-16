pub trait Kernel {
    fn value(&self, x: f64) -> f64;
    fn integral(&self, _x: f64) -> f64 {
        todo!("Implement numerical integration as default")
    }
}

pub struct Rectangular;

impl Kernel for Rectangular {
    fn value(&self, x: f64) -> f64 {
        if x.abs() <= 1.0 {
            0.5
        } else {
            0.0
        }
    }

    fn integral(&self, x: f64) -> f64 {
        if x < -1.0 {
            0.0
        } else if x.abs() <= 1.0 {
            0.5 * (1.0 + x)
        } else {
            1.0
        }
    }
}

pub struct Triangular;

impl Kernel for Triangular {
    fn value(&self, x: f64) -> f64 {
        if x.abs() <= 1.0 {
            1.0 - x.abs()
        } else {
            0.0
        }
    }

    fn integral(&self, x: f64) -> f64 {
        if x < -1.0 {
            0.0
        } else if x >= -1.0 && x <= 0.0 {
            0.5 * (1.0 + x).powi(2)
        } else if x > 0.0 && x <= 1.0 {
            0.5 + x - 0.5 * x.powi(2)
        } else {
            1.0
        }
    }
}

pub struct Gaussian;

impl Kernel for Gaussian {
    fn value(&self, x: f64) -> f64 {
        1.0 / (2.0 * std::f64::consts::PI).sqrt() * (-0.5 * x * x).exp()
    }

    fn integral(&self, x: f64) -> f64 {
        0.5 * unsafe { crate::erf(x / 2.0_f64.sqrt()) }
    }
}

pub struct DensityEsimator<K: Kernel> {
    k: K,
    data: Vec<f64>,
    h: f64,
}

impl<K: Kernel> DensityEsimator<K> {
    pub fn new(k: K, data: &[f64], h: f64) -> Self {
        DensityEsimator {
            k,
            data: data.to_vec(),
            h,
        }
    }

    pub fn value(&self, x: f64) -> f64 {
        let mut sum = 0.0;
        let h_recip = 1.0 / self.h;
        let n_recip = 1.0 / self.data.len() as f64;

        for &d in &self.data {
            sum += h_recip * n_recip * self.k.value(h_recip * (x - d));
        }
        sum
    }

    pub fn average(&self, x_min: f64, x_max: f64) -> f64 {
        let mut sum = 0.0;
        let h_recip = 1.0 / self.h;
        let n_recip = 1.0 / self.data.len() as f64;
        let dx_recip = 1.0 / (x_max - x_min);

        for &d in &self.data {
            // A factor of `h` appears from changing integration variables and this
            // cancels with the factor of `1/h`.
            sum += dx_recip * n_recip *
                (self.k.integral(h_recip * (x_max - d)) - self.k.integral(h_recip * (x_min - d)));
        }
        sum
    }
}
