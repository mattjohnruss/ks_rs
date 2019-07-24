/// Trait for finite-difference stencils
pub trait Stencil {
    /// Get the stencil
    fn stencil() -> &'static [(isize, f64)];

    /// Apply the stencil to a function `f`, centring at index `i`
    fn apply<F>(i: usize, f: F) -> f64
    where
        F: Fn(usize) -> f64,
    {
        Self::stencil().iter().fold(0.0, |acc, (k, w)| {
            let cell_offset = (i as isize + k) as usize;
            acc + w * f(cell_offset)
        })
    }
}

macro_rules! stencil_impl {
    ($s:ident, $e:expr) => {
        #[allow(dead_code)]
        pub struct $s;
        impl crate::stencil::Stencil for $s {
            fn stencil() -> &'static [(isize, f64)] {
                &$e
            }
        }
    }
}

/// First-order accurate stencils
pub mod first_order {
    stencil_impl!(Forward1, [(0, -1.0), (1, 1.0)]);
    stencil_impl!(Backward1, [(0, 1.0), (-1, -1.0)]);
}

/// Second-order accurate stencils
pub mod second_order {
    stencil_impl!(Central1, [(-1, -0.5), (1, 0.5)]);
    stencil_impl!(Forward1, [(0, -1.5), (1, 2.0), (2, -0.5)]);
    stencil_impl!(Backward1, [(0, 1.5), (-1, -2.0), (-2, 0.5)]);
    stencil_impl!(Central2, [(-1, 1.0), (0, -2.0), (1, 1.0)]);
}
