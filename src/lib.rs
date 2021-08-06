pub mod keller_segel;
mod stencil;
pub mod utilities;
pub mod timestepping;
pub mod adr;
pub mod kde;
pub mod models;

extern crate ndarray;

use libc::c_double;

// Pull in the erf() function from libm, as Rust's `std` doesn't have it
#[link(name = "m")]
extern "C" {
    fn erf(x: c_double) -> c_double;
}
