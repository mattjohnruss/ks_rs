pub mod adr;
pub mod kde;
pub mod keller_segel;
pub mod models;
pub mod steady_state;
mod stencil;
pub mod timestepping;
pub mod utilities;

use libc::c_double;

// Pull in the erf() function from libm, as Rust's `std` doesn't have it
#[link(name = "m")]
extern "C" {
    fn erf(x: c_double) -> c_double;
}
