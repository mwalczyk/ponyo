extern crate image;

mod helpers;
mod fluid_quantity;
mod fluid_solver;

use helpers::Dimension;
use fluid_quantity::FluidQuantity;
use fluid_solver::FluidSolver;

// References:
// https://pdfs.semanticscholar.org/9d47/1060d6c48308abcc98dbed850a39dbfea683.pdf
// https://github.com/ethanjli/dye-transport-simulation
// https://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid-EulerParticle.pdf

// TODO:
// 1. Run `incremental-fluids` repo to compare results
// 2. Fix advection routine so that a separate particle
//    is traced for each component of the velocity field
// 3. Based on (2) and Bridson, revise the interpolation
//    routine: is trilinear interpolation needed?
// 4. Review and understand the Gauss-Seidel algorithm
// 5. Based on (3), revise the projection routine, as
//    needed
// 6. Get rid of Vector struct, if deemed unnecessary
// 7. Fix density parameter: shouldn't this be 1000 kg/m^3?
// 8. Implement inflow method and/or several initialization
//    options
// 9. Port to the GPU using compute shaders or otherwise

static PRELUDE: &'static str = "Ponyo - a 2D, semi-Lagrangian fluid solver";
const ITERATIONS: usize = 100;

fn main() {
    println!("{}", PRELUDE);

    // Create and initialize a new solver.
    let mut solver = FluidSolver::new(512, 512);
    solver.init();

    // Set up some simulation variables.
    let mut total_t = 0.0;
    let delta_t = 0.005;
    let mut iter = 0;

    // Update the solver: run the simulation for 8.0 seconds.
    while total_t < 8.0 {

        // Take 4 separate sub-steps during each iteration.
        for _ in 0..4 {
            solver.update(delta_t);
            total_t += delta_t;
        }

        // Save the current frame to disk.
        solver.to_image(&format!("images/iter_{}.png", iter));
        iter += 1;

        println!("Total run time: {} seconds", total_t);
    }
}
