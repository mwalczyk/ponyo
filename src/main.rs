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
// Gauss-Siedel: https://www.youtube.com/watch?v=ajJD0Df5CsY

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
const ITERATIONS: usize = 6000;

fn main() {
    println!("{}", PRELUDE);

    // Create and initialize a new solver.
    let mut solver = FluidSolver::new(128, 128);
    solver.init();

    for i in 0..ITERATIONS {
        solver.update();
        solver.to_image(&format!("images/iter_{}.png", i));

        if i % 500 == 0 && i != ITERATIONS {
            solver.init();
        }
        println!("Completed iteration: {}", i);
    }
}
