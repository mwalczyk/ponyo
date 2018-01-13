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

static PRELUDE: &'static str = "Ponyo - a 2D, semi-Lagrangian fluid solver";
const ITERATIONS: usize = 10;

fn main() {
    println!("{}", PRELUDE);

    // Create solver
    let mut solver = FluidSolver::new(512, 512);
    solver.init();

    // Update solver
    for i in 0..ITERATIONS {
        //solver.update();
        //solver.to_image(&format!("iter_{}.png", i))
    }
}
