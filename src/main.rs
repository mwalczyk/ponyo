mod helpers;
mod fluid_quantity;
mod fluid_solver;

use helpers::Dimension;
use fluid_quantity::FluidQuantity;
use fluid_solver::FluidSolver;

// References:
// https://pdfs.semanticscholar.org/9d47/1060d6c48308abcc98dbed850a39dbfea683.pdf

static INTRO_TEXT: &'static str = "Ponyo - a 2D, semi-Lagrangian fluid solver";
const NUM_ITERATIONS: usize = 10;

fn main() {
    println!("{}", INTRO_TEXT);

    // Create solver
    let solver_dims = Dimension { nx: 512, ny: 512 };
    let mut solver = FluidSolver::new(solver_dims);

    // Update solver
    for i in 0..NUM_ITERATIONS {
        // ...
    }
}
