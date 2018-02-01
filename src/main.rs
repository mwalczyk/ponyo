extern crate image;

mod helpers;
mod fluid_quantity;
mod fluid_solver;

use helpers::Dimension;
use fluid_quantity::FluidQuantity;
use fluid_solver::FluidSolver;

static PRELUDE: &'static str = "Ponyo - a 2D, semi-Lagrangian fluid solver";
const TOTAL_FRAMES: usize = 500;
const ADD_FLUID_EVERY: usize = 500;

fn main() {
    println!("{}", PRELUDE);

    // Create and initialize a new solver.
    let w = 128_usize;
    let h = 128_usize;
    let mut solver = FluidSolver::new(w, h);

    // Parameters for adding new fluid.
    let block_w = 30_usize;
    let block_h = 30_usize;
    let upper_left_x = w / 2 - block_w / 2;
    let upper_left_y = h / 2 - block_h / 2;
    let u = 6.0;
    let v = 0.0;

    // Start the simulation.
    for i in 0..TOTAL_FRAMES {
        // Add some fluid, except on the last frame.
        if i % ADD_FLUID_EVERY == 0 && i != TOTAL_FRAMES {
            solver.add_source(upper_left_x,
                              upper_left_y,
                              block_w,
                              block_h,
                              u,
                              v);
        }

        // Update the solver and save out frames to disk.
        solver.update();
        solver.to_image(&format!("images/frame_{}.png", i));

        println!("Completed iteration: {}", i);
    }
}
