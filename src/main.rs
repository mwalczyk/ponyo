extern crate image;
extern crate pbr;

mod helpers;
mod fluid_quantity;
mod fluid_solver;
mod render;

use std::env;
use fluid_quantity::FluidQuantity;
use fluid_solver::FluidSolver;
use pbr::ProgressBar;

static IMAGE_DIR: &'static str = "images";
const TOTAL_FRAMES: usize = 6000;
const ADD_FLUID_EVERY: usize = 500;

fn help() {
    println!("Ponyo - a 2D, semi-Lagrangian fluid solver");
}

fn main() {
    help();

    let sp = render::shader_program::ShaderProgram { id: 32 };
    sp.test();

    // Parse command-line arguments:
    //
    // --output_dir -o
    // --resolution -r
    // --frame_count -f
    // --density -d
    // ...
    let args: Vec<String> = env::args().collect();
    match args.len() {
        1 => (),
        2 => (),
        _ => help()
    }

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

    // Initialize the progress bar.
    let mut progress = ProgressBar::new(TOTAL_FRAMES as u64);

    // Start the simulation.
    for i in 0..TOTAL_FRAMES {
        // Increment the progress bar.
        progress.inc();

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
        solver.to_image(&format!("{}/frame_{}.png", IMAGE_DIR, i));
        //println!("Completed iteration: {}", i);
    }
}
