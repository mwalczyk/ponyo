extern crate image;
extern crate pbr;
extern crate clap;

mod helpers;
mod fluid_quantity;
mod fluid_solver;
mod render;

use std::env;
use std::str::FromStr;

use fluid_quantity::FluidQuantity;
use fluid_solver::FluidSolver;

use pbr::ProgressBar;
use clap::{Arg, App, SubCommand};

// TODO:
// https://en.wikipedia.org/wiki/Stencil_code
// http://developer.download.nvidia.com/books/HTML/gpugems/gpugems_ch38.html

const ADD_FLUID_EVERY: usize = 500;

fn main() {
    let sp = render::shader_program::ShaderProgram { id: 32 };
    sp.test();

    // Parse command-line arguments.
    let matches = App::new("Ponyo")
        .version("1.0")
        .author("Mike Walczyk <mwalczyk2@gmail.com>")
        .about("A semi-Lagrangian fluid solver")
        .arg(Arg::with_name("output_directory")
            .short("o")
            .long("output_directory")
            .value_name("DIRECTORY")
            .help("Sets the directory where images will be rendered to")
            .takes_value(true)
        )
        .arg(Arg::with_name("frame_count")
            .short("f")
            .long("frame_count")
            .value_name("COUNT")
            .help("Sets the number of frames that will be rendered")
            .takes_value(true)
        )
        .arg(Arg::with_name("resolution")
            .short("r")
            .long("resolution")
            .value_name("WIDTHxHEIGHT")
            .help("Sets the resolution of the simulation, in pixels")
            .takes_value(true)
            .use_delimiter(true)
            .value_delimiter("x")
        ).get_matches();

    let v = Arg::with_name("v");

    // Print the final program arguments.
    let output_directory = matches.value_of("output_directory").unwrap_or("images");
    println!("Setting output directory to: {}", output_directory);

    let frame_count = usize::from_str(matches.value_of("frame_count").unwrap_or("6000")).unwrap();
    println!("Setting frame count to: {}", frame_count);

    let resolution = matches.value_of("resolution").unwrap_or("128x128");
    let dims: Vec<&str> = resolution.split("x").collect();
    let w = usize::from_str(dims[0]).unwrap();
    let h = usize::from_str(dims[1]).unwrap();
    println!("Setting resolution to: {}x{}", w, h);

    // Create and initialize a new solver.
    let mut solver = FluidSolver::new(w, h);

    // Parameters for adding new fluid.
    let block_w = 30_usize;
    let block_h = 30_usize;
    let upper_left_x = w / 2 - block_w / 2;
    let upper_left_y = h / 2 - block_h / 2;
    let u = 6.0;
    let v = 0.0;

    // Initialize the progress bar.
    let mut progress = ProgressBar::new(frame_count as u64);

    // Start the simulation.
    for i in 0..frame_count {
        // Increment the progress bar.
        progress.inc();

        // Add some fluid, except on the last frame.
        if i % ADD_FLUID_EVERY == 0 && i != frame_count {
            solver.add_source(upper_left_x,
                              upper_left_y,
                              block_w,
                              block_h,
                              u,
                              v);
        }

        // Update the solver and save out frames to disk.
        solver.update();
        solver.to_image(&format!("{}/frame_{}.png", output_directory, i));
        println!("Completed iteration: {}", i);
    }
}
