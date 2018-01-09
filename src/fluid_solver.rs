use std::f64;

use helpers::{Dimension, Vector};
use fluid_quantity::{Staggered, FluidQuantity};

pub struct FluidSolver {
    dims: Dimension,
    p: FluidQuantity,   // Pressure field
    u: FluidQuantity,   // i-component of velocity field
    v: FluidQuantity,   // j-component of velocity field
    d: FluidQuantity,   // The density field of the fluid
    grid_cell_size: f64 // The size of each grid cell
}

impl FluidSolver {
    pub fn new(dims: Dimension) -> FluidSolver {
        println!("Initializing a new FluidSolver...");

        // The u and v-components of the velocity field are one unit
        // larger along their respective dimensions to form a staggered
        // MAC (marker-and-cell) grid.
        FluidSolver {
            dims,
            p: FluidQuantity::new(dims, Staggered::None),
            u: FluidQuantity::new(dims.expand(1, 0), Staggered::OffsetX),
            v: FluidQuantity::new(dims.expand(0, 1), Staggered::OffsetY),
            d: FluidQuantity::new(dims, Staggered::None),
            grid_cell_size: 1.0 / (dims.nx as f64)
        }
    }

    pub fn to_image(&self) {
        // Saves out an image derived from the solver's density
        // field.
        // TODO
    }

    fn init() {
        // Create an initial divergence-free velocity field with
        // components `u` and `v`. Initialize the pressure field
        // `p` and density field `d`.
        // TODO
    }

    fn velocity_at(&self, i: usize, j: usize) {
        // Returns the velocity at grid cell (i, j). Since we are
        // using a staggered grid, the velocity vector must be
        // reconstructed by performing linear interpolations of
        // the surrounding grid cells.
        // TODO
    }

    fn project(&mut self, delta_t: f64) {
        // Calculate and apply just the right amount of pressure to
        // make `u` divergence-free.
        // TODO
    }


    fn advect(&mut self, delta_t: f64, q: &mut FluidQuantity) {
        // Advect quantity `q` through the velocity field `u` for
        // a time interval `delta_t`. This should ONLY be called
        // with a divergence-free velocity field `u`, i.e. one that
        // meets the incompressibility constraint.
        //
        // Here, we take a semi-Lagrangian approach. We run time
        // "backwards" to find the start point of a particle that
        // ends up at each grid cell. We do this as follows:
        //
        //              x_p = x_g - delta_t * u(x_g)
        //
        // which gives us `x_p`, the previous position of the
        // hypothetical particle. Because this point will not be
        // on the grid, we simply interpolate the old value of `q`
        // from the old values on the grid around `x_p`.
        //
        // Note that we need to use the appropriate averaged
        // velocity to estimate particle trajectories.
        //
        // This semi-Lagrangian approach is unconditionally
        // stable: we can't create larger or smaller values of `q`
        // than were already present in the previous time step
        // since we are performing linear/bilinear/trilinear
        // interpolations.
        // TODO
    }

    fn determine_time_step(&self) -> f64 {
        // Find the length of the largest velocity vector in our
        // field `u`. Is this supposed to be interpolated?
        // TODO
        let mut u_max = f64::MIN_POSITIVE;

        for i in 0..self.u.dims.nx {
            for j in 0..self.u.dims.ny {
                let u_component = self.u.at(i, j);
                let v_component = self.u.at(i, j);
                let velocity = Vector::new(u_component, v_component);

                if velocity.length() > u_max {
                    u_max = velocity;
                }
            }
        }

        // The fluid should not move more than MAX_GRID_CELL_TRAVERSAL
        // grid cells per iteration.
        const MAX_GRID_CELL_TRAVERSAL: usize = 5;
        let delta_t = (MAX_GRID_CELL_TRAVERSAL as f64 * self.grid_cell_size) / u_max;

        delta_t
    }

    pub fn update(&mut self) {

        // 0. Update the hash table of marker cells (i.e. cells that
        // currently contain fluid): can be ignored initially
        // TODO

        // 1. Determine a good time step `delta_t`
        let delta_t = self.determine_time_step();

        // 2. Update the velocity field (self-advection) via backwards
        // particle trace
        // TODO

        // 3. Apply body forces (i.e. gravity): can be ignored initially
        // TODO

        // 4. Project the velocity field to obey the incompressibility condition:
        //      a. Setup the matrix A of coefficients
        //      b. Setup the vector b, which contains `div(u)` for each grid cell
        //      c. Use the Gauss-Siedel method to solve for the pressures
        //      d. Apply the pressure:
        //         u' = u - delta_t / (density * grid_cell_size) * grad(pressure)
        //                             ^^^^^^^
        //         Note that the `density` term above was a constant (1.0) before.
        //         How does this work when the density is a field?
        // TODO
    }
}
