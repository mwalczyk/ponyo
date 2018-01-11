extern crate image;

use std::f64;

use helpers::{Dimension, Vector, lerp};
use fluid_quantity::{Staggered, FluidQuantity};

use image::{GenericImage, ImageBuffer};

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
        //
        // The MAC grid method discretizes space into a grid of cells.
        // Each cell has a pressure `p` defined at its center. It also
        // has a velocity with components `u` and `v`, but the
        // components are placed at the centers of 2 of the cell faces:
        // `u` on the x-min face and `v` on the y-min face.
        //
        // However, in code, we use integer indices so that:
        //      p(i, j) = P(i + 0.0, j + 0.0)
        //      u(i, j) = U(i - 0.5, j + 0.0)
        //      v(i, j) = V(i + 0.0, j - 0.5)
        FluidSolver {
            dims,
            p: FluidQuantity::new(dims, Staggered::None),
            u: FluidQuantity::new(dims.expand(1, 0), Staggered::OffsetX),
            v: FluidQuantity::new(dims.expand(0, 1), Staggered::OffsetY),
            d: FluidQuantity::new(dims, Staggered::None),
            grid_cell_size: 1.0 / (dims.nx.min(dims.ny) as f64)
        }
    }

    pub fn to_image(&self, path: &str) {
        // Saves out an image derived from the solver's density
        // field.
        let center = Vector::new(self.d.dims.nx as f64 * 0.5,
                                        self.d.dims.ny as f64 * 0.5);

        let image = ImageBuffer::from_fn(self.d.dims.nx as u32, self.d.dims.ny as u32, |x, y| {
            let position = Vector::new(x as f64, y as f64);

            if position.distance(center) < 100.0 {
                image::Luma([0u8])
            } else {
                image::Luma([255u8])
            }
        });

        image.save(path).unwrap();
    }

    fn init() {
        // Create an initial divergence-free velocity field with
        // components `u` and `v`. Initialize the pressure field
        // `p` and density field `d`.
        //
        // Maybe set a region to 1.0 density, with a vertical
        // velocity of 3.0 or so?
        // TODO
    }

    fn get_interpolated_velocity(&self, i: usize, j: usize, oi: f64, oj: f64) -> Vector {
        // This function returns the velocity vector within a particular
        // grid cell (i, j). Since we are using a staggered grid, the
        // velocity vector must be reconstructed by performing a linear
        // interpolation of the surrounding grid cells. To retrieve the
        // velocity vector at the cell center, `oi` and `oj` should be
        // 0.5.
        //
        // Note that in the text, we calculate the value of a quantity
        // that is staggered along the x-axis as follows:
        //      u(i, j) = (u(i - 0.5, j) + u(i + 0.5, j)) * 0.5
        //
        // In code, we represent u(i, j) as:
        //      U(i - 0.5, j + 0.0)
        //
        // Which leads to the calculations below.
        let u_inter = self.u.at(i, j) * (1.0 - oi) + self.u.at(i + 1, j) * oi;
        let v_inter = self.v.at(i, j) * (1.0 - oj) + self.v.at(i, j + 1) * oj;

        Vector::new(u_inter, v_inter)
    }

    fn project(&mut self, delta_t: f64) {
        // Calculate and apply just the right amount of pressure to
        // make the velocity field divergence-free.
        // TODO
    }

    fn advect(&mut self, u: &mut FluidQuantity, delta_t: f64, q: &mut FluidQuantity) {
        // TODO
    }

    fn determine_time_step(&self) -> f64 {
        // Find the length of the largest velocity vector.
        let mut velocity_norm_max = f64::NEG_INFINITY;

        for i in 0..self.u.dims.nx {
            for j in 0..self.u.dims.ny {
                // Get the velocity vector at the center of cell (i, j).
                let velocity = self.get_interpolated_velocity(i, j, 0.5, 0.5);

                if velocity.length() > velocity_norm_max {
                    velocity_norm_max = velocity.length();
                }
            }
        }

        // The fluid should not move more than MAX_GRID_CELL_TRAVERSAL
        // grid cells per iteration.
        const MAX_GRID_CELL_TRAVERSAL: usize = 5;
        let delta_t = (MAX_GRID_CELL_TRAVERSAL as f64 * self.grid_cell_size) / velocity_norm_max;

        delta_t
    }

    fn backwards_trace(&mut self, delta_t: f64) {
        // Advect quantity `q` through the velocity field for
        // a time interval `delta_t`. This should ONLY be called
        // with a divergence-free velocity field, i.e. one that
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

        // Advection cannot happen in-place: we need to create new buffers
        let mut u_next = FluidQuantity::new(self.dims.expand(1, 0), Staggered::OffsetX);
        let mut v_next = FluidQuantity::new(self.dims.expand(0, 1), Staggered::OffsetY);

        // Iterate over all grid centers
        for i in 0..self.dims.nx {
            for j in 0..self.dims.ny {

                // The starting position of the particle
                let position = Vector::new(i as f64, j as f64);

                // Get the velocity vector at the center of cell (i, j).
                let velocity = self.get_interpolated_velocity(i, j, 0.5, 0.5);

                // Trace backwards using Runge-Kutta order two (RK2) interpolation.
                let position_midpoint = position - velocity * 0.5 * delta_t;

                // Get the fractional part of the particle's position vector
                let mut oi = position_midpoint.x - position_midpoint.x.floor();
                let mut oj = position_midpoint.y - position_midpoint.y.floor();

                // The (integer) indices of the grid cell that this particle was traced from
                let i_mid = position_midpoint.x.floor() as usize;
                let j_mid = position_midpoint.y.floor() as usize;

                // Fluid quantities that are staggered must be correctly handled here
                let staggered = Staggered::None;

                match staggered {
                    Staggered::OffsetX => oi += 0.5,
                    Staggered::OffsetY => oj += 0.5,
                    _ => ()
                }

                // Clamp
                oi = oi.max(0.0).min(self.dims.nx as f64);
                oj = oj.max(0.0).min(self.dims.ny as f64);

                // Linear interpolation of the fluid quantity via the four surrounding
                // grid cells:
                //
                //             |
                //      x_01   |   x_11
                //             |
                //      ---------------
                //             |
                //      x_00   |   x_10
                //             |
                //
                // Based on this diagram:
                //      a = lerp(x_00, x_10, .. )   <-- horizontal, bottom
                //      b = lerp(x_01, x_11, .. )   <-- horizontal, top
                //      c = lerp(a, b, .. )         <-- vertical
                let a = lerp(self.d.at(i_mid, j_mid), self.d.at(i_mid + 1, j_mid), oi);
                let b = lerp(self.d.at(i_mid, j_mid + 1), self.d.at(i_mid + 1, j_mid + 1), oi);
                let c = lerp(a, b, oj);

                // Set the velocity components at cell (i, j) in the new FluidQuantity buffers.
                // TODO

                // Set the fluid quantity at cell (i, j) in the new FluidQuantity buffer.
                // TODO
            }
        }

        self.u = u_next;
        self.v = v_next;
    }

    pub fn update(&mut self) {

        // 0. Update the hash table of marker cells (i.e. cells that
        // currently contain fluid): can be ignored initially.
        // TODO

        // 1. Determine a good time step `delta_t`
        let delta_t = self.determine_time_step();

        // 2. Update the velocity field (self-advection) via backwards
        // particle trace
        self.backwards_trace(delta_t);

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
