extern crate image;

use std::f64;

use helpers::{Dimension, Vector, lerp};
use fluid_quantity::{Staggered, FluidQuantity};

use image::{GenericImage, ImageBuffer};

pub struct FluidSolver {
    // The width of the solver grid
    w: usize,

    // The height of the solver grid
    h: usize,

    // The pressure field
    p: FluidQuantity,

    // The i-component of velocity field
    u: FluidQuantity,

    // The j-component of velocity field
    v: FluidQuantity,

    // The dye field, for visualization purposes
    dye: FluidQuantity,

    // The density of the fluid
    density: f64,

    // The size of each grid cell
    grid_cell_size: f64
}

impl FluidSolver {
    pub fn new(w: usize, h: usize) -> FluidSolver {
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
        let dims = Dimension::new(w, h);

        FluidSolver {
            w,
            h,
            p: FluidQuantity::new(w, h, Staggered::None),   // Pressure field starts as all 0.0
            u: FluidQuantity::new(w + 1, h, Staggered::OffsetX),
            v: FluidQuantity::new(w, h + 1, Staggered::OffsetY),
            dye: FluidQuantity::new(w, h, Staggered::None),
            density: 1.0, // TODO: should be 1000 kg/m^3
            grid_cell_size: 1.0 / (w.min(h) as f64)
        }
    }

    pub fn to_image(&self, path: &str) {
        let image = ImageBuffer::from_fn(self.w as u32,
                                                                                  self.h as u32,
                                                                                  |x, y| {
            let color = (1.0 - self.dye.get(x as usize, y as usize)) * 255.0;

            image::Luma([color as u8])
        });

        image.save(path).unwrap();
    }

    pub fn init(&mut self) {
        // Create an initial divergence-free velocity field with
        // components `u` and `v`. Initialize the pressure field
        // `p` and dye `d`.
        let block_size = 50;
        let block_center_x = 200;
        let block_center_y = self.h / 2;
        let u_start = 3.0;

        for i in (block_center_x - block_size / 2)..(block_center_x + block_size / 2) {
            for j in (block_center_y - block_size / 2)..(block_center_y + block_size / 2) {
                self.u.set(i, j, u_start);
                self.dye.set(i, j, 1.0);
            }
        }
    }

    fn add_source(&mut self) {
        // Adds fluid density at some user-specified region.
        // TODO
    }

    fn get_interpolated_quantity(&self, q: &FluidQuantity, x: f64, y: f64) -> f64 {
        // Returns the value of quantity `q` at real-valued coordinates `x` and `y`.
        // A bilinear interpolation will be performed in order to reconstruct values
        // that do not lie on exact grid centers.
        let mut oi = x - x.floor();
        let mut oj = y - y.floor();

        // The (integer) indices of the grid cell
        let mut i = x.floor() as usize;
        let mut j = y.floor() as usize;

        // Fluid quantities that are staggered must be correctly handled here
        match q.staggered {
            Staggered::OffsetX => oi += 0.5,
            Staggered::OffsetY => oj += 0.5,
            _ => ()
        }

        // Clamp
        i = i.max(0).min(self.w - 2); // TODO: this is not correct
        j = j.max(0).min(self.h - 2); // TODO: this is not correct

        // Bilinear interpolation of the fluid quantity via the four surrounding
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
        let x_00 = q.get(i, j);
        let x_10 = q.get(i + 1, j);
        let x_01 = q.get(i, j + 1);
        let x_11 = q.get(i + 1, j + 1);

        let a = lerp(x_00, x_10, oi);
        let b = lerp(x_01, x_11, oi);
        let c = lerp(a, b, oj);

        c
    }




    fn determine_time_step(&self) -> f64 {
        // Find the length of the largest velocity vector.
        let mut velocity_norm_max = f64::NEG_INFINITY;

        for i in 0..self.w {
            for j in 0..self.h {
                // Get the velocity vector at the center of cell (i, j).
                let velocity = Vector::new(self.get_interpolated_quantity(&self.u, i as f64, j as f64),
                                                  self.get_interpolated_quantity(&self.v, i as f64, j as f64));

                if velocity.length() > velocity_norm_max {
                    velocity_norm_max = velocity.length();
                }
            }
        }

        // The fluid should not move more than MAX_GRID_CELL_TRAVERSAL
        // grid cells per iteration.
        const MAX_GRID_CELL_TRAVERSAL: usize = 1;
        let delta_t = (MAX_GRID_CELL_TRAVERSAL as f64 * self.grid_cell_size) / velocity_norm_max;

        delta_t
    }

    // self-advection:
    // self.u = self.advect(&self.u);
    // self.v = self.advect(&self.v);
    //
    // self.dye = self.advect(&self.dye);
    fn advect_quantity(&self, q: &FluidQuantity) -> FluidQuantity {
        let mut q_next = FluidQuantity::new(q.w, q.h, q.staggered);

        for i in 0..q.w {
            for j in 0..q.h {

                // The starting position of the particle is simply (i, j).
                let position = Vector::new(i as f64, j as f64);

                // Get the velocity vector at position (i, j).
                let velocity = Vector::new(self.get_interpolated_quantity(&self.u, position.x, position.y),
                                                  self.get_interpolated_quantity(&self.v, position.x, position.y));

                // Trace backwards using Runge-Kutta order two (RK2) interpolation.
                // TODO
                let position_prev = position - velocity * delta_t;

                // Set the value of the fluid quantity in the new buffer.
                q_next.set(i, j, self.get_interpolated_quantity(q, position.x, position.y));
            }
        }

        q_next
    }

    fn advect(&mut self, delta_t: f64) {
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
        let mut u_next = self.u.clone();
        let mut v_next = self.v.clone();
        let mut dye_next = self.dye.clone();

        // Iterate over all grid centers
        for i in 0..self.w {
            for j in 0..self.h {
                {
                    // The starting position of the particle
                    let position = Vector::new(i as f64, j as f64);

                    // Get the velocity vector at the center of cell (i, j).
                    let velocity = Vector::new(self.get_interpolated_quantity(&self.u, position.x, position.y),
                                               self.get_interpolated_quantity(&self.v, position.x, position.y));

                    // Trace backwards using Runge-Kutta order two (RK2) interpolation.
                    // TODO: here we need to do a separate particle trace for EACH component of the velocity
                    let position_prev = position - velocity * delta_t;

                    // Advect dye
                    dye_next.set(i, j, self.get_interpolated_quantity(&self.dye, position_prev.x, position_prev.y));
                }

                // Advect velocity, component-wise
                u_next.set(i, j, self.get_interpolated_quantity(&self.u, position_prev.x, position_prev.y));
                v_next.set(i, j, self.get_interpolated_quantity(&self.v, position_prev.x, position_prev.y));


            }
        }

        // Swap buffers
        self.u = u_next;
        self.v = v_next;
        self.dye = dye_next;
    }

    fn apply_body_forces(&mut self, delta_t: f64) {
        // TODO
    }

    fn project(&mut self, delta_t: f64, epsilon: f64, max_iters: usize) {
        // Calculate the right amount of pressure to make the velocity
        // field divergence-free.

        // Build RHS of pressure equation.
        let mut rhs = FluidQuantity::new(self.w, self.h, Staggered::None);
        let scale = 1.0 / self.grid_cell_size;
        for i in 0..self.w {
            for j in 0..self.h {
                let div_u_x = self.u.get(i + 1, j) - self.u.get(i, j);
                let div_v_y = self.v.get(i, j + 1) - self.v.get(i, j);
                let div = div_u_x + div_v_y;

                rhs.set(i, j, -scale * div);
            }
        }

        // Project
        let scale = delta_t / (self.density * self.grid_cell_size * self.grid_cell_size);
        let mut max_delta = 0.0_f64;
        for iter in 0..max_iters {

            max_delta = 0.0;

            for i in 0..self.w {
                for j in 0..self.h {
                    let mut diag = 0.0;
                    let mut off_diag = 0.0;

                    if i > 0 {
                        diag += scale;
                        off_diag -= scale * self.p.get(i - 1, j);
                    }

                    if j > 0 {
                        diag += scale;
                        off_diag -= scale * self.p.get(i, j - 1);
                    }

                    if i < (self.w - 1) {
                        diag += scale;
                        off_diag -= scale * self.p.get(i + 1, j);
                    }

                    if j < (self.h - 1) {
                        diag += scale;
                        off_diag -= scale * self.p.get(i, j + 1);
                    }

                    let p_next = (rhs.get(i, j) - off_diag) / diag;

                    // Can we exit the solver?
                    let abs_diff = (self.p.get(i, j) - p_next).abs();
                    max_delta = max_delta.max(abs_diff);

                    // Update pressure at cell (i, j).
                    self.p.set(i, j, p_next);
                }
            }

            if max_delta < epsilon {
                println!("Exiting solver after {} iterations, maximum change is: {}", iter, max_delta);
                break;
            }
        }


        // Update each component of the velocity field based on the
        // pressure gradient. Note that this update should ONLY be
        // applied to components of the velocity that border a grid
        // cell that contains fluid.
        //
        // Some updates may require the pressure of grid cells that
        // lie either outside of the grid or outside of the fluid.
        // We must specify boundary conditions to handle this.
        // a) Dirichlet: at free surface boundaries, the pressure
        //    is zero. We handle this automatically when we initialize
        //    the solver.
        // b) Neumann: at solid walls, we substitute in the solid's
        //    velocity, which is zero for static solids. We handle this
        //    by setting the velocity components along the border
        //    to zero.
        let scale = delta_t  / (self.density * self.grid_cell_size);

        for i in 0..self.w {
            for j in 0..self.h {


                // Set the velocity field along the borders to zero.
                if i == 0 {
                    self.u.set(i, j, 0.0);
                }
                if j == 0 {
                    self.v.set(i, j, 0.0);
                }

                // At the interior cells, apply the pressure gradient.
                if i > 0 && j > 0 {
                    // TODO: check these indices
                    let (grad_p_x, grad_p_y) = self.p.grad(i, j);

                    // Update u-component of velocity
                    let u_next = self.u.get(i, j) - scale * grad_p_x;
                    self.u.set(i, j, u_next);

                    // Update v-component of velocity
                    let v_next = self.v.get(i, j) - scale * grad_p_y;
                    self.v.set(i, j, v_next);
                }
            }
        }



    }




    pub fn update(&mut self) {

        const FRAME_TIME: f64 = 1.0 / 60.0;
        let mut total_t = 0.0;
        let mut total_iters = 0;

        while total_t < FRAME_TIME {
            // 0. Update the hash table of marker cells (i.e. cells that
            // currently contain fluid).
            // TODO

            // 1. Determine a good time step `delta_t`
            let delta_t = self.determine_time_step() * 1000.0; // TODO: this multiplier shouldn't be here

            // 4. Project the velocity field to obey the incompressibility condition:
            //      a. Setup the matrix A of coefficients
            //      b. Setup the vector b, which contains the negative divergence
            //         at each grid cell
            //      c. Use the Gauss-Siedel method to solve for the pressures
            //      d. Apply the pressure (update velocities)
            self.project(delta_t, 1e-5, 6);

            // 2. Update the velocity field (self-advection) and other quantities
            // via backwards particle trace RK2 scheme.
            self.advect(delta_t);

            // 3. Apply body forces (i.e. gravity).
            self.apply_body_forces(delta_t);

            total_t += delta_t;
            total_iters += 1;
        }

        //self.init();
        println!("Update complete - {} solver iterations ran", total_iters);
    }
}
