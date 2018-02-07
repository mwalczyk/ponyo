extern crate image;

use std::f64;

use helpers::{Vector, lerp};
use fluid_quantity::{Staggered, FluidQuantity};

use image::{GenericImage, ImageBuffer};

pub enum InterpolationScheme {
    /// First-order Runge Kutta interpolation
    RK1,

    /// Second-order Runge Kutta interpolation
    RK2
}

pub struct FluidSolver {
    /// The width of the solver grid
    w: usize,

    /// The height of the solver grid
    h: usize,

    /// The pressure field
    p: FluidQuantity,

    /// The i-component of velocity field
    u: FluidQuantity,

    /// The j-component of velocity field
    v: FluidQuantity,

    /// The dye field, for visualization purposes
    dye: FluidQuantity,

    /// The right-hand side of the pressure solve
    rhs: FluidQuantity,

    /// The density of the fluid
    density: f64,

    /// The duration (in seconds) of each frame of the simulation
    frame_time: f64,

    /// The maximum number of grid cells that a chunk of fluid is
    /// allowed to traverse during a given timestep
    max_grid_cell_traversal: usize,

    /// The size of each grid cell
    grid_cell_size: f64,

    /// A body force that will be applied globally to the entire
    /// fluid
    gravity: (f64, f64),

    /// The interpolation algorithm that will be used during the
    /// backwards particle trace
    interpolation_scheme: InterpolationScheme
}

impl FluidSolver {

    /// Constructs and returns a new solver whose simulation
    /// domain will be `w`x`h` pixels in resolution.
    ///
    /// The u and v-components of the velocity field are one unit
    /// larger along their respective dimensions to form a staggered
    /// MAC (marker-and-cell) grid.
    ///
    /// The MAC grid method discretizes space into a grid of cells.
    /// Each cell has a pressure `p` defined at its center. It also
    /// has a velocity with components `u` and `v`, but the
    /// components are placed at the centers of 2 of the cell faces:
    /// `u` on the x-min face and `v` on the y-min face.
    ///
    /// However, in code, we use integer indices so that:
    ///      p(i, j) = P(i + 0.0, j + 0.0)
    ///      u(i, j) = U(i - 0.5, j + 0.0)
    ///      v(i, j) = V(i + 0.0, j - 0.5)
    pub fn new(w: usize, h: usize) -> FluidSolver {
        println!("Initializing a new FluidSolver...");

        FluidSolver {
            w,
            h,
            p: FluidQuantity::new(w, h, Staggered::None),
            u: FluidQuantity::new(w + 1, h, Staggered::OffsetX),
            v: FluidQuantity::new(w, h + 1, Staggered::OffsetY),
            dye: FluidQuantity::new(w, h, Staggered::None),
            rhs: FluidQuantity::new(w, h, Staggered::None),
            density: 0.1,
            frame_time: 1.0 / 60.0,
            max_grid_cell_traversal: 10,
            grid_cell_size: 1.0 / (w.min(h) as f64),
            gravity: (0.0, 0.0),
            interpolation_scheme: InterpolationScheme::RK2
        }
    }

    /// Sets the density of the fluid.
    pub fn density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// Sets the duration of each frame. For example, to run the simulation at 60
    /// frames per second, set `frame_time` to `1.0 / 60.0`.
    pub fn frame_time(mut self, frame_time: f64) -> Self {
        self.frame_time = frame_time;
        self
    }

    /// Sets the maximum number of grid cells that a chunk of fluid is allowed to
    /// traverse during a given timestep
    pub fn max_grid_cell_traversal(mut self, max_grid_cell_traversal: usize) -> Self {
        self.max_grid_cell_traversal = max_grid_cell_traversal;
        self
    }

    /// Sets the gravity force vector in both the `x` and `y` directions to
    /// `gravity.0` and `gravity.1`, respectively.
    pub fn gravity(mut self, gravity: (f64, f64)) -> Self {
        self.gravity = gravity;
        self
    }

    /// Saves the current state of the solver to disk.
    pub fn to_image(&self, path: &str) {
        let image = ImageBuffer::from_fn(self.w as u32,
                                                                                  self.h as u32,
                                                                                  |x, y| {
            let color = (1.0 - self.dye.get(x as usize, y as usize)) * 255.0;
            image::Luma([color as u8])
        });

        image.save(path).unwrap();
    }

    /// Adds fluid density to some user-specified rectangular region with upper-left
    /// corner (`upper_left_x`, `upper_left_y`) and dimensions `w`x`h`. The fluid
    /// in this region will have horizontal velocity `u` and vertical velocity `v`.
    pub fn add_source(&mut self,
                      mut upper_left_x: usize,
                      mut upper_left_y: usize,
                      mut w: usize,
                      mut h: usize,
                      u: f64,
                      v: f64) {
        // Clamp the coordinates of the upper-left corner, if necessary.
        upper_left_x = upper_left_x.min(self.w).max(0);
        upper_left_y = upper_left_y.min(self.h).max(0);

        for i in upper_left_x..(upper_left_x + w) {
            for j in upper_left_y..(upper_left_y + h) {
                *self.u.get_mut(i, j) = u;
                *self.v.get_mut(i, j) = v;
                *self.dye.get_mut(i, j) = 1.0;
            }
        }
    }

    /// Returns the value of quantity `q` at real-valued coordinates `x` and `y`.
    /// A bilinear interpolation will be performed in order to reconstruct values
    /// that do not lie on exact grid centers.
    ///
    /// Note that `x` and `y` are treated as coordinates on a non-staggered grid.
    /// In other words, `x` = `y` = 0.0 refers to the center of the cell in the
    /// bottom-left corner of the grid.
    fn get_interpolated_quantity(&self, q: &FluidQuantity, mut x: f64, mut y: f64) -> f64 {
        // Clamp the coordinates onto the grid domain.
        x = x.max(0.0).min((self.w - 2) as f64);
        y = y.max(0.0).min((self.h - 2) as f64);

        let mut oi = x - x.floor();
        let mut oj = y - y.floor();

        // Calculate the (integer) indices of the grid cell.
        let mut i = x.floor() as usize;
        let mut j = y.floor() as usize;

        // Fluid quantities that are staggered must be correctly handled here.
        match q.staggered {
            Staggered::OffsetX => oi += 0.5,
            Staggered::OffsetY => oj += 0.5,
            _ => ()
        }

        // Bilinear interpolation of the fluid quantity via the four
        // surrounding grid cells:
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

    /// Update the hash table of marker cells (i.e. cells that actually contain
    /// fluid).
    fn update_marker_cells(&mut self) {
        // TODO
    }

    /// Determines an appropriate timestep, given the current
    /// state of the velocity field.
    fn determine_time_step(&self) -> f64 {
        // Find the length of the largest velocity vector.
        let mut velocity_norm_max = f64::NEG_INFINITY;

        for i in 0..self.w {
            for j in 0..self.h {
                // Get the velocity vector at the center of cell (i, j).
                let velocity_x = self.get_interpolated_quantity(&self.u, i as f64, j as f64);
                let velocity_y = self.get_interpolated_quantity(&self.v, i as f64, j as f64);
                let length = (velocity_x * velocity_x + velocity_y * velocity_y).sqrt();

                if length > velocity_norm_max {
                    velocity_norm_max = length;
                }
            }
        }

        let delta_t = (self.max_grid_cell_traversal as f64 * self.grid_cell_size) / velocity_norm_max;

        delta_t
    }

    /// Advect quantity `q` through the velocity field for
    /// a time interval `delta_t`. This should only be called
    /// with a divergence-free velocity field, i.e. one that
    /// meets the incompressibility constraint.
    fn advect_quantity(&self, delta_t: f64, q: &FluidQuantity) -> FluidQuantity {
        let mut q_next = FluidQuantity::new(q.w, q.h, q.staggered);

        for i in 0..q.w {
            for j in 0..q.h {

                // The starting position of the particle is simply (i, j).
                let mut position = Vector::new(i as f64, j as f64);

                match q.staggered {
                    Staggered::OffsetX => position.x -= 0.5,
                    Staggered::OffsetY => position.y -= 0.5,
                    _ => ()
                }

                // Get the velocity vector at the point where quantity `q` is
                // sampled. For example, the u-component of the velocity field
                // is sampled at the center of the vertical face the grid cell,
                // so we want to trace a particle from this position backwards.
                let velocity = Vector::new(self.get_interpolated_quantity(&self.u, position.x, position.y),
                                                  self.get_interpolated_quantity(&self.v, position.x, position.y));

                // Trace backwards using Runge-Kutta interpolation.
                let position_prev = match self.interpolation_scheme {
                  InterpolationScheme::RK1 => position - velocity * delta_t,
                  InterpolationScheme::RK2 => {
                      let position_mid = position - velocity * delta_t * 0.5;
                      let velocity_mid = Vector::new(self.get_interpolated_quantity(&self.u, position_mid.x, position_mid.y),
                                                            self.get_interpolated_quantity(&self.v, position_mid.x, position_mid.y));
                      position_mid - velocity_mid * delta_t
                  }
                };

                // Set the value of the fluid quantity in the new buffer.
                *q_next.get_mut(i, j) = self.get_interpolated_quantity(q, position_prev.x, position_prev.y);
            }
        }

        q_next
    }

    /// Advect all of the relevant fluid quantities.
    fn advect(&mut self, delta_t: f64) {
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
        //
        // Advection cannot happen in-place: we need to create new
        // buffers and swap them.
        self.u = self.advect_quantity(delta_t, &self.u);
        self.v = self.advect_quantity(delta_t, &self.v);
        self.dye = self.advect_quantity(delta_t, &self.dye);
    }

    /// Apply body forces (such as gravity) to the entire fluid.
    fn apply_body_forces(&mut self) {
        for i in 0..self.w {
            for j in 0..self.h {
                *self.u.get_mut(i, j) += self.gravity.0;
                *self.v.get_mut(i, j) += self.gravity.1;
            }
        }
    }

    /// Build the right-hand side of pressure equation.
    fn build_rhs(&mut self) {
        let scale = 1.0 / self.grid_cell_size;

        // Technically, we only need to do this calculation for cells that
        // contain fluid. Since we aren't currently using markers, we just
        // iterate over the entire field.
        for i in 0..self.w {
            for j in 0..self.h {
                // Calculate the divergence of the velocity field at (i, j).
                let div_u_x = self.u.get(i + 1, j) - self.u.get(i, j);
                let div_v_y = self.v.get(i, j + 1) - self.v.get(i, j);
                let div = div_u_x + div_v_y;

                // The negative divergence becomes the right-hand side of the
                // pressure equation.
                self.rhs.set(i, j, -scale * div);
            }
        }
    }

    /// Project using the Gauss-Seidel method, which is an iterative
    /// approach to solving systems of linear equations.
    fn solve_pressure(&mut self, delta_t: f64, tolerance: f64, max_iterations: usize) {
        let delta_x_squared = self.grid_cell_size * self.grid_cell_size;
        let scale = delta_t / (self.density * delta_x_squared);

        let mut max_delta = 0.0_f64;

        for _ in 0..max_iterations {
            max_delta = 0.0;
            for i in 0..self.w {
                for j in 0..self.h {
                    let mut diagonal = 0.0;
                    let mut off_diagonal = 0.0;

                    // In this function, we are essentially solving an equation
                    // `Ap = b`, where `p` contains the unknown pressure at each
                    // grid cell location.
                    //
                    // First, we assemble the coefficients for the row of `A`
                    // corresponding to cell (i, j). Note that because each cell
                    // has at most 4 neighbors, `A` will be extremely sparse.
                    //
                    // The entry corresponding to p(i, j) will have a coefficient
                    // of the form:
                    //
                    //          [Δt / (ρ * Δx * Δx)] * n
                    //
                    // where `n` is the number of neighboring cells that contain
                    // either fluid or air. This entry would appear along the
                    // diagonal of the matrix `A`, which is why we accumulate
                    // this value in a variable named `diag`.
                    //
                    // Entries corresponding to the cells neighboring cell (i, j)
                    // will have coefficients of the form:
                    //
                    //          -[Δt / (ρ * Δx * Δx)]
                    //
                    // For now, we assume that all solid boundaries are stationary.

                    // Check if this cell is on the left, vertical edge.
                    if i > 0 {
                        diagonal += scale;
                        off_diagonal -= scale * self.p.get(i - 1, j);
                    }
                    // Check if this cell is on the bottom, horizontal edge.
                    if j > 0 {
                        diagonal += scale;
                        off_diagonal -= scale * self.p.get(i, j - 1);
                    }
                    // Check if this cell is on the right, vertical edge.
                    if i < (self.w - 1) {
                        diagonal += scale;
                        off_diagonal -= scale * self.p.get(i + 1, j);
                    }
                    // Check if this cell is on the top, horizontal edge.
                    if j < (self.h - 1) {
                        diagonal += scale;
                        off_diagonal -= scale * self.p.get(i, j + 1);
                    }

                    let p_next = (self.rhs.get(i, j) - off_diagonal) / diagonal;

                    // Can we exit the solver?
                    let abs_diff = (self.p.get(i, j) - p_next).abs();
                    max_delta = max_delta.max(abs_diff);

                    // Update the pressure at cell (i, j).
                    *self.p.get_mut(i, j) = p_next;
                }
            }

            // If none of the values drastically changed, we exit
            // the outer loop.
            if max_delta < tolerance {
                break;
            }
        }
    }

    /// Update each component of the velocity field based on the
    /// pressure gradient. Note that this update should only be
    /// applied to components of the velocity that border a grid
    /// cell that contains fluid.
    ///
    /// Some updates may require the pressure of grid cells that
    /// lie either outside of the grid or outside of the fluid.
    /// We must specify boundary conditions to handle this:
    ///
    /// a) Dirichlet: at free surface boundaries, the pressure
    ///    is zero. We handle this automatically when we initialize
    ///    the solver.
    /// b) Neumann: at solid walls, we substitute in the solid's
    ///    velocity, which is zero for static solids. We handle this
    ///    by setting the velocity components along the border
    ///    to zero.
    fn apply_pressure(&mut self, delta_t: f64) {
        let scale = delta_t  / (self.density * self.grid_cell_size);
        for i in 1..self.w {
            for j in 1..self.h {
                let (grad_p_x, grad_p_y) = self.p.grad(i, j);

                // Update velocity
                *self.u.get_mut(i, j) -= scale * grad_p_x;
                *self.v.get_mut(i, j) -= scale * grad_p_y;
            }
        }

        // Set the velocity field along the borders to zero.
        for j in 0..self.h {
            *self.u.get_mut(0, j) = 0.0;          // Left column
            *self.u.get_mut(self.w - 1, j) = 0.0; // Right column
        }
        for i in 0..self.w {
            *self.v.get_mut(i, 0) = 0.0;          // Bottom row
            *self.v.get_mut(i, self.h - 1) = 0.0; // Top row
        }
    }

    /// Calculate the right amount of pressure to make the velocity
    /// field divergence-free.
    fn project(&mut self, delta_t: f64, tolerance: f64, max_iterations: usize) {
        self.build_rhs();
        self.solve_pressure(delta_t, tolerance, max_iterations);
        self.apply_pressure(delta_t);
    }

    /// Moves the entire fluid simulation forward in time.
    pub fn update(&mut self) {
        let mut total_time = 0.0_f64;

        while total_time < self.frame_time {
            // (0) Update the hash table of marker cells.
            self.update_marker_cells();

            // (1) Determine a good time step, based on the current velocity field.
            let delta_t = self.determine_time_step();

            // (2) Force the fluid to obey the incompressibility condition.
            self.project(delta_t, 1e-5, 600);

            // (3) Update the velocity field (self-advection) and other quantities.
            self.advect(delta_t);

            // (4) Apply body forces.
            self.apply_body_forces();

            total_time += delta_t;
        }

    }
}
