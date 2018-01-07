use helpers::Dimension;
use fluid_quantity::FluidQuantity;

pub struct FluidSolver {
    dims: Dimension,
    p: FluidQuantity, // Pressure field
    u: FluidQuantity, // i-component of velocity field
    v: FluidQuantity, // j-component of velocity field
}

impl FluidSolver {
    pub fn new(dims: Dimension) -> FluidSolver {
        println!("Initializing a new FluidSolver...");

        // The u and v-components of the velocity field are one unit
        // larger along their respective dimensions to form a staggered
        // MAC (marker-and-cell) grid.
        //
        // Syntactically, this means:
        // p(i, j) = P(i, j)
        // u(i, j) = U(i - 1/2, j)
        // v(i, j) = V(i, j - 1/2)
        FluidSolver {
            dims,
            p: FluidQuantity::new(dims),
            u: FluidQuantity::new(Dimension::new(dims.nx + 1, dims.ny)),
            v: FluidQuantity::new(Dimension::new(dims.nx, dims.ny + 1))
        }
    }

    fn init() {
        // Create an initial divergence-free velocity field `u`.
        // TODO
    }

    fn velocity_at(&self, i: usize, j: usize) {
        // Returns the velocity at grid cell (i, j). Since we are
        // using a staggered grid, the velocity vector must be
        // reconstructed by performing linear interpolations of
        // the surrounding grid cells.
        // TODO
    }

    fn project(&mut self, delta_t: f64, u: &mut FluidQuantity) {
        // Calculate and apply just the right amount of pressure to
        // make `u` divergence-free.
        // TODO
    }


    fn advect(&mut self, u: &FluidQuantity, delta_t: f64, q: &mut FluidQuantity) {
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
        // field `u`.
        // TODO
        let u_max = 0.0;

        // Treat each pixel as a grid cell: here `delta_x` represents
        // the width of an individual cell.
        let delta_x = 1.0 / (self.dims.nx as f64);

        // The fluid should not move more than 5 grid cells per
        // iteration.
        const MAX_GRID_CELL_TRAVERSAL: usize = 5;
        let delta_t = (MAX_GRID_CELL_TRAVERSAL as f64 * delta_x) / u_max;

        delta_t
    }

    pub fn update(&mut self) {
        // 1. Determine a good time step `delta_t`
        let delta_t = self.determine_time_step();


        // 2. Update the velocity field (self-advection)
        // TODO

        // 3. Add body forces (i.e. gravity)
        // TODO

        // 4. Project the velocity field to obey the incompressibility condition
        // TODO
    }
}
