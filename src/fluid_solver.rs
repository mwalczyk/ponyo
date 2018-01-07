use helpers::Dimension;
use fluid_quantity::FluidQuantity;

pub struct FluidSolver {
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
            p: FluidQuantity::new(dims),
            u: FluidQuantity::new(Dimension::new(dims.nx + 1, dims.ny)),
            v: FluidQuantity::new(Dimension::new(dims.nx, dims.ny + 1))
        }
    }

    fn init() {
        // Create an initial divergence-free velocity field `u`
    }

    fn project(&mut self, delta_t: f64, u: &mut FluidQuantity) {
        // Calculate and apply just the right amount of pressure to
        // make `u` divergence-free
    }

    fn advect(&mut self, u: &FluidQuantity, delta_t: f64, q: &mut FluidQuantity) {
        // Advect quantity `q` through the velocity field `u` for
        // a time interval `delta_t`. This should ONLY be called
        // with a divergence-free velocity field `u`, i.e. one that
        // meets the incompressibility constraint.
        //
        // Here, we take a semi-Lagrangian approach.
    }

    pub fn update(&mut self) {
        // 1. Determine a good time step `delta_t`
        let delta_t = 0.0;


        // 2. Update the velocity field (self-advection)


        // 3. Add body forces (i.e. gravity)


        // 4. Project the velocity field to obey the incompressibility condition
    }
}
