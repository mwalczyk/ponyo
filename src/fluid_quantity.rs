use helpers::{Dimension, Vector};

#[derive(Copy, Clone)]
pub enum Staggered {
    None,
    OffsetX,
    OffsetY
}

impl Staggered {
    pub fn as_offset(&self) -> Vector {
        match *self {
            Staggered::None => Vector::new(0.0, 0.0),
            Staggered::OffsetX => Vector::new(-0.5, 0.0),
            Staggered::OffsetY => Vector::new(0.0, -0.5)
        }
    }
}

#[derive(Copy, Clone)]
pub enum Cell {
    Fluid{ q: f64 },     // A grid cell containing a fluid quantity `q`
    Solid,               // A grid cell at a solid boundary (i.e. wall)
    Empty{ p: f64 }      // A grid cell containing air, with pressure `p`
}

// Consider making this quantity a generic that
// must implement a trait like "Differentiable"
// ...
#[derive(Clone)]
pub struct FluidQuantity {
    // The width and height of the grid
    pub dims: Dimension,

    // The underlying data store
    data: Vec<f64>,

    // An enum controlling where the quantity is sampled
    // inside each grid cell
    pub staggered: Staggered
}

impl FluidQuantity {
    pub fn new(dims: Dimension, staggered: Staggered) -> FluidQuantity {
        FluidQuantity {
            dims,
            data: vec![0.0; dims.nx * dims.ny],
            staggered
        }
    }

    pub fn from_fn<F>(dims: Dimension, staggered: Staggered, f: F) -> FluidQuantity
        where F: Fn(usize, usize) -> f64 {

        let mut fluid_quantity = FluidQuantity {
            dims,
            data: vec![0.0; dims.nx * dims.ny],
            staggered
        };

        for i in 0..fluid_quantity.dims.nx {
            for j in 0..fluid_quantity.dims.ny {
                fluid_quantity.set(i, j, f(i, j));
            }
        }

        fluid_quantity
    }

    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        assert!(i > 0 && j > 0 && i < self.dims.nx && j < self.dims.ny);
        self.data[i + self.dims.nx * j] = v;
    }

    // Retrieve the quantity at grid cell (i, j): by convention, (0, 0)
    // is the bottom left corner of the grid
    pub fn at(&self, i: usize, j: usize) -> f64 {
        assert!(i > 0 && j > 0 && i < self.dims.nx && j < self.dims.ny);
        self.data[i + self.dims.nx * j]
    }

    // Returns the gradient of the fluid quantity at grid cell (i, j)
    pub fn grad(&self, i: usize, j: usize) -> Vector {
        Vector::new(self.at(i, j) - self.at(i - 1, j), // Partial w.r.t. x
                    self.at(i, j) - self.at(i, j - 1)) // Partial w.r.t. y
    }

    // Returns the divergence of the fluid quantity at grid cell (i, j)
    pub fn div(&self, i: usize, j: usize) -> f64 {
        self.at(i + 1, j) - self.at(i, j) +
        self.at(i, j + 1) - self.at(i, j)
    }

    // Returns the laplacian of the fluid quantity at grid cell (i, j)
    pub fn lap(&self, i: usize, j: usize) -> f64 {
        self.at(i + 1, j) + self.at(i - 1, j) +
        self.at(i, j + 1) + self.at(i, j - 1) -
        self.at(i, j) * 4.0
    }
}