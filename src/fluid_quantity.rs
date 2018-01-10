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
    staggered: Staggered
}

impl FluidQuantity {
    pub fn new(dims: Dimension, staggered: Staggered) -> FluidQuantity {
        FluidQuantity {
            dims,
            data: vec![0.0; dims.nx * dims.ny],
            staggered
        }
    }

    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        assert!(i > 0 && j > 0 && i < self.dims.nx && j < self.dims.ny);
        self.data[i + self.dims.nx * j] = v;
    }

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