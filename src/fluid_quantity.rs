use helpers::{Dimension, Vector};

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

    pub fn at(&self, i: usize, j: usize) -> f64 {
        assert!(i > 0 && j > 0 && i < self.dims.nx && j < self.dims.ny);
        // The MAC grid method discretizes space into a grid of cells.
        // Each cell has a pressure `p` defined at its center. It also
        // has a velocity `u` (with components `u` and `v`), but the
        // components are placed at the centers of 2 of the cell faces:
        // `u` on the x-min face and `v` on the y-min face.
        //
        // Syntactically, this means:
        // p(i, j) = P(i + 0.0, j + 0.0)
        // u(i, j) = U(i - 0.5, j + 0.0)
        // v(i, j) = V(i + 0.0, j - 0.5)
        //
        // For a quantity that is staggered along the x-axis, we can
        // reconstruct the value at the grid center as follows:
        // u(i, j) = (u(i - 0.5, j) + u(i + 0.5, j)) / 2
        //
        // A similar process applies to a quantity that is staggered
        // along the y-axis.
        let offset = self.staggered.as_offset();

        self.data[i + self.dims.nx * j]
    }

    pub fn vector_at(&self, i: usize, j: usize)

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