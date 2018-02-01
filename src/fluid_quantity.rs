use helpers::{Dimension, Vector};

#[derive(Copy, Clone)]
pub enum Staggered {
    /// The quantity is sampled at the center the grid cells
    None,

    /// The quantity is staggered along the x-axis: it will be
    /// sampled along the vertical faces of the grid cells
    OffsetX,

    /// The quantity is staggered along the y-axis: it will be
    /// sampled along the horizontal faces of the grid cells
    OffsetY
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
    /// The width of the grid containing this quantity
    pub w: usize,

    /// The height of the grid containing this quantity
    pub h: usize,

    /// The underlying data store
    data: Vec<f64>,

    /// An enum controlling where the quantity is sampled
    /// inside each grid cell
    pub staggered: Staggered
}

impl FluidQuantity {
    pub fn new(w: usize, h: usize, staggered: Staggered) -> FluidQuantity {
        FluidQuantity {
            w,
            h,
            data: vec![0.0; w * h],
            staggered
        }
    }

    pub fn from_fn<F>(w: usize, h: usize, staggered: Staggered, f: F) -> FluidQuantity
        where F: Fn(usize, usize) -> f64 {

        let mut fluid_quantity = FluidQuantity {
            w,
            h,
            data: vec![0.0; w * h],
            staggered
        };

        for i in 0..fluid_quantity.w {
            for j in 0..fluid_quantity.h {
                fluid_quantity.set(i, j, f(i, j));
            }
        }

        fluid_quantity
    }

    fn check_bounds(&self, i: usize, j: usize) {
        assert!(i >= 0 && j >= 0 && i < self.w && j < self.h);
    }

    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.check_bounds(i, j);
        self.data[i + self.w * j] = v;
    }

    /// Retrieve the quantity at grid cell (i, j): by convention, (0, 0)
    /// is the bottom left corner of the grid
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.check_bounds(i, j);
        self.data[i + self.w * j]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut f64 {
        self.check_bounds(i, j);
        &mut self.data[i + self.w * j]
    }

    /// Returns the gradient of the fluid quantity at grid cell (i, j)
    pub fn grad(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.get(i, j) - self.get(i - 1, j); // Partial w.r.t. x
        let dy = self.get(i, j) - self.get(i, j - 1); // Partial w.r.t. y
        (dx, dy)
    }

    /// Returns the divergence of the fluid quantity at grid cell (i, j)
    pub fn div(&self, i: usize, j: usize) -> f64 {
        self.get(i + 1, j) - self.get(i, j) +
        self.get(i, j + 1) - self.get(i, j)
    }

    /// Returns the laplacian of the fluid quantity at grid cell (i, j)
    pub fn lap(&self, i: usize, j: usize) -> f64 {
        self.get(i + 1, j) + self.get(i - 1, j) +
        self.get(i, j + 1) + self.get(i, j - 1) -
        self.get(i, j) * 4.0
    }
}