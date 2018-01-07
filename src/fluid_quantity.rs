use helpers::Dimension;

// Consider making this quantity a generic that
// must implement a trait like "Differentiable"
// ...

pub struct FluidQuantity {
    dims: Dimension,
    data: Vec<f64>
}

impl FluidQuantity {
    pub fn new(dims: Dimension) -> FluidQuantity {
        // Initialize the data buffer
        let data = vec![0.0; dims.nx * dims.ny];

        FluidQuantity { dims, data }
    }

    pub fn at(&self, x: usize, y: usize) -> Option<f64> {
        if x < 0 || y < 0 || x > self.dims.nx || y > self.dims.ny {
            return None
        }

        Some(self.data[x + self.dims.nx * y])
    }
}