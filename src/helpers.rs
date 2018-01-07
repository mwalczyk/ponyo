#[derive(Copy, Clone)]
pub struct Dimension {
    pub nx: usize,
    pub ny: usize
}

impl Dimension {
    pub fn new(nx: usize, ny: usize) -> Dimension {
        Dimension{ nx, ny }
    }
}

pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    0.0
}