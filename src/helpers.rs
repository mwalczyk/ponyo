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

pub fn hash(i: usize, j: usize, k: usize) -> usize {
    // Hash function presented in [Worley 1996]
    541 * i + 79 * j + 31 * k
}