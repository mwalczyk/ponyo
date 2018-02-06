use std::ops::{Add, Sub, Mul};

#[derive(Copy, Clone)]
pub struct Vector {
    pub x: f64,
    pub y: f64
}

impl Vector {
    pub fn new(x: f64, y: f64) -> Vector {
        Vector{ x, y }
    }

    pub fn ones() -> Vector {
        Vector::new(1.0, 1.0)
    }

    pub fn zeros() -> Vector {
        Vector::new(0.0, 0.0)
    }

    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    pub fn dot(&self, other: &Vector) -> f64 {
        self.x * other.x + self.y * other.y
    }

    pub fn distance(&self, b: Vector) -> f64 {
        (*self - b).length()
    }
}

// Vector-vector addition
impl Add for Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Vector {
        Vector::new(self.x + other.x, self.y + other.y)
    }
}

// Vector-vector subtraction
impl Sub for Vector {
    type Output = Vector;

    fn sub(self, other: Vector) -> Vector {
        Vector::new(self.x - other.x, self.y - other.y)
    }
}

// Vector-scalar multiplication
impl Mul<f64> for Vector {
    type Output = Vector;

    fn mul(self, scalar: f64) -> Vector {
        Vector::new(self.x * scalar, self.y * scalar)
    }
}

/// Returns the linear interpolation of `a` and `b`, given
/// an interpolant `c`.
///
/// In particular, if `c` is 0.0, `a` will be returned. If
/// `c` is 1.0, `b` will be returned.
pub fn lerp(a: f64, b: f64, mut c: f64) -> f64 {
    // Clamp `c`, if necessary.
    c = c.min(1.0).max(0.0);

    (1.0 - c) * a + c * b
}

/// Hash function presented in [Worley 1996].
pub fn hash(i: usize, j: usize, k: usize) -> usize {
    541 * i + 79 * j + 31 * k
}