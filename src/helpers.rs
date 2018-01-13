use std::ops::{Add, Sub, Mul};

#[derive(Copy, Clone)]
pub struct Dimension {
    pub nx: usize,
    pub ny: usize
}

impl Dimension {
    pub fn new(nx: usize, ny: usize) -> Dimension {
        Dimension{ nx, ny }
    }

    pub fn expand(&self, ix: usize, iy: usize) -> Dimension {
        Dimension::new(self.nx + ix, self.ny + iy)
    }
}

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

pub fn lerp(a: f64, b: f64, mut c: f64) -> f64 {
    // Returns the linear interpolation of `a` and `b`, given
    // an interpolant `c`. In particular, if `c` is 0.0, `a`
    // will be returned. If `c` is 1.0, `b` will be returned.
    c = c.min(1.0).max(0.0);

    (1.0 - c) * a + c * b
}

pub fn hash(i: usize, j: usize, k: usize) -> usize {
    // Hash function presented in [Worley 1996]
    541 * i + 79 * j + 31 * k
}

pub fn gauss_seidel(tolerance: f64, max_iterations: usize) -> f64 {
    // Source: https://www3.nd.edu/~zxu2/acms40390F12/Lec-7.3.pdf

    // 1. Set all elements of the vector x to 0
    let mut x_i = 0.0;
    let mut x_k = x_i;

    // 2. Outer loop
    for i in 0..max_iterations {
        // 3. Update x_i using x_k
        // TODO

        // 4. Calculate the length of the vector ||x_i - x_k||
        // TODO
        let length = 0.0;
        if length < tolerance {
            return x_i;
        }

        // 5. Set x_k to x_i
        x_k = x_i;
    }

    x_i
}