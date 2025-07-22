use std::ops::{Add, Mul, Sub, Div, Neg};
use nalgebra::{Vector2,};

use crate::consts::*;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DV2{
    pub x : f64,
    pub y : f64,
}

impl DV2{
    pub fn new(x : f64, y : f64) -> Self{
        DV2 { x , y }
    }
    pub fn from_car(xy_c : Vector2<f64>) -> Self{
        let x = xy_c.dot(&Vector2::new(0.0, -2.0/SQRT_3));
        let y = xy_c.dot(&Vector2::new(1.0,  1.0/SQRT_3));

        DV2 { x , y }
    }
    pub fn to_car(self) -> Vector2<f64>{
        self.x * Vector2::new(0.5,-SQRT_3/2.) + self.y * Vector2::new(1.0,0.0)
    }
    pub fn from_car2(xy_c : Vector2<f64>) -> Self{
        let x = xy_c.dot(&Vector2::new(1./SQRT_3, -1.));
        let y = xy_c.dot(&Vector2::new(1./SQRT_3,  1.));

        DV2 { x , y }
    }
    pub fn to_car2(self) -> Vector2<f64>{
        self.x * Vector2::new(SQRT_3/2., -0.5) + self.y * Vector2::new(SQRT_3/2., 0.5)
    }
}

impl Add for DV2 {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Mul<f64> for DV2 {
    type Output = Self; 

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Sub for DV2 {
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Div<f64> for DV2 {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Mul<DV2> for f64 {
    type Output = DV2;
    fn mul(self, rhs: DV2) -> Self::Output {
        rhs * self
    }
}

impl Neg for DV2 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}