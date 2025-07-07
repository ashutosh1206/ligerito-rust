use crate::binary_field::BinaryField;
use std::ops::{Add, Mul};

pub struct QuadraticEvals<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    e0: F,
    e1: F,
    e2: F,
}

pub struct QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    a: F,
    b: F,
    c: F,
}

impl<F> QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    // TODO: do we really need this function if we are not performing
    // any operation on a,b,c and simply initializing the struct?
    pub fn new(a: F, b: F, c: F) -> Self {
        Self { a, b, c }
    }

    pub fn eval_quadratic(&self, r: F) -> F {
        self.a * r * r + self.b * r + self.c
    }
}

pub fn quadratic_from_evals<F>(at0: F, at1: F, atx: F, x: Option<F>) -> QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    // TODO: Can we unwrap x in a simpler way from Option<F>?
    let x = x.unwrap_or_else(|| F::new(F::ValueType::try_from(3).unwrap()));
    let numerator = atx + at0 + x * (at1 + at0);
    let denominator = x * x + x;
    // TODO: is there a better way to do this than .unwrap()?
    let a = numerator + denominator.inverse().unwrap();
    let b = at1 + at0 + a;
    return QuadraticPoly { a, b, c: at0 };
}

pub fn fold_quadratic<F>(p1: QuadraticPoly<F>, p2: QuadraticPoly<F>, alpha: F) -> QuadraticPoly<F>
where
    F: Copy + BinaryField + Add<Output = F> + Mul<Output = F>,
    F::ValueType: TryFrom<usize>,
    <F::ValueType as TryFrom<usize>>::Error: std::fmt::Debug,
{
    QuadraticPoly::new(
        p1.a + alpha * p2.a,
        p1.b + alpha * p2.b,
        p1.c + alpha * p2.c,
    )
}
