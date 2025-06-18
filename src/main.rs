use std::ops::{Add, Div, Mul};

// x^16 + x^5 + x^3 + x^2 + 1
const IRREDUCIBLE: u32 = 65581;

#[derive(Debug, Clone, Copy)]
struct BinaryElem16 {
    value: u16,
}

fn mod_irreducible(a: u32) -> u16 {
    let field_size: u32 = 16;
    let mut result = a;

    for i in (field_size..32).rev() {
        if result & (1 << i) != 0 {
            result ^= IRREDUCIBLE << (i - field_size);
        }
    }

    result as u16
}

fn poly_multiply_gf2_u16(a: u16, b: u16) -> u32 {
    let mut result: u32 = 0;

    // Karatsuba multiplication algorithm,
    // except that we use XOR instead of addition
    let mut tempa = a as u32;
    let mut tempb = b as u32;
    while tempb != 0 {
        if tempb & 1 != 0 {
            result ^= tempa;
        }
        tempa <<= 1;
        tempb >>= 1;
    }

    result
}

impl BinaryElem16 {
    fn pow(&self, mut exponent: u16) -> Self {
        if exponent == 0 {
            return Self { value: 1 };
        }

        let mut result = Self { value: 1 };
        let mut base = self.clone();

        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exponent >>= 1;
        }

        return result;
    }

    fn inverse(&self) -> Self {
        self.pow(65534)
    }
}

impl Add for BinaryElem16 {
    type Output = BinaryElem16;

    fn add(self, other: BinaryElem16) -> Self::Output {
        BinaryElem16 {
            value: self.value ^ other.value,
        }
    }
}

impl Mul for BinaryElem16 {
    type Output = BinaryElem16;

    fn mul(self, other: BinaryElem16) -> Self::Output {
        let result = poly_multiply_gf2_u16(self.value, other.value);
        let reduced = mod_irreducible(result);

        BinaryElem16 { value: reduced }
    }
}

impl Div for BinaryElem16 {
    type Output = BinaryElem16;

    fn div(self, other: BinaryElem16) -> Self::Output {
        return self * other.inverse();
    }
}

fn main() {
    println!("Hello, world!");
    let a = BinaryElem16 { value: 15000 };
    let b = BinaryElem16 { value: 14198 };
    let c = a + b;
    let d = a * b;
    println!("a: {:?}", a);
    println!("b: {:?}", b);
    println!("addition: {:?}", c);
    println!("mult: {:?}", d);
    println!("inverse: {:?}", a.inverse());
}
