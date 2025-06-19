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

    fn inverse(&self) -> Result<Self, &'static str> {
        if self.value == 0 {
            Err("Cannot compute inverse of zero in GF(2^16)")
        } else {
            // By Fermat's little theorem we have: self^{p-1} = 1 (mod p)
            // Hence we have: self^{-1} = self^{2^16-2}
            // TODO: see if there are faster alternatives to computing inverses
            Ok(self.pow(65534))
        }
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
        return self * other.inverse().expect("Cannot divide by zero");
    }
}

impl PartialEq for BinaryElem16 {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

#[cfg(test)]
mod tests {
    use std::u16;

    use crate::BinaryElem16;
    use rand;

    const ONE: BinaryElem16 = BinaryElem16 { value: 1 };
    const ZERO: BinaryElem16 = BinaryElem16 { value: 0 };
    const MAX_VAL: BinaryElem16 = BinaryElem16 { value: u16::MAX };

    #[test]
    fn test_add() {
        let random1 = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        let random2 = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        let random3 = BinaryElem16 {
            value: rand::random::<u16>(),
        };

        assert_eq!(ONE + ONE, ZERO);
        assert_eq!(random1 + ZERO, random1);
        // Ensure no overflows when adding in GF(2^16)
        assert_eq!(
            MAX_VAL + ONE,
            BinaryElem16 {
                value: u16::MAX - 1
            }
        );
        // Commutative property
        assert_eq!(random1 + random2, random2 + random1);
        // Associative property
        assert_eq!(random1 + (random2 + random3), (random1 + random2) + random3);
    }

    #[test]
    fn test_mul() {
        let random1 = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        let random2 = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        let random3 = BinaryElem16 {
            value: rand::random::<u16>(),
        };

        assert_eq!(ZERO * ZERO, ZERO);
        assert_eq!(random1 * ZERO, ZERO);
        assert_eq!(ONE * ONE, ONE);
        assert_eq!(random1 * ONE, random1);
        assert_eq!(MAX_VAL * MAX_VAL, BinaryElem16 { value: 0x5419 });

        // Commutative property
        assert_eq!(random1 * random2, random2 * random1);
        // Associative property
        assert_eq!(random1 * (random2 * random3), (random1 * random2) * random3);
    }

    #[test]
    fn test_inverse() {
        assert_eq!(ONE.inverse().unwrap(), ONE);

        let result = ZERO.inverse();
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Cannot compute inverse of zero in GF(2^16)"
        );

        assert_eq!(
            BinaryElem16 { value: 15000 }.inverse().unwrap(),
            BinaryElem16 { value: 0xc3a0 }
        );

        let random = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        assert_eq!(random * random.inverse().unwrap(), ONE);
    }

    #[test]
    fn test_div() {
        let random1 = BinaryElem16 {
            value: rand::random::<u16>(),
        };

        assert_eq!(random1 / ONE, random1);

        let panic_result = std::panic::catch_unwind(|| random1 / ZERO);
        assert!(panic_result.is_err(), "Division by 0 should have panicked");
    }

    #[test]
    fn test_ops_combination() {
        let random1 = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        let random2 = BinaryElem16 {
            value: rand::random::<u16>(),
        };
        let random3 = BinaryElem16 {
            value: rand::random::<u16>(),
        };

        // Distributive property
        assert_eq!(
            random1 * (random2 + random3),
            (random1 * random2) + (random1 * random3)
        );
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
