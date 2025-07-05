use std::ops::{Add, Div, Mul};

pub trait BinaryField: Sized {
    type ValueType;
    fn zero() -> Self;
    fn one() -> Self;
    fn new(value: Self::ValueType) -> Self;
    fn inverse(&self) -> Result<Self, &'static str>;
}

macro_rules! define_galois_field {
    ($struct_name: ident, $value_type: ty, $irreducible:expr, $max_exp:expr) => {
        #[derive(Debug, Clone, Copy)]
        pub struct $struct_name {
            pub value: $value_type,
        }

        impl $struct_name {
            pub fn pow(&self, mut exponent: $value_type) -> Self {
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

            pub fn inverse(&self) -> Result<Self, &'static str> {
                if self.value == 0 {
                    Err("Cannot compute inverse of zero")
                } else {
                    // By Fermat's little theorem we have: self^{p-1} = 1 (mod p)
                    // Which gives us: self^{-1} = self^{2^16-2} for BinaryElem16, for example
                    // TODO: see if there are faster alternatives to computing inverses
                    // max_exp is 65534 for BinaryElem16
                    Ok(self.pow($max_exp))
                }
            }

            fn mod_irreducible(a_hi: $value_type, a_lo: $value_type) -> $value_type {
                let field_size_bits = std::mem::size_of::<$value_type>() * 8;

                let mut result_hi = a_hi;
                let mut result_lo = a_lo;
                let irreducible = $irreducible as $value_type;

                for i in (field_size_bits..(2 * field_size_bits)).rev() {
                    let bit_set = if i >= field_size_bits {
                        result_hi & (1 << (i - field_size_bits)) != 0
                    } else {
                        result_lo & (1 << i) != 0
                    };

                    if bit_set {
                        let shift_amount = i - field_size_bits;
                        if shift_amount == 0 {
                            result_lo ^= irreducible;
                        } else if shift_amount < field_size_bits {
                            let shifted_hi = irreducible >> (field_size_bits - shift_amount);
                            let shifted_lo = irreducible << shift_amount;

                            result_hi ^= shifted_hi;
                            result_lo ^= shifted_lo;
                        }
                    }
                }

                result_lo
            }

            fn poly_multiply_gf2(a: $value_type, b: $value_type) -> ($value_type, $value_type) {
                let mut result_hi: $value_type = 0;
                let mut result_lo: $value_type = 0;

                let mut tempa_hi: $value_type = 0;
                let mut tempa_lo: $value_type = a;
                let mut tempb = b;
                let bits = std::mem::size_of::<$value_type>() * 8;

                while tempb != 0 {
                    if tempb & 1 != 0 {
                        result_lo ^= tempa_lo;
                        result_hi ^= tempa_hi;
                    }

                    tempa_hi = (tempa_hi << 1) | (tempa_lo >> (bits - 1));
                    tempa_lo = tempa_lo << 1;

                    tempb >>= 1;
                }

                (result_hi, result_lo)
            }
        }

        impl Add for $struct_name {
            type Output = Self;

            fn add(self, other: Self) -> Self::Output {
                Self {
                    value: self.value ^ other.value,
                }
            }
        }

        impl PartialEq for $struct_name {
            fn eq(&self, other: &Self) -> bool {
                self.value == other.value
            }
        }

        impl Div for $struct_name {
            type Output = Self;

            fn div(self, other: Self) -> Self::Output {
                return self * other.inverse().expect("Cannot divide by zero");
            }
        }

        impl Mul for $struct_name {
            type Output = Self;

            fn mul(self, other: Self) -> Self::Output {
                let (result_hi, result_lo) = Self::poly_multiply_gf2(self.value, other.value);
                let reduced = Self::mod_irreducible(result_hi, result_lo);
                Self { value: reduced }
            }
        }

        impl BinaryField for $struct_name {
            type ValueType = $value_type;

            fn zero() -> Self {
                Self { value: 0 }
            }

            fn one() -> Self {
                Self { value: 1 }
            }

            fn new(value: Self::ValueType) -> Self {
                Self { value }
            }

            fn inverse(&self) -> Result<Self, &'static str> {
                self.inverse()
            }
        }
    };
}

// irreducible polynomial: x^16 + x^5 + x^3 + x^2 + 1 (stored without x^16 term)
// max_exp: 2^16-2
define_galois_field!(BinaryElem16, u16, 0x2D, 65534);
// irreducible polynomial: x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1 (stored without x^32 term)
// max_exp: 2^32-2
define_galois_field!(BinaryElem32, u32, 0x8299, 4294967294);
// irreducible polynomial: x^128 + x^7 + x^2 + x + 1 (stored without x^128 term)
// max_exp: 2^128-2
define_galois_field!(
    BinaryElem128,
    u128,
    0x87,
    340282366920938463463374607431768211454
);

#[cfg(test)]
mod tests {
    use super::*;
    /// Macro to generate comprehensive tests for all field sizes
    macro_rules! test_galois_field {
        ($field_type:ty, $value_type:ty, $test_suffix:ident) => {
            paste::paste! {
                #[test]
                fn [<test_add_ $test_suffix>]() {
                    const ONE: $field_type = $field_type { value: 1 };
                    const ZERO: $field_type = $field_type { value: 0 };
                    const MAX_VAL: $field_type = $field_type { value: $value_type::MAX };

                    let random1 = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    let random2 = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    let random3 = $field_type {
                        value: rand::random::<$value_type>(),
                    };

                    assert_eq!(ONE + ONE, ZERO);
                    assert_eq!(random1 + ZERO, random1);
                    // Ensure no overflows when adding in GF(2^16)
                    assert_eq!(
                        MAX_VAL + ONE,
                        $field_type {
                            value: $value_type::MAX - 1
                        }
                    );
                    // Commutative property
                    assert_eq!(random1 + random2, random2 + random1);
                    // Associative property
                    assert_eq!(random1 + (random2 + random3), (random1 + random2) + random3);
                }

                #[test]
                fn [<test_mul_ $test_suffix>]() {
                    const ONE: $field_type = $field_type { value: 1 };
                    const ZERO: $field_type = $field_type { value: 0 };

                    let random1 = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    let random2 = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    let random3 = $field_type {
                        value: rand::random::<$value_type>(),
                    };

                    assert_eq!(ZERO * ZERO, ZERO);
                    assert_eq!(random1 * ZERO, ZERO);
                    assert_eq!(ONE * ONE, ONE);
                    assert_eq!(random1 * ONE, random1);

                    // Commutative property
                    assert_eq!(random1 * random2, random2 * random1);
                    // Associative property
                    assert_eq!(random1 * (random2 * random3), (random1 * random2) * random3);
                }

                #[test]
                fn [<test_inverse_ $test_suffix>]() {
                    const ONE: $field_type = $field_type { value: 1 };
                    const ZERO: $field_type = $field_type { value: 0 };
                    assert_eq!(ONE.inverse().unwrap(), ONE);

                    let result = ZERO.inverse();
                    assert!(result.is_err());
                    assert_eq!(
                        result.unwrap_err(),
                        "Cannot compute inverse of zero"
                    );

                    let random = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    assert_eq!(random * random.inverse().unwrap(), ONE);
                }

                #[test]
                fn [<test_div_ $test_suffix>]() {
                    const ONE: $field_type = $field_type { value: 1 };
                    const ZERO: $field_type = $field_type { value: 0 };
                    let random1 = $field_type {
                        value: rand::random::<$value_type>(),
                    };

                    assert_eq!(random1 / ONE, random1);

                    let panic_result = std::panic::catch_unwind(|| random1 / ZERO);
                    assert!(panic_result.is_err(), "Division by 0 should have panicked");
                }

                #[test]
                fn [<test_ops_combination_ $test_suffix>]() {
                    let random1 = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    let random2 = $field_type {
                        value: rand::random::<$value_type>(),
                    };
                    let random3 = $field_type {
                        value: rand::random::<$value_type>(),
                    };

                    // Distributive property
                    assert_eq!(
                        random1 * (random2 + random3),
                        (random1 * random2) + (random1 * random3)
                    );
                }
            }
        };
    }
    test_galois_field!(BinaryElem16, u16, binary_elem_16);
    test_galois_field!(BinaryElem32, u32, binary_elem_32);
    test_galois_field!(BinaryElem128, u128, binary_elem_128);
}
