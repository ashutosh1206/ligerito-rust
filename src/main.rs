use cryptoutils::BinaryElem16;

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
