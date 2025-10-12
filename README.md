# Ligerito Rust implementation

Pure Rust implementation of Ligerito: https://angeris.github.io/papers/ligerito.pdf

Currently, the project mirrors Julia reference: https://github.com/bcc-research/Ligerito.jl and contains other zk utilities from: https://github.com/bcc-research/ligerito-impl

## Installation

- Build the code: `cargo build --release`
- Run the code for polynomials of coefficient size $2^{20}$, with release builds enabled: `cargo run --release`

> [!WARNING]
> This is an academic prototype, and is NOT ready for production use.