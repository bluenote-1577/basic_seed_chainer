# basic_seed_chainer

Experiments for TODO by implementing a basic seed-chain-extend algorithm with sketching in rust. The program generates a random string,
a mutated version of the random string with error rate theta, and aligns them via seed-chain-extend. 

**Seeding**: using open syncmer or minimizer seeds. NOTE: We use fxhash for k-mer hashing instead of a rolling hash so it's not as fast as it could be. 

**Chaining**: using a linear gap cost with a MRQ data structure as described in the minigraph paper.

**Extension**: using rust-bio's simple dynamic programming extension algorithm or wavefront alignment (WFA). 

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
### Install
```
git clone https://github.com/bluenote-1577/basic_seed_chainer
cd basic_seed_chainer
cargo build --release
./target/release/basic_seed_chainer
```

### Parameters

TODO

### Results

TODO
