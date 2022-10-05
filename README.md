# basic_seed_chainer

Experiments for TODO by implementing a basic seed-chain-extend algorithm with sketching in rust. The program generates a random string,
a mutated version of the random string with error rate theta, and aligns them via seed-chain-extend. Running times and recoverability (as defined in our paper) are output after the experiments are run. 

The value of the k-mers is increasing as `k = C log n` where `n` is the sequence length, and `C` is defined to be `2/(1 - 2 * alpha)` where `alpha = -log (1 - theta)` with log base 4. We simulate alignments for `k = 10, 11, 12, 13, ...` up to a user specified value. Other simulation parameters can be specified and are outlined below. For a quick overview of the algorithm:

**Seeding**: using open syncmer or minimizer seeds. NOTE: We don't use any sort of bitwise algorithms for representing k-mers. We have not optimized for k-mer matching, seeding, etc so it will be slow.  

**Chaining**: using a linear gap cost with a MRQ data structure as described in the minigraph paper.

**Extension**: using rust-bio's simple dynamic programming extension algorithm or wavefront alignment (WFA). 

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
### Install
```
git clone https://github.com/bluenote-1577/basic_seed_chainer
cd basic_seed_chainer
cargo build --release
./target/release/basic_seed_chainer -h
#10 iterations at theta 0.05 for k = 11,...,19
./taget/release/basic_seed_chainer 10 0.05 
```

### Parameters

One can specify a few parameters for the simulation described by the `basic_seed_chainer -h` menu. Here we'll list a few of importance.

1. `-s`: performs sketching with open syncmers with `c`, the reciprocal of fraction of selected k-mers, equal to `c = k - 6`. 
2. `--wfa`: use wavefront alignment (WFA) extension instead of slow generic DP extension. 
3. `-k INT`: maximum k-mer size to iterate up to. 

### Results
Outputs are of the following form:
``[src/main.rs:400] extend_cumulative = [
    9.04046e-5,
    0.0001350113,
    0.00023640381,
    0.0003151813,
    0.00067358697,
    0.0014464213,
    0.003459784,
    0.008586039,
    0.018265765,
    0.041963242,
]``

where the mean runtimes over all iterations of `k` are output. Cumulative chaining times and also recoverability are output. 
