# Multiple Sequence Alignment with Speculative Parallelism
This repository contains a C++ implementation of a speculative Berger-Munson MSA algorithm.


## Benchmarking script
To quickly run the code after building, run `./bench.sh` from the `/code` directory. This should take about 5 minute to complete.


Outputs can be seen in `/data/output`.


## Usage
To build the project, run `make` from the `/code` directory.


To run the sequential code,
```
./bm_seq -i input_file -o output_file -r random_mode
```


To run the parallel code,
```
mpirun -np num_procs ./bm_par -i input_file -o output_file -r random_mode
```

The `-r` flag allows for pseudorandomness, to facilitate benchmarking and debugging. To use pseudorandomness, pass `-r P`. To use "true" randomness, pass `-r R` (or omit the flag).

# Documentation
Run `make docs`. Results are in `/code/docs/`.
