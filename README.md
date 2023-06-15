# haspac
Search for **Ha**damard **S**ubmatrices in **Pa**ley **C**ores

**Summary**
HasPac is a simple set of routines to search for Hadamard matrices inside of Paley cores.
For now it is only able to do a randomized search.

HasPac makes use of the library Cliquer, created by Niskanen and Östergård: https://users.aalto.fi/~pat/cliquer.html

To compile run
make haspac CL_DIR=/path/to/cliquer

The program haspac has the following usage:
haspac <hadamard_order> <n_tries> <input_file>

The program then will read a list of prime numbers from <input_file>,
and check if there is an Hadamard matrix of order <hadamard_order> in the p-th paley matrix `paley(p)`, for each prime p in <input_file>.
If `r` is an array of length n, and `s` is a subset of {1,2,...,n}, denote by `r[s]` the subarray of s consisting of those indices taken from `s` in order. Haspac carries a randomized search, consisting of three steps:

1. Randomly select a subset of columns, say `cols` , of size <hadamard_order> from `paley(p)`
2. For each row `r_i` of `paley(p)`, create the orthogonality graph `g` of the restricted rows `r_i[cols]`
3. Call cliquer to find a clique of size <hadamard_order> in `g`

