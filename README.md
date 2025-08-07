# Deterministic Hashing – Bachelor Thesis Code

This repository contains the source code developed for my bachelor thesis named "Deterministic Hashing".

## Source Files (`src/`)

The repository contains four C++ implementations of static dictionary construction algorithms:

- `dd_rand.cpp` – Double displacement procedure (randomized version).  
  Hashes keys without collisions.  
  Expected construction time: **O(n)**.

- `dd_derand.cpp` – Double displacement procedure (derandomized version).  
  Hashes keys without collisions.  
  Construction time: **O(n log n)**.

- `ruz_1.cpp` – Ružić's slower construction.  
  Allows up to *n* collisions if |U| ≤ 2<sup>⌈2log N - 2loglogN⌉</sup>.  
  Construction time: **O(n log² n)**.

- `ruz_2.cpp` – Ružić's faster construction.  
  Allows up to *Φ<sup>2</sup>n* collisions under the same universe constraint.  
  Construction time: **O(n log log n)**.

---

## Test Datasets (`tests/`)

The `tests` directory includes three datasets:

- `correctness/`  
  Test cases with large keys, used to verify `dd_rand` and `dd_derand`.

- `correctness_small/`  
  Keys are small enough that universe reduction is not required.  
  Used for testing `ruz_1` and `ruz_2`.

- `performance/`  
  Large-scale benchmarks, possibly with large key values.  
  Used to assess scalability and performance.

> Test case filenames use a numeric prefix:  
> A file named `10_s0.in` contains *2<sup>10</sup>* keys, `11_s0.in` contains *2<sup>11</sup>* keys, and so on.

---

##  Running the Code

To run any algorithm on a dataset, use the following command format:

```bash
make check-<dataset>-<algorithm_name>
```
The following command runs the selected algorithm on all test cases from the specified dataset and displays the number of collisions this hashing yielded.
