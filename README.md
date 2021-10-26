Simple R package that exports a single C++ function called `wcCpp()`. It calculates overall and per-locus F-statistics (F<sub>IT</sub>, F<sub>ST</sub>, and F<sub>IS</sub>) as well as pairwise F<sub>ST</sub>. It uses the formulas given in Weir and Cockerham (1984).

I wrote this because existing functions in R (i.e. `Fst()` from [`pegas`](https://cran.r-project.org/web/packages/pegas/index.html); `wc()` and `pairwise.WCfst()` from 
[`hierfstat`](https://cran.r-project.org/web/packages/hierfstat/index.html)) are not fast - this function runs considerably faster.

Below are the results of benchmarking the functions using a dataset of 1937 individuals with 20 loci, 289 alleles, and 5 populations.

```
  expression                      median `itr/sec` mem_alloc n_itr total_time
1 wcCpp(alls, pp)                35.94ms    27.7     465.8KB    14   504.83ms
2 pegas::Fst(loc)                  2.36s     0.424    64.7MB     1      2.36s
3 hierfstat::wc(hf)                 2.1s     0.477    86.9MB     1       2.1s
4 hierfstat::pairwise.WCfst(hf)    9.88s     0.101   442.8MB     1      9.88s
```
Weir, B. S. and Cockerham, C. C. (1984) Estimating F-statistics for the analysis of population structure. Evolution, 38, 1358–1370.
