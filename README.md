# MonetDB-Bindex

Bindex is a two-layered index structure based on binned bitmaps that can be used to significantly accelerate the scan operations for in-memory column stores.[1]

We intergrate Bindex into MonetDB-11.33.11 (Apr2019-SP1). For more information on MonetDB, we refer to the original [README](README.monetdb.rst) and [website](https://www.monetdb.org).

Currently we support:

- single table selection
- integer value selection
- sequential pipeline optimizer[2]
- MAL[3] interface

## Build & Install

Please refer to [HowToStart.rst](HowToStart.rst) for build prerequisites and other details.

```bash
git clone https://github.com/zingdle/SIGMOD2020-MonetDB-Bindex.git

cd MonetDB-Bindex

./bootstrap

mkdir build && cd build

../configure --prefix=<your-install-path> --disable-strict --disable-assert --disable-debug --enable-optimize

make -j

make install
```

## Tests

We add serveral Bindex tests in `monetdb5/modules/kernel/Tests`:

```
abbdxselect.malc                      # anti between
btbdxselect.malc                      # between
gebdxselect.malc                      # >=
gtbdxselect.malc                      # >
lebdxselect.malc                      # <=
ltbdxselect.malc                      # <
sk*bdxselect.malc                     # skew data
```

Each test is written in `MAL`:

```
b := bat.new(:int);                   # Create a BAT[4]
bat.append(b, 27645);                 # Append data
...
bat.bindex(b);                        # Build Bindex
s := algebra.bdxselect(b, args);      # Selection
res := algebra.bdxprojection(s, b);   # Projection
io.print(res);                        # Print result
```

To run all these tests:

```bash
Mtest.py <absolute-path-to-this-repo>/monetdb5/modules/kernel/Tests
```

## TPC-H Benchmark

Currently we support TPC-H query 1 and 6.

### TPC-H Data

We use `tpch-tools`[5] to generate and load TPC-H data (SF 10) into MonetDB:

```bash
git clone https://github.com/eyalroz/tpch-tools.git --recursive

cd tpch-tools

# -h, help info
# -s, scale factor
# -d, database name
# -f, datafarm path
# -D, tpch data path, if already generated
# -v, verbose
# -k, keep data
# -G, use generated tpch data
# -P, daemon port

# one for naive
scripts/setup-tpch-db -s 10 -d tpch-naive-10 -f <path-to-dbfarm> -v -k

# one for bindex
scripts/setup-tpch-db -s 10 -d tpch-bindex-10 -f <path-to-dbfarm> -v -k -G

# 2 databases running with state `R`
monetdb status
```

### Query Execution

We set `sql_optimizer` to `sequential_pipe`[2] and use `EXPLAIN <SQL>` to get MAL instructions (`bench/mal/naive`). Then we modified it to work with Bindex API (`bench/mal/bindex`).

To run the benchmark:

```bash
cd bench

# run benchmark
./run.sh

# show result
./show.sh
```

## References

1. Linwei Li, Kai Zhang, Jiading Guo, Wen He, Zhenying He, Yinan Jing, Weili Han, X. Sean Wang. BinDex: A Two-Layered Index for Fast and Robust Scans. Proceedings of the 39th ACM International Conference on Management of Data (SIGMOD), Portland OR, USA, June 14-19, 2020.
2. https://www.monetdb.org/Documentation/mserver5-man-page
3. https://www.monetdb.org/Documentation/Manuals/MonetDB/MALreference
4. https://www.monetdb.org/Documentation/Manuals/MonetDB/Kernel/Modules/BAT
5. https://github.com/eyalroz/tpch-tools
