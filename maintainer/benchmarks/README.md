# ESPResSo Benchmark suite

This folder contains the benchmark suite of ESPResSo implemented using the [Reframe](https://reframe-hpc.readthedocs.io/en/stable/#) framework. The benchmark suite is responsible for pulling and building ESPResSo (alternatively a prebuilt ESPResSo can be used), running a set of benchmark tests and summarising, as well as visualising the results.

## Usage

The benchmark suite is meant to be executed from two different systems:  
(1) Local execution from different workstations.  
(2) Execution on the [ant](https://pages.icp.uni-stuttgart.de/pages/hpc/_pages/_hpc_facilities/ant.html) HPC-cluster.  

Reframe will detect the current system and automatically choose the correct mode of execution. To run the benchmark suite, execute the script [suite.sh](suite.sh). The script will execute the reframe benchmark suite and write all artefacts into a folder with the name `<PREFIX>_dd_mm_yyyy_<n>`, with a specified prefix, dd\_mm\_yyyy the current date and \<n\> counting the number of runs on the same day.


## Structure

| File | Description |
| ---- | ----------- |
| [ferrofluid.py](ferrofluid.py) / [lb.py](lb.py) / [lj.py](lj.py) / [mc_acid_base_reservoir.py](mc_acid_base_reservoir.py) / [p3m.py](p3m.py) | Implemenentation of benchmark tests for various systems |
| [benchmarks.py](benchmarks.py) | Helper functions used in the benchmark files |
| [espresso_benchmarks.py](espresso_benchmarks.py) | Reframe testsuite, defining how different test configurations are wired together |
| [suite.sh](suite.sh) | Script to execute the reframe benchmark suite, exposing different options |
| [reframe_config.py](reframe_config.py) | Configuration of reframe to our system architecture |
| [benchmark_utils.py](benchmark_utils.py) | Definition of all benchmark test parameterisations |
| [plot_benchmarks.py](plot_benchmarks.py) | Helper script to generate an interactive timeline of past benchmark results |
| [clean_stage.py](clean_stage.py) | Helper script to clean up left-over espresso builds that are created during the benchmark execution |
| [benchmark_history.py](benchmark_history.py) | Helper script to generate a html file from benchmark run results, including the interactive figure from plot_benchmarks.py |
| [publish_pages.sh](publish_pages.sh) | Helper script to automate pushing the generated summary html to a gitlab pages repository (Meant to be executed within a cron-job) |
