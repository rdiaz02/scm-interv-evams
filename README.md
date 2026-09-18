R code for [A structural causal framework for interventions on evolutionary accumulation models](http://arxiv.org/abs/2606.12597), by Diaz-Uriarte, Arroyo, and Johnston.


To run the code and all tests you'll need to install:
- [evamtools](https://github.com/rdiaz02/EvAM-Tools). Because of its dependencies, that package is not available from CRAN or BioConductor. Installation is detailed in  [Installing and running EvAM Tools](https://github.com/rdiaz02/EvAM-Tools#installing-and-running). If you are on Windows, the fastest route might be to use a Docker image (explained also in that link). If you are on macOS (at least macOS on Apple silicon), the only system dependency is graphviz from Homebrew, which Rgraphviz and then OncoSimulR need.
- [hyperhmm](https://github.com/StochasticBiology/hyperhmm) package. From R, you can do `remotes::install_github("StochasticBiology/hyperhmm")` (if you do not have the `remotes` package installed, install it first by doing `install.packages("remotes")`).
- You will also need to install the following packages, available from CRAN: `uuid`, `expm`, `markovchain` and `testthat`. From R you can type `install.packages(c("uuid", "expm", "markovchain", "testthat"))`.


[RhpcBLASctl](https://cran.r-project.org/web/packages/RhpcBLASctl/index.html) is a required package (it is also imported by `evamtools`). The tests use `mclapply` over all available cores, and several operations (markovchain, the Armadillo/BLAS solves, `genots_from_trm`) call into the BLAS, which with many builds is itself multithreaded. Without limiting the BLAS to one thread, each forked worker would start its own pool of threads and could lead to the machine being heavily oversubscribed. (With the single-threaded reference BLAS that R uses by default on macOS and Windows, none of this is needed. On a macOS R linked against Apple's Accelerate framework the calls to `RhpcBLASctl`'s functions are expected to have no effect so you might want to adjust the number of threads yourself before launching the tests; this has not been tested here.)

## Running the tests

The recommended way is to change to the test directory (`code_intervention/tests`) and do, from a terminal, `R --vanilla -f run-all-tests.R &> run-all-tests.Rout`. The tests take between 35 and 45 minutes, depending on hardware.
