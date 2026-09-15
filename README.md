R code for [A structural causal framework for interventions on evolutionary accumulation models](http://arxiv.org/abs/2606.12597), by Diaz-Uriarte, Arroyo, and Johnston.


To run the code you'll need to install the [evamtools](https://github.com/rdiaz02/EvAM-Tools) R package. Because of its dependencies, that package is not available from CRAN or BioConductor. Installation is detailed in  [Installing and running EvAM Tools](https://github.com/rdiaz02/EvAM-Tools#installing-and-running). If you are Windows, the fastest, most expedite route might be to use a Docker image (explained also in that link). If you are on macOS (at least macOS on Apple silicon), the only system dependency is graphviz from Homebrew, which Rgraphviz and then OncoSimulR need.

You will also need to install the [hyperhmm](https://github.com/StochasticBiology/hyperhmm) package. From R, you can do `remotes::install_github("StochasticBiology/hyperhmm")` (if you do not have the `remotes` package installed, install it first by doing `install.packages("remotes")`).
