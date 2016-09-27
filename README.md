### Perceptron
An OCaml implementation of the classic Perceptron algorithm and its counterpart the Von Neumann algorithm, as well as Pál Ruján's Perceptron with Maximal Stability (https://hal.archives-ouvertes.fr/jpa-00246719/document). The latter combines the first 2 algorithms. If the data is linearly separable, it finds the linear classifier with the highest margin; if the data isn't separable, it can prove it, by finding a solution to an alternate linear problem.

### Summaries
**PerceptronThesis.pdf** is a paper about the code, with mathematical notation, theorems and proofs. It's all in French, as I wrote it for a competitive exam I took in France.

**PerceptronGist.pdf** is a summary in English containing the algorithms and important theorems.

### Running the algorithms
This repository contains a demo of the 3 algorithms on the Iris Data Set from the UCI Machine Learning repository (https://archive.ics.uci.edu/ml/datasets/Iris/). You need Ocaml to compile the .ml files.

1. Run compile.sh.

2. Run egs/xxx_vs_yyy to run all 3 algorithms on the given set of data, attempting to separate iris xxx from iris yyy. Results will be printed to the console.

To run the algorithms on your own data, have a look at the XxxVsYyy.ml files and change the appropriate parts. You will need to write an interface to read your data.
