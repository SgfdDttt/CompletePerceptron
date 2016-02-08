### Perceptron
An OCaml implementation of the classic Perceptron algorithm and its counterpart the Von Neumann algorithm, as well as Pál Ruján's Perceptron with Maximal Stability (https://hal.archives-ouvertes.fr/jpa-00246719/document).

### Summaries
**PerceptronThesis.pdf** contains the full version of the paper about the code, with mathematical notation, theorems and proofs. It's all in French, as I wrote it for a competitive exam I took in France.

**PerceptronGist.pdf** is a summary in English containing the algorithms and important theorems.

### Running the algorithms
You need ocaml to run compile the .ml files.

Run the bash script compile.sh

Run egs/xxx_vs_yyy to run all 3 algorithms on the given set of data, attempting to separate xxx from yyy. Results will be printed to the console.

To run the algorithms on your own data, have a look at the XxxVsYyy.ml files and change the appropriate parts. You might need to write an interface to read your data.
