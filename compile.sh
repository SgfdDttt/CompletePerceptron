mkdir egs
cd src/
ocamlopt -o ../egs/setosa_vs_virginica LinearAlgebra.ml AdvancedPerceptron.ml VonNeumann.ml PalRujanPerceptron.ml SetosaVsVirginica.ml 
ocamlopt -o ../egs/setosa_vs_versicolor LinearAlgebra.ml AdvancedPerceptron.ml VonNeumann.ml PalRujanPerceptron.ml SetosaVsVersicolor.ml 
ocamlopt -o ../egs/virginica_vs_versicolor LinearAlgebra.ml AdvancedPerceptron.ml VonNeumann.ml PalRujanPerceptron.ml VirginicaVsVersicolor.ml 

