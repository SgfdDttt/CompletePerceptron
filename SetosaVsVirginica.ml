open LinearAlgebra
open AdvancedPerceptron
open VonNeumann
open PalRujanPerceptron

(*The data: Setosa and Virginica, from UCI Machine Learning Repository (http://archive.ics.uci.edu/ml/datasets/Iris).*)
(*The coefficients 1 to 4 are the features, the 5-th coefficient is the label.*)
let data = [|
[|5.1;3.5;1.4;0.2;1.|];
[|4.9;3.0;1.4;0.2;1.|];
[|4.7;3.2;1.3;0.2;1.|];
[|4.6;3.1;1.5;0.2;1.|];
[|5.0;3.6;1.4;0.2;1.|];
[|5.4;3.9;1.7;0.4;1.|];
[|4.6;3.4;1.4;0.3;1.|];
[|5.0;3.4;1.5;0.2;1.|];
[|4.4;2.9;1.4;0.2;1.|];
[|4.9;3.1;1.5;0.1;1.|];
[|5.4;3.7;1.5;0.2;1.|];
[|4.8;3.4;1.6;0.2;1.|];
[|4.8;3.0;1.4;0.1;1.|];
[|4.3;3.0;1.1;0.1;1.|];
[|5.8;4.0;1.2;0.2;1.|];
[|5.7;4.4;1.5;0.4;1.|];
[|5.4;3.9;1.3;0.4;1.|];
[|5.1;3.5;1.4;0.3;1.|];
[|5.7;3.8;1.7;0.3;1.|];
[|5.1;3.8;1.5;0.3;1.|];
[|5.4;3.4;1.7;0.2;1.|];
[|5.1;3.7;1.5;0.4;1.|];
[|4.6;3.6;1.0;0.2;1.|];
[|5.1;3.3;1.7;0.5;1.|];
[|4.8;3.4;1.9;0.2;1.|];
[|5.0;3.0;1.6;0.2;1.|];
[|5.0;3.4;1.6;0.4;1.|];
[|5.2;3.5;1.5;0.2;1.|];
[|5.2;3.4;1.4;0.2;1.|];
[|4.7;3.2;1.6;0.2;1.|];
[|4.8;3.1;1.6;0.2;1.|];
[|5.4;3.4;1.5;0.4;1.|];
[|5.2;4.1;1.5;0.1;1.|];
[|5.5;4.2;1.4;0.2;1.|];
[|4.9;3.1;1.5;0.1;1.|];
[|5.0;3.2;1.2;0.2;1.|];
[|5.5;3.5;1.3;0.2;1.|];
[|4.9;3.1;1.5;0.1;1.|];
[|4.4;3.0;1.3;0.2;1.|];
[|5.1;3.4;1.5;0.2;1.|];
[|5.0;3.5;1.3;0.3;1.|];
[|4.5;2.3;1.3;0.3;1.|];
[|4.4;3.2;1.3;0.2;1.|];
[|5.0;3.5;1.6;0.6;1.|];
[|5.1;3.8;1.9;0.4;1.|];
[|4.8;3.0;1.4;0.3;1.|];
[|5.1;3.8;1.6;0.2;1.|];
[|4.6;3.2;1.4;0.2;1.|];
[|5.3;3.7;1.5;0.2;1.|];
[|5.0;3.3;1.4;0.2;1.|];
[|6.3;3.3;6.0;2.5;-1.|];
[|5.8;2.7;5.1;1.9;-1.|];
[|7.1;3.0;5.9;2.1;-1.|];
[|6.3;2.9;5.6;1.8;-1.|];
[|6.5;3.0;5.8;2.2;-1.|];
[|7.6;3.0;6.6;2.1;-1.|];
[|4.9;2.5;4.5;1.7;-1.|];
[|7.3;2.9;6.3;1.8;-1.|];
[|6.7;2.5;5.8;1.8;-1.|];
[|7.2;3.6;6.1;2.5;-1.|];
[|6.5;3.2;5.1;2.0;-1.|];
[|6.4;2.7;5.3;1.9;-1.|];
[|6.8;3.0;5.5;2.1;-1.|];
[|5.7;2.5;5.0;2.0;-1.|];
[|5.8;2.8;5.1;2.4;-1.|];
[|6.4;3.2;5.3;2.3;-1.|];
[|6.5;3.0;5.5;1.8;-1.|];
[|7.7;3.8;6.7;2.2;-1.|];
[|7.7;2.6;6.9;2.3;-1.|];
[|6.0;2.2;5.0;1.5;-1.|];
[|6.9;3.2;5.7;2.3;-1.|];
[|5.6;2.8;4.9;2.0;-1.|];
[|7.7;2.8;6.7;2.0;-1.|];
[|6.3;2.7;4.9;1.8;-1.|];
[|6.7;3.3;5.7;2.1;-1.|];
[|7.2;3.2;6.0;1.8;-1.|];
[|6.2;2.8;4.8;1.8;-1.|];
[|6.1;3.0;4.9;1.8;-1.|];
[|6.4;2.8;5.6;2.1;-1.|];
[|7.2;3.0;5.8;1.6;-1.|];
[|7.4;2.8;6.1;1.9;-1.|];
[|7.9;3.8;6.4;2.0;-1.|];
[|6.4;2.8;5.6;2.2;-1.|];
[|6.3;2.8;5.1;1.5;-1.|];
[|6.1;2.6;5.6;1.4;-1.|];
[|7.7;3.0;6.1;2.3;-1.|];
[|6.3;3.4;5.6;2.4;-1.|];
[|6.4;3.1;5.5;1.8;-1.|];
[|6.0;3.0;4.8;1.8;-1.|];
[|6.9;3.1;5.4;2.1;-1.|];
[|6.7;3.1;5.6;2.4;-1.|];
[|6.9;3.1;5.1;2.3;-1.|];
[|5.8;2.7;5.1;1.9;-1.|];
[|6.8;3.2;5.9;2.3;-1.|];
[|6.7;3.3;5.7;2.5;-1.|];
[|6.7;3.0;5.2;2.3;-1.|];
[|6.3;2.5;5.0;1.9;-1.|];
[|6.5;3.0;5.2;2.0;-1.|];
[|6.2;3.4;5.4;2.3;-1.|];
[|5.9;3.0;5.1;1.8;-1.|]
|] ;;

data_preprocessing data ;;

print_string "===== ADVANCED PERCEPTRON =====\n";

let (weights, margin) = advanced_perceptron data 1 in
print_string "advanced perceptron, 1 iteration:\n";
print_weight_and_margin weights margin "    "; print_string "\n";

let (weights, margin) = advanced_perceptron data 10 in
print_string "advanced perceptron, 10 iterations:\n";
print_weight_and_margin weights margin "    "; print_string "\n";

let (weights, margin) = advanced_perceptron data 100 in
print_string "advanced perceptron, 100 iterations:\n";
print_weight_and_margin weights margin "    "; print_string "\n";

print_string "===== VON NEUMANN ALGORITHM =====\n";

let (weights, margin) = von_neumann_algorithm data 1 in
print_string "von_neumann_algorithm, 1 iteration:\n";
print_coefficients_and_module weights margin "    "; print_string "\n";

let (weights, margin) = von_neumann_algorithm data 10 in
print_string "von_neumann_algorithm, 10 iterations:\n";
print_coefficients_and_module weights margin "    "; print_string "\n";


let (weights, margin) = von_neumann_algorithm data 100 in
print_string "von_neumann_algorithm, 100 iterations:\n";
print_coefficients_and_module weights margin "    "; print_string "\n";

print_string "===== PAL RUJAN ALGORITHM =====\n";

let (weights, margin, iterations) = pal_rujan_algorithm data in
print_string "finished in "; print_int iterations; print_string " iterations.\n";
if (margin = 0.) then
	begin
	  print_string "no solution to perceptron problem.\n";
	  print_coefficients_and_module weights margin "    "; print_string "\n";
	end
else
	begin
	  print_string "solution to perceptron problem:\n";
	  print_weight_and_margin weights margin "    "; print_string "\n";
	end;
;;