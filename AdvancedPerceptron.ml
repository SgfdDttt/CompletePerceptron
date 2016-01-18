open LinearAlgebra;;


(*runs the perceptron for n interations on matrix tt, keeping in memory the best solution found so far and its margin.*)
(*the data are the columns of tt, with an extra entry at the bottom equal to 1, and the whole vector multiplied by 1 or -1 depending on the class it belongs to*)
(*if at some point, the weight vector is 0, then there is no solution (Gordan Theorem)*)

let advanced_perceptron tt n =
  let rec aux best_weight_vector best_margin tt current_weight_vector remaining_iterations = match remaining_iterations with
    |0 -> begin
			      if all_zero current_weight_vector then
							begin
							  print_string "NO SOLUTION;";
							end;
							normalisation best_weight_vector; (best_weight_vector, best_margin)
          end;    
		|n -> let j = mini_scal tt current_weight_vector in
          add current_weight_vector tt.(j) ;       
              if best_margin < margin tt current_weight_vector then
								let new_best = copy current_weight_vector in
			          normalisation current_weight_vector ;
	              aux new_best (margin tt new_best) tt current_weight_vector (remaining_iterations-1);
		          else
							  aux best_weight_vector best_margin tt current_weight_vector (remaining_iterations-1);
	in let w = Array.make (Array.length tt.(0)) 1. in
  normalisation w;
  aux w (margin tt w) tt w n;
;;
