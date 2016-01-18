open LinearAlgebra;;

(*minimum of (Bx | e_i) over i*)
let mini_index tt x =
  let p = Array.length tt in
  let j = ref 0 in
    for k = 0 to p-1 do
      if bilin_prod tt x k < bilin_prod tt x !j then j := k
    done ;
  !j;
;;

(*checks if a weight vector is a solution to the perceptron problem*)
let test_solution tt x =
  let p = Array.length tt in
  let bool = ref true in
  for k = 0 to p-1 do
    bool := (!bool)&&((scal tt.(k) x) > 0.)
  done ;
  !bool
;;

(*Von Neumann algorithm. Returns the n-th iteration and the minimum value attained.*)
let von_neumann_algorithm tt n =
  let rec iteration matrix current_solution it =
		match it with
    |0 -> begin
			      match test_solution matrix current_solution with
		          | true -> print_string "NO SOLUTION\n";
					    | false -> ();
					end;
				  (current_solution, sqrt(quad_prod matrix current_solution));
    |int when int > 0 -> let j = mini_index matrix current_solution in
										let a = (quad_prod matrix current_solution) -. (bilin_prod matrix current_solution j) in
										let b = (quad_prod matrix current_solution) +. sqrt(scal matrix.(j) matrix.(j)) -. 2. *. (bilin_prod matrix current_solution j) in
										let lambda = a /. b in
									  let column_j = base matrix j in
                    scal_mult (1. -. lambda) current_solution ;
                    scal_mult lambda column_j ;
                    add current_solution column_j ;
					           begin
                       match test_solution matrix current_solution with
					               |true -> iteration matrix current_solution (it-1);
					               |false -> iteration matrix current_solution (it-1);
					           end;
		|_ -> raise (Failure "von_neumann_algorithm:iteration: unexpected match case");
	in
	let v = base tt 1 in
	iteration tt v n;
;;


(*let mat = [| [|1.;1.;1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]|] in
let weights, margin = von_neumann_algorithm mat 1 in
print_string "1 iteration:\n";
print_coefficients_and_module weights margin "";
print_string"\n\n";
;;

let mat = [| [|1.;1.;1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]|] in
let weights, margin = von_neumann_algorithm mat 10 in
print_string "10 iterations:\n";
print_coefficients_and_module weights margin "";
print_string"\n\n";
;;

let mat = [| [|1.;1.;1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]; [|-1.; -1.; -1.|]|] in
let weights, margin = von_neumann_algorithm mat 100 in
print_string "100 iterations:\n";
print_coefficients_and_module weights margin "";
print_string"\n\n";
;;
*)