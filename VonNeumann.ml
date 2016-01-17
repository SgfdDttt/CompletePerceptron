open LinearAlgebra;;

(*minimum of (Bx | e_i) over i*)
let mini_index tt x =
let p = (Array.length tt) - 1 in
let j = ref 0 in
  for k = 0 to p do
    if bilin_prod tt x k < bilin_prod tt x !j then j := k
  done ;
!j ;;

(*checks if a weight vector is a solution to the perceptron problem*)
let test_solution tt x =
let p = Array.length tt - 1 in
let bool = ref true in
for k = 0 to p do
bool := (!bool)&&((scal tt.(k) x) > 0.)
done ; !bool ;;

(*Von Neumann algorithm. Returns the n-th iteration and the minimum value attained.*)
let algorithme_von_neumann tt n =
  let rec aux matrix current_solution n =
		match n with
  |0 -> (current_solution, sqrt(quad_prod matrix current_solution))
  |n -> let j = mini_index matrix current_solution in
      let a = (quad_prod matrix current_solution) -. (bilin_prod matrix current_solution j) in
      let b = (quad_prod matrix current_solution) +. sqrt(scal matrix.(j) matrix.(j)) -. 2. *. (bilin_prod matrix current_solution j) in
      let lambda = a /. b in
      let column_j = base matrix j in
      scal_mult (1. -. lambda) current_solution ;
      scal_mult lambda column_j ;
      add current_solution column_j ;
      if test_solution matrix current_solution then
				  raise (Failure "no solution")
				else
					aux matrix current_solution (n-1)
in let v = base tt 1 in
aux tt v n ;;