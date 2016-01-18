open LinearAlgebra;;

(*Pal Rujan algorithm
Matrices are float vect vect, vectors are float vect. M.(i) is the i-th column of M*)
(*The active set A_m is an array of length p+1. The slot of index 0 contains the cardinal of A_m. The slot i contains 1 if c_i belongs to A_m, and 0 otherwise.*)

(*Presence of an element in a list*)
let rec presence x l = match l with
| [] -> false
| t::q -> (x=t)||(presence x q) ;;

(*Within the columns of tt, looks for a vector whose scalar product with x is the smallest.*)
(*If that scalar product is smaller than alpha (the margin), it returns the number of that vector (and not the index of the column in tt), otherwise it returns 0*)
let smaller_scalar_product x alpha tt =
  let p = Array.length tt in
  let k = ref 0 in
    for i = 1 to p do
      if scal x tt.(i-1) < alpha
      then match !k with
        |0 -> k := i
        |n -> if scal x tt.(i-1) < scal x tt.(n-1) then k := i ;
    done ;
  !k ;;

(*Computing A_m.(0). This algorithm is useful in practice.*)
let updating_cardinal aa =
let n = Array.length aa in
let j = ref 0 in
  for k = 1 to n-1 do
    if aa.(k) = 1 then j := !j + 1
  done ;
aa.(0) <- !j ;;

(*Matrix containing c_i where i belongs to A_m, and matrix that is on the left side of the equation used to find the weights.*)
let helper_matrices a_m tt qq =
  let n = Array.length tt in
  let k = a_m.(0) in
  let p = Array.length qq in
  let zz = Array.make_matrix (k+1) (k+1) 0. in
    for i = 1 to k do
      zz.(i).(0) <- 1. ; zz.(0).(i) <- -1.
    done;
  let x = ref 1 in
  for a = 1 to p do
    if a_m.(a) = 1 then
    begin
    let y = ref 1 in
      for b = 1 to p do
        if a_m.(b) = 1 then
        begin
        zz.(!x).(!y) <- qq.(a-1).(b-1) ; y := !y + 1
        end ;
      done ;
    x := !x + 1 ;
   end ;
  done ;
  let ss = Array.make_matrix k n 0. in
  let x = ref 0 in
    for a = 1 to p do
      if a_m.(a) = 1 then
      begin
        ss.(!x) <- tt.(a-1) ; x := !x + 1 ;
      end ;
    done ;
  ss,zz;;

(*updating the active set : A_m <- A_m+1*)
let active_set a_m v j =
  let x = ref 0 in
  for i = 0 to Array.length v - 1 do
    if (v.(i) < v.(!x))&&(i <> j) then x := i
  done ;
  if v.(!x) > 0. then
		0
  else
    let m = ref (!x + 1) in let k = ref 1 in
    while !m > 0 do
      if a_m.(!k) = 1 then m := !m - 1 ;
        k := !k + 1 ;
      done ;
    a_m.(!k-1) <- 0 ; updating_cardinal a_m ;
!k-1 ;;

(*Take a vector out of A_m that is not vector number k*)
let take_vector_out a_m k =
let bool = ref true in
let x = ref 1 in
  while !bool do
    if (a_m.(!x) = 0)||(!x = k) then
		  x := !x + 1
	  else
		  bool := false
  done ;
a_m.(!x) <- 0 ; updating_cardinal a_m ;;

(*finds coefficients lambda_i >= 0 such that sum lambda_i c_i = 0 based on u, the solution of the characteristic equation*)
let convex_coefficients a_m u =
  let p = Array.length a_m in
  let v = Array.make (p-1) 0. in
  let x = ref 0 in
    for k = 1 to p-1 do
      match a_m.(k) with
        | 1 -> v.(k-1) <- u.(!x) ; x := !x + 1
        | 0 -> ()
		    | _ -> raise (Failure "convex_coefficients: unexpected value")
    done ;
v ;;

(*Margin of weight vector u on active set a_m*)
let margin_active_set a_m u tt =
  let number_lowest_scal_product = ref 1 in
  while a_m.(!number_lowest_scal_product) = 0 do
    if a_m.(!number_lowest_scal_product) = 0 then number_lowest_scal_product := !number_lowest_scal_product + 1;
	  if (!number_lowest_scal_product=Array.length a_m) then raise (Failure "empty active set");
  done ;
  for k = 1 to Array.length a_m - 1 do
    if (a_m.(k) = 1)&&(scal u tt.(k-1) <= scal u tt.(!number_lowest_scal_product-1))
    then number_lowest_scal_product := k
  done ;
scal u tt.(!number_lowest_scal_product-1) ;;

(*If there is a solution to the characteristic equation, returns*)
(*whether the characteristic equation has a solution,*)
(*the vecteur obtained by removing the first coefficient of the solution,*)
(*that first coefficient,*)
(*and s_m*)
let temporary_solution tt qq a_m =
  let ss,zz = helper_matrices a_m tt qq in
  let v = Array.make (a_m.(0)+1) 0. in
  v.(0) <- 1. ;
  let (bool, lambda) = solve_gaussian_pivot zz v in
  let u = Array.make a_m.(0) 0. in
  for k = 0 to a_m.(0)-1 do
    u.(k) <- lambda.(k+1)
  done ;
  bool, u, lambda.(0), ss
;;

(*updates the active set until all lambda_i are strictly positive*)
let rec update_active_set tt qq a_m j =
  match temporary_solution tt qq a_m with
    | bool, _, z, _ when (not bool)||(z < 0.) -> take_vector_out a_m j ; updating_cardinal a_m ; update_active_set tt qq a_m j;
    | true, u, _, _ -> begin match active_set a_m u j with
                         | 0 -> updating_cardinal a_m;
                         | k -> updating_cardinal a_m ; update_active_set tt qq a_m j;
											 end
		|_,_,_,_ -> raise (Failure "update_active_set: unexpected configuration");
;;

(*One iteration of the algorithm. Returns the active set.*)
let rec aux_pal_rujan tt qq a_m l number_iterations =
  if presence a_m l then
    false, a_m, number_iterations
	else
    let bb = copy a_m in
    let (bool, u, z, ss) = temporary_solution tt qq a_m in
    match z with
      | 0. -> true, a_m, (number_iterations + 1)
      | _ -> let w = matrix_mult ss u in
             normalisation w ;
             let c = margin_active_set a_m w tt in
             match smaller_scalar_product w c tt with
               |0 -> true, a_m, (number_iterations + 1)
               |k -> a_m.(k) <- 1 ; updating_cardinal a_m ; update_active_set tt qq a_m k ; updating_cardinal a_m ; aux_pal_rujan tt qq a_m (bb::l) (number_iterations+1)
;;

(*Pal Rujan algorithm. Returns*)
(*a vector of weights, which is either a solution to the perceptron problem, or to the alternate problem,*)
(*a positive real number, which is either the margin of the weight vector, or 0. (alternate problem),*)
(*the number of iterations before completion of the algorithm.*)
let pal_rujan_algorithm tt =
  let rec aux tt k acc = 
    let p = Array.length tt in
    let qq = Array.make_matrix p p 0. in
      for a = 0 to p-1 do
        for b = 0 to p-1 do
          qq.(a).(b) <- scal tt.(a) tt.(b)
        done ;
      done ;
    let aa = Array.make (p+1) 0 in
    aa.(0) <- 1 ; aa.(k) <- 1 ;
    let (bool, bb, n) = aux_pal_rujan tt qq aa [] 0 in
    match bool, temporary_solution tt qq bb with
    |false, (_, _, _, _) -> aux tt (k+1) (acc+n)
    |_, (_, u, 0., _) -> (convex_coefficients bb u), 0., (acc+n)
    |_, (_, u, _, ss) -> let w = matrix_mult ss u in normalisation w ;
                         match margin_active_set bb w tt with
                           |c when c < 0. -> aux tt (k+1) (acc+n)
                           |c -> (w, c, acc+n)
  in aux tt 1 0
;;