open Printf;;

(*copy a vector*)
let copy t =
let n = Array.length t in
let u = Array.make n t.(0) in
  for k = 0 to n-1 do
    u.(k) <- t.(k)
  done ;
u ;;

(*print an array - that's not algebra but still useful*) 
let print_array a =
  let print_with_semicolon x = print_float x; print_string "; ";
	in Array.iter (print_with_semicolon) a;
;;

(*print weights and margin - that's not algebra but still useful*) 
let print_weight_and_margin weights margin prefix =
	print_string prefix; print_string "weights:\n    ";
	print_string prefix; print_array weights; print_string "\n";
	print_string prefix; print_string "margin:\n    ";
	print_string prefix; print_float margin;
;;

(*print convex_coefficients and module - that's not algebra but still useful*) 
let print_coefficients_and_module weights margin prefix =
	print_string prefix; print_string "convex coefficients:\n    ";
	print_string prefix; print_array weights; print_string "\n";
	print_string prefix; print_string "module:\n    ";
	print_string prefix; print_float margin;
;;

(*add vector b to vector a, in place*)
let add a b =
for k = 0 to Array.length a - 1 do
a.(k) <- a.(k) +. b.(k)
done;;

(*multiplication of vector a by lambda, in place*)
let scal_mult lambda a =
let n = Array.length a - 1 in
  for k = 0 to n do
    a.(k) <- lambda *. a.(k)
  done
;;

(*scalar product of vectors a and b*)
let scal a b =
  let length = min (Array.length b) (Array.length a) in
  let result = ref 0. in
  for k = 0 to length-1 do
    result := !result +. (a.(k))*.(b.(k))
  done;
  !result;
;;

(*normalisation of vector x, in place*)
let normalisation x =
let length = Array.length x in
let module_x = ref 0. in
for k = 0 to length - 1 do
module_x:= !module_x +. (x.(k))*.(x.(k))
done;
module_x := sqrt (!module_x) ;
for k = 0 to length - 1 do
x.(k) <- x.(k)/.(!module_x)
done ;;

(*checks to see if all the coefficients in a vector are 0*)
let all_zero m =
let n = Array.length m and bool = ref true in
for k=0 to n-1 do
if m.(k) <> 0. then bool := false
done ;
!bool;;

(*looking for the index of the column of tt that has the smallest scalar product with w. tt.(k) is the k-th column of t*)
let mini_scal tt w =
let p = Array.length tt - 1 in
let mini_index = ref 0 in
for ii = 0 to p do
if scal tt.(ii) w < scal tt.(!mini_index) w then mini_index:=ii done;
!mini_index ;;

(*margin of a weight vector*)
let margin tt w =
let j = mini_scal tt w and module_w = sqrt(scal w w) in (scal tt.(j) w)/.module_w ;;

(*matrix multiplication*)
let matrix_mult tt x =
  let p = min (Array.length tt) (Array.length x) in
	match p with
	| 0 -> [||];
	|_ -> let n = Array.length tt.(0) in
        let result = Array.make n 0. in
        for kk = 0 to n-1 do
          let c = ref 0. in
          for ii = 0 to p-1 do
            c := !c +. tt.(ii).(kk) *. x.(ii)
          done ;
          result.(kk) <- !c ;
        done ;
        result;
;;

(*making i-th canonical base vector of length the number of columns that tt has, with a 1 at i-th index and 0 elsewhere*)
let base tt i =
	let n = Array.length tt in
  let v = Array.make n 0. in
  v.(i) <- 1. ;
  v;
;;

(*bilinear product of x and e_i (i-th base vector) based on tt: (tt x|tt e_i)*)

let bilin_prod tt x i =
	let mat = matrix_mult tt x in
  scal mat tt.(i);
;;

(*quadratic bilinear product of tt and x: (tt x | tt x)*)

let quad_prod tt x =
	let vect = matrix_mult tt x in
  scal vect vect ;;

(*Gaussian pivot, solving matrix equation A*X = B*)

(*elementary operation of Gaussian pivot*)
(*substracting b times the line i from the line j*)
let elementary_operation aa xx i j b =
  let p = Array.length aa in
    for k = 0 to p-1 do
      aa.(k).(j) <- aa.(k).(j) -. b*.aa.(k).(i)
    done ;
  xx.(j) <- xx.(j) -. b*.xx.(i) ;
;;

(*Exchanging 2 lines in the system*)
let element_switch aa i j k =
  let mem = aa.(k).(i) in
  aa.(k).(i) <- aa.(k).(j) ;
  aa.(k).(j) <- mem ;
;;

let line_switch aa x i j =
  let mem = x.(i) in
  x.(i) <- x.(j) ;
  x.(j) <- mem ;
  let p = Array.length aa in
    for k = 0 to p-1 do
      element_switch aa i j k
    done ;
;;

(*Solving A*X = B when A is a upper triangular matrix. A assumed to be square.*)
(* If there is no solution, it returns false and some vector. Otherwise, true and the solution.*)

let rec solve_triangular aa bb =
	match Array.length aa with
	|0 -> (true, [||]);
  |1 -> if aa.(0).(0) = 0. then
          begin
            if bb.(0) <> 0. then (false, bb) else (true, [|0.|])
          end
        else (true, [| bb.(0) /. aa.(0).(0)|])
  |n when (aa.(n-1).(n-1) = 0.)&&(bb.(n-1) <> 0.) -> (false, bb)
  |n when aa.(n-1).(n-1) = 0. -> let bb_1 = Array.make (n-1) 0. in
                                 for k = 0 to (n-2) do
                                   bb_1.(k) <- bb.(k)
                                 done ;
                                 let aa_1 = Array.make_matrix (n-1) (n-1) 0. in
                                 for a = 0 to n-2 do
                                   for b = 0 to n-2 do
                                     aa_1.(a).(b) <- aa.(a).(b)
                                   done ;
                                 done ;
                                 let (bool, xx_1) = solve_triangular aa_1 bb_1 in 
                                 let xx = Array.make n 0. in
                                 for k = 0 to n-2 do
                                   xx.(k) <- xx_1.(k)
                                 done ;
                                   (bool, xx)
  |n -> let y = bb.(n-1) /. aa.(n-1).(n-1) in
        let bb_1 = Array.make (n-1) 0. in
        for k = 0 to (n-2) do
          bb_1.(k) <- bb.(k) -. y *. aa.(n-1).(k)
        done ;
        let aa_1 = Array.make_matrix (n-1) (n-1) 0. in
        for a = 0 to n-2 do
          for b = 0 to n-2 do
            aa_1.(a).(b) <- aa.(a).(b)
          done ;
        done ;
        let bool, xx_1 = solve_triangular aa_1 bb_1 in 
        let xx = Array.make n 0. in
        for k = 0 to n-2 do
          xx.(k) <- xx_1.(k)
        done ;
				xx.(n-1) <- y;
        (bool, xx)
;;

(*Solving A*X = B in the general case. A assumed to be square.*)
(* If there is no solution, it returns false and some vector. Otherwise, true and the solution.*)
let solve_gaussian_pivot aa bb =
  let rec aux matrix_a matrix_b n = match n with
    |n when n = Array.length matrix_a -> solve_triangular matrix_a matrix_b
    |n -> let p = Array.length matrix_a in
          let pivot = ref p in
          for i = p-1 downto n do
            if matrix_a.(n).(i) <> 0. then pivot := i
          done ;
          match !pivot with
            |x when x=p -> aux matrix_a matrix_b (n+1)
            |_ -> line_switch matrix_a matrix_b !pivot n ;
                  for j = n+1 to p-1 do
                    elementary_operation matrix_a matrix_b n j (matrix_a.(n).(j)/.matrix_a.(n).(n))
                  done ;
                  aux matrix_a matrix_b (n+1)
  in aux aa bb 0
;;

(*Preprocessing the data*)

(*changes c_i into y_i*c_i*)
let g v =
  let n = Array.length v in
    for k=0 to n-2 do
	    v.(k) <- v.(k)*.v.(n-1);
    done;
;;

(*changes all c_i into y_i*c_i*)
let treat tt =
  let n = Array.length tt in
    for k = 0 to n-1 do
      g tt.(k)
    done;
;;

let vector_normalisation v =
  let n = Array.length v in
  let denominator = ref 0. in
  for k = 0 to n-1 do
    denominator := !denominator +. v.(k) *. v.(k)
  done ;
  denominator := sqrt(!denominator) ;
  for k = 0 to n-1 do
    v.(k) <- v.(k) /. !denominator
  done ;
;;

let data_normalisation tt =
let p = Array.length tt in
let n = Array.length tt.(0) in
let normalized = Array.make_matrix p n 0. in
  for k = 0 to p-1 do
    vector_normalisation tt.(k);
  done;
normalized;
;;

let data_preprocessing tt =
	treat tt;
	data_normalisation tt;
;;