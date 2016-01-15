(*copy a vector*)
let copy t =
let n = Array.length t in
let u = Array.make n t.(0) in
  for k = 0 to n-1 do
    u.(k) <- t.(k)
  done ;
u ;;

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
!result;;

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
let p = min (Array.length tt - 1) (Array.length x - 1) and n = Array.length tt.(0) - 1 in
let result = Array.make p 0. in
  for kk = 0 to n do
    let c = ref 0. in
      for ii = 0 to p do
        c := !c +. tt.(ii).(kk) *. x.(ii)
      done ;
    result.(kk) <- !c ;
  done ;
result ;;

(*making i-th canonical base vector of same size as the columns of tt, with a 1 at i-th index and 0 elsewhere*)
let base tt i =
	let n = Array.length tt.(0) in
let v = Array.make n 0. in
v.(i) <- 1. ; v ;;

(*bilinear product of x and e_i (i-th base vector) based on tt: (tt x|tt e_i)*)

let bilin_prod tt x i = scal (matrix_mult tt x) tt.(i) ;;

(*quadratic bilinear product of tt and x: (tt x | tt x)*)

let quad_prod tt x =
	let vect = matrix_mult tt x in
  scal vect vect ;;