(*Addition du vecteur b au vecteur a*)

let add a b =
let n = vect_length a in
for k = 0 to n-1 do
a.(k) <- a.(k) +. b.(k)
done;;

(*Multiplication du vecteur a par lambda*)

let mult_vecteur lambda a =
let n = vect_length a - 1 in
  for k = 0 to n do
    a.(k) <- lambda *. a.(k)
  done
;;

(*Produit scalaire canonique de a et b*)

let scal a b =
let N = vect_length b in
let s = ref 0. in
for k = 0 to N-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(*Multiplication matricielle*)

let mult T x =
let p = vect_length T - 1 and n = vect_length T.(0) - 1 in
let v = make_vect p 0. in
  for k = 0 to n do
    let c = ref 0. in
      for i= 0 to p do
        c := !c +. T.(i).(k) *. x.(i)
      done ;
    v.(k) <- !c ;
  done ;
v ;;

(*e_i*)

let base T i =
let n = vect_length T in
let v = make_vect n 0. in
v.(i) <- 1. ; v ;;

(*(Tx|Te_i)*)

let B T x i = scal (mult T x) T.(i) ;;

(*q(x)*)

let q T x = scal (mult T x) (mult T x) ;;

(*minimum de B x e_i pour i*)

let miniB T x =
let p = vect_length T -1 in
let j = ref 0 in
  for k = 0 to p do
    if B T x k < B T x !j then j:=k
  done ;
!j ;;

(*Teste si un vecteur est solution*)

let test_solution T x =
let p = vect_length T - 1 in
let bool = ref true in
for k = 0 to p do
bool := (!bool)&&((B T x k) > 0.)
done ; !bool ;;

exception classifiable ;;

(*Algorithme de Von Neumann. Renvoie la nieme itération de l'algorithme et ||Tx_n||*)

let algorithme_von_neumann T n =
  let rec aux T x n = match n with
  |0 -> (x, sqrt(q T x))
  |n -> let j = miniB T x in
      let a = (q T x) -. (B T x j) in
      let b = (q T x) +. sqrt(scal T.(j) T.(j)) -. 2. *. (B T x j) in
      let lambda = a /. b in
      let v = base T j in
      mult_vecteur (1. -. lambda) x ;
      mult_vecteur lambda v ;
      add x v ;
      if test_solution T x then raise classifiable else aux T x (n-1)
in let v = base T 1 in
aux T v n ;;





