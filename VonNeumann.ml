open LinearAlgebra;;

(*minimum de bigB x e_i pour i*)

let minibigB bigT x =
let p = Array.length bigT -1 in
let j = ref 0 in
  for k = 0 to p do
    if bigB bigT x k < bigB bigT x !j then j:=k
  done ;
!j ;;

(*bigTeste si un vecteur est solution*)

let test_solution bigT x =
let p = Array.length bigT - 1 in
let bool = ref true in
for k = 0 to p do
bool := (!bool)&&((bigB bigT x k) > 0.)
done ; !bool ;;

(*Algorithme de Von Neumann. Renvoie la nieme itération de l'algorithme et ||bigTx_n||*)

let algorithme_von_neumann bigT n =
  let rec aux bigT x n = match n with
  |0 -> (x, sqrt(q bigT x))
  |n -> let j = minibigB bigT x in
      let a = (q bigT x) -. (bigB bigT x j) in
      let b = (q bigT x) +. sqrt(scal bigT.(j) bigT.(j)) -. 2. *. (bigB bigT x j) in
      let lambda = a /. b in
      let v = base bigT j in
      scal_mult (1. -. lambda) x ;
      scal_mult lambda v ;
      add x v ;
      if test_solution bigT x then raise (Failure "classifiable") else aux bigT x (n-1)
in let v = base bigT 1 in
aux bigT v n ;;
