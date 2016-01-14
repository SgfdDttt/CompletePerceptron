(*Addition du vecteur b au vecteur a*)

let add a b =
let n = Array.length a in
for k = 0 to n-1 do
a.(k) <- a.(k) +. b.(k)
done;;

(*Multiplication du vecteur a par lambda*)

let mult_vecteur lambda a =
let n = Array.length a - 1 in
  for k = 0 to n do
    a.(k) <- lambda *. a.(k)
  done
;;

(*Produit scalaire canonique de a et b*)

let scal a b =
let bigN = Array.length b in
let s = ref 0. in
for k = 0 to bigN-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(*Multiplication matricielle*)

let mult bigT x =
let p = Array.length bigT - 1 and n = Array.length bigT.(0) - 1 in
let v = Array.make p 0. in
  for k = 0 to n do
    let c = ref 0. in
      for i= 0 to p do
        c := !c +. bigT.(i).(k) *. x.(i)
      done ;
    v.(k) <- !c ;
  done ;
v ;;

(*e_i*)

let base bigT i =
let n = Array.length bigT in
let v = Array.make n 0. in
v.(i) <- 1. ; v ;;

(*(bigTx|bigTe_i)*)

let bigB bigT x i = scal (mult bigT x) bigT.(i) ;;

(*q(x)*)

let q bigT x = scal (mult bigT x) (mult bigT x) ;;

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
      mult_vecteur (1. -. lambda) x ;
      mult_vecteur lambda v ;
      add x v ;
      if test_solution bigT x then raise (Failure "inclassifiable") else aux bigT x (n-1)
in let v = base bigT 1 in
aux bigT v n ;;





