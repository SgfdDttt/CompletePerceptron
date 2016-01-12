(*Copie conforme d'un tableau*)

let copie t =
let n = vect_length t in
let u = make_vect n t.(0) in
  for k = 0 to n-1 do
    u.(k) <- t.(k)
  done ;
u ;;

(*Addition en place du vecteur b au vecteur a*)
let add a b =
let n = vect_length a in
for k = 0 to n-1 do
a.(k) <- a.(k) +. b.(k)
done;;

(*Produit scalaire canonique de a et b*)

let scal a b =
let N = vect_length b in
let s = ref 0. in
for k = 0 to N-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(*Normalisation du vecteur x*)

let normalisation x =
let N = vect_length x in
let s = ref 0. in
for k = 0 to N-1 do
s:= !s +. (x.(k))*.(x.(k))
done;
s := sqrt (!s) ;
for k = 0 to N-1 do
x.(k) <- x.(k)/.(!s)
done ;;

(*Teste si u est le vecteur nul*)

let test_nul u =
let n = vect_length u - 1 in
let bool = ref true in
for k = 0 to n do
bool := (!bool)&&(u.(k) = 0.)
done ;
!bool;;

(*T.(k) désigne la kieme colonne de T. Recherche du plus petit produit scalaire*)

let miniscal T w =
let p = vect_length T - 1 in
let k = ref 0 in
for i = 0 to p do
if scal T.(i) w < scal T.(!k) w then k:=i done;
!k ;;

(*Marge d'un vecteur*)

let marge T w =
let j = miniscal T w and N = sqrt(scal w w) in (scal T.(j) w)/.N ;;

exception inclassifiable ;;

(*Cet algorithme renvoie le résultat de n itérations du perceptron avec la marge obtenue, en gardant en mémoire la meilleure solution rencontrée jusque là. Si le vecteur solution passe par 0, c'est inclassifiable d'après Gordan.*)

let perceptron_avancé T n =
  let rec aux y a T x n = match n with
  |0 -> (y, a)
  |n ->
      let j = miniscal T x in
      add x T.(j) ;
      if test_nul x
      then raise inclassifiable
      else
      begin
      if a < marge T x
      then let z = copie x in normalisation z ; aux z (marge T z) T x (n-1) else aux y a T x (n-1)
      end ;
in let N = vect_length T.(0) in
let w = make_vect N 0. in
w.(0) <- 1. ; aux w (marge T w) T w n ;;



