(*Copie conforme d'un tableau*)

let copie t =
let n = Array.length t in
let u = Array.make n t.(0) in
  for k = 0 to n-1 do
    u.(k) <- t.(k)
  done ;
u ;;

(*Addition en place du vecteur b au vecteur a*)
let add a b =
let n = Array.length a in
for k = 0 to n-1 do
a.(k) <- a.(k) +. b.(k)
done;;

(*Produit scalaire canonique de a et b*)

let scal a b =
let bigN = Array.length b in
let s = ref 0. in
for k = 0 to bigN-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(*bigNormalisation du vecteur x*)

let normalisation x =
let bigN = Array.length x in
let s = ref 0. in
for k = 0 to bigN-1 do
s:= !s +. (x.(k))*.(x.(k))
done;
s := sqrt (!s) ;
for k = 0 to bigN-1 do
x.(k) <- x.(k)/.(!s)
done ;;

(*bigTeste si u est le vecteur nul*)

let test_nul u =
let n = Array.length u - 1 in
let bool = ref true in
for k = 0 to n do
bool := (!bool)&&(u.(k) = 0.)
done ;
!bool;;

(*bigT.(k) désigne la kieme colonne de bigT. Recherche du plus petit produit scalaire*)

let miniscal bigT w =
let p = Array.length bigT - 1 in
let k = ref 0 in
for i = 0 to p do
if scal bigT.(i) w < scal bigT.(!k) w then k:=i done;
!k ;;

(*Marge d'un vecteur*)

let marge bigT w =
let j = miniscal bigT w and bigN = sqrt(scal w w) in (scal bigT.(j) w)/.bigN ;;

(*Cet algorithme renvoie le résultat de n itérations du perceptron avec la marge obtenue, en gardant en mémoire la meilleure solution rencontrée jusque là. Si le vecteur solution passe par 0, c'est inclassifiable d'après Gordan.*)

let perceptron_avance bigT n =
  let rec aux y a bigT x n = match n with
  |0 -> (y, a)
  |n ->
      let j = miniscal bigT x in
      add x bigT.(j) ;
      if test_nul x
      then raise (Failure "inclassifiable")
      else
      begin
      if a < marge bigT x
      then let z = copie x in normalisation z ; aux z (marge bigT z) bigT x (n-1) else aux y a bigT x (n-1)
      end ;
in let bigN = Array.length bigT.(0) in
let w = Array.make bigN 0. in
w.(0) <- 1. ; aux w (marge bigT w) bigT w n ;;



