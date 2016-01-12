(*addition du vecteur b au vecteur a*)
let add a b =
let n = vect_length a in
for k = 0 to n-1 do
a.(k) <- a.(k) + b.(k)
done;;

(*produit scalaire canonique de a et b*)

let scal a b =
let N = vect_length b in
let s = ref 0 in
for k = 0 to N-1 do
s:= !s + (a.(k))*(b.(k))
done;
!s;;

(*teste si tous les coefficients d'une matrice colonne valent 0*)

let all_zero m =
let n = vect_length m and bool = ref true in
for k=0 to n-1 do
if m.(k) <> 0 then bool := false
done ;
!bool;;

(*T.(k) désigne la kieme colonne de T
Cet algorithme s'il converge renvoie un vecteur poids convenable*)

let perceptron T =
let N = vect_length T.(0) and p = vect_length T in
let W = make_vect N 0 in
let bool = ref false in
while not (!bool) do
let D = make_vect N 0 in
for k = 0 to p-1 do
if (scal T.(k) W) <= 0 then (add D T.(k))
done;
if all_zero D then bool := true else add W D
done;
W;;

