(*addition du vecteur b au vecteur a*)
let add a b =
let n = Array.length a in
for k = 0 to n-1 do
a.(k) <- a.(k) + b.(k)
done;;

(*produit scalaire canonique de a et b*)

let scal a b =
let n = Array.length b in
let s = ref 0 in
for k = 0 to n-1 do
s:= !s + (a.(k))*(b.(k))
done;
!s;;

(*teste si tous les coefficients d'une matrice colonne valent 0*)

let all_zero m =
let n = Array.length m and bool = ref true in
for k=0 to n-1 do
if m.(k) <> 0 then bool := false
done ;
!bool;;

(*T.(k) désigne la kieme colonne de T
Cet algorithme s'il converge renvoie un vecteur poids convenable*)

let perceptron t =
let n = Array.length t.(0) and p = Array.length t in
let w = Array.make n 0 in
let bool = ref false in
while not (!bool) do
let d = Array.make n 0 in
for k = 0 to p-1 do
if (scal t.(k) w) <= 0 then (add d t.(k))
done;
if all_zero d then bool := true else add w d
done;
w;;

