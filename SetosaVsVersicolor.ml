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






(*traitement des donnees*)

let f v =
let n = Array.length v in
let u = Array.make n 0. in
  for k= 0 to n-1 do
    u.(k) = v.(k)
  done ;
let s = ref 0. in
  for k = 0 to n-1 do
    s := !s +. u.(k) *. u.(k)
  done ;
s := sqrt(!s) ;
for k = 0 to n-1 do
  u.(k) = u.(k) /. !s
done ;
u ;;

let normalisation_donnees bigA =
let p = Array.length bigA in
let n = Array.length bigA.(0) in
let bigB = Array.make_matrix p n 0. in
  for k = 0 to p-1 do
    bigB.(k) = f bigA.(k)
  done ;
bigB;;

(*Donnees pretraitees*)

let g v = match v.(4) with
| 1. -> ()
| -1. -> for k = 0 to 3 do
           v.(k) = -.v.(k)
         done ;;

let treat bigT =
let n = Array.length bigT in
  for k = 0 to n-1 do
    g bigT.(k)
  done ;;

let bigT =
[| [|5.1;3.5;1.4;0.2;1.|];
[|4.9;3.0;1.4;0.2;1.|];
[|4.7;3.2;1.3;0.2;1.|];
[|4.6;3.1;1.5;0.2;1.|];
[|5.0;3.6;1.4;0.2;1.|];
[|5.4;3.9;1.7;0.4;1.|];
[|4.6;3.4;1.4;0.3;1.|];
[|5.0;3.4;1.5;0.2;1.|];
[|4.4;2.9;1.4;0.2;1.|];
[|4.9;3.1;1.5;0.1;1.|];
[|5.4;3.7;1.5;0.2;1.|];
[|4.8;3.4;1.6;0.2;1.|];
[|4.8;3.0;1.4;0.1;1.|];
[|4.3;3.0;1.1;0.1;1.|];
[|5.8;4.0;1.2;0.2;1.|];
[|5.7;4.4;1.5;0.4;1.|];
[|5.4;3.9;1.3;0.4;1.|];
[|5.1;3.5;1.4;0.3;1.|];
[|5.7;3.8;1.7;0.3;1.|];
[|5.1;3.8;1.5;0.3;1.|];
[|5.4;3.4;1.7;0.2;1.|];
[|5.1;3.7;1.5;0.4;1.|];
[|4.6;3.6;1.0;0.2;1.|];
[|5.1;3.3;1.7;0.5;1.|];
[|4.8;3.4;1.9;0.2;1.|];
[|5.0;3.0;1.6;0.2;1.|];
[|5.0;3.4;1.6;0.4;1.|];
[|5.2;3.5;1.5;0.2;1.|];
[|5.2;3.4;1.4;0.2;1.|];
[|4.7;3.2;1.6;0.2;1.|];
[|4.8;3.1;1.6;0.2;1.|];
[|5.4;3.4;1.5;0.4;1.|];
[|5.2;4.1;1.5;0.1;1.|];
[|5.5;4.2;1.4;0.2;1.|];
[|4.9;3.1;1.5;0.1;1.|];
[|5.0;3.2;1.2;0.2;1.|];
[|5.5;3.5;1.3;0.2;1.|];
[|4.9;3.1;1.5;0.1;1.|];
[|4.4;3.0;1.3;0.2;1.|];
[|5.1;3.4;1.5;0.2;1.|];
[|5.0;3.5;1.3;0.3;1.|];
[|4.5;2.3;1.3;0.3;1.|];
[|4.4;3.2;1.3;0.2;1.|];
[|5.0;3.5;1.6;0.6;1.|];
[|5.1;3.8;1.9;0.4;1.|];
[|4.8;3.0;1.4;0.3;1.|];
[|5.1;3.8;1.6;0.2;1.|];
[|4.6;3.2;1.4;0.2;1.|];
[|5.3;3.7;1.5;0.2;1.|];
[|5.0;3.3;1.4;0.2;1.|];
[|7.0;3.2;4.7;1.4;-1.|];
[|6.4;3.2;4.5;1.5;-1.|];
[|6.9;3.1;4.9;1.5;-1.|];
[|5.5;2.3;4.0;1.3;-1.|];
[|6.5;2.8;4.6;1.5;-1.|];
[|5.7;2.8;4.5;1.3;-1.|];
[|6.3;3.3;4.7;1.6;-1.|];
[|4.9;2.4;3.3;1.0;-1.|];
[|6.6;2.9;4.6;1.3;-1.|];
[|5.2;2.7;3.9;1.4;-1.|];
[|5.0;2.0;3.5;1.0;-1.|];
[|5.9;3.0;4.2;1.5;-1.|];
[|6.0;2.2;4.0;1.0;-1.|];
[|6.1;2.9;4.7;1.4;-1.|];
[|5.6;2.9;3.6;1.3;-1.|];
[|6.7;3.1;4.4;1.4;-1.|];
[|5.6;3.0;4.5;1.5;-1.|];
[|5.8;2.7;4.1;1.0;-1.|];
[|6.2;2.2;4.5;1.5;-1.|];
[|5.6;2.5;3.9;1.1;-1.|];
[|5.9;3.2;4.8;1.8;-1.|];
[|6.1;2.8;4.0;1.3;-1.|];
[|6.3;2.5;4.9;1.5;-1.|];
[|6.1;2.8;4.7;1.2;-1.|];
[|6.4;2.9;4.3;1.3;-1.|];
[|6.6;3.0;4.4;1.4;-1.|];
[|6.8;2.8;4.8;1.4;-1.|];
[|6.7;3.0;5.0;1.7;-1.|];
[|6.0;2.9;4.5;1.5;-1.|];
[|5.7;2.6;3.5;1.0;-1.|];
[|5.5;2.4;3.8;1.1;-1.|];
[|5.5;2.4;3.7;1.0;-1.|];
[|5.8;2.7;3.9;1.2;-1.|];
[|6.0;2.7;5.1;1.6;-1.|];
[|5.4;3.0;4.5;1.5;-1.|];
[|6.0;3.4;4.5;1.6;-1.|];
[|6.7;3.1;4.7;1.5;-1.|];
[|6.3;2.3;4.4;1.3;-1.|];
[|5.6;3.0;4.1;1.3;-1.|];
[|5.5;2.5;4.0;1.3;-1.|];
[|5.5;2.6;4.4;1.2;-1.|];
[|6.1;3.0;4.6;1.4;-1.|];
[|5.8;2.6;4.0;1.2;-1.|];
[|5.0;2.3;3.3;1.0;-1.|];
[|5.6;2.7;4.2;1.3;-1.|];
[|5.7;3.0;4.2;1.2;-1.|];
[|5.7;2.9;4.2;1.3;-1.|];
[|6.2;2.9;4.3;1.3;-1.|];
[|5.1;2.5;3.0;1.1;-1.|];
[|5.7;2.8;4.1;1.3;-1.|] |] ;;

treat bigT ;;
let bigA = normalisation_donnees bigT ;;

let w,a = perceptron_avance bigA 10 ;;
print_string "perceptron 10\n";;
print_float a;;
print_string "\n";;
print_string "[| ";;
for ii = 0 to Array.length w - 1 do
print_float w.(ii);
print_string " ; ";
done;
print_string "|]\n";;

