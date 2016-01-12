(*Algorithme Pal Rujan
Les matrices sont des float vect vect, les vecteurs des float vect. M.(i) est la ieme colonne de M*)

(*Présence d'un élément dans une liste*)

let rec présence x l = match l with
| [] -> false
| t::q -> (x=t)||(présence x q) ;;

(*Copie conforme d'un tableau*)

let copie t =
let n = vect_length t in
let u = make_vect n t.(0) in
  for k = 0 to n-1 do
    u.(k) <- t.(k)
  done ;
u ;;

(*Préliminaires d'algèbre linéaire*)

(*Produit scalaire de a et b*)

let scal a b =
let N = vect_length b in
let s = ref 0. in
for k = 0 to N-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(*Multiplication d'un vecteur par une matrice*)

let multiplication A X =
let n = vect_length A.(0) in
let p = vect_length A in
let w = make_vect n 0. in
  for a = 0 to n-1 do
    for b = 0 to p-1 do
      w.(a) <- w.(a) +. (X.(b)) *. (A.(b).(a))
    done ;
  done ;
 w ;;

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

(*Teste si m est le vecteur nul*)

let test_nul m =
let n = vect_length m and bool = ref true in
for k=0 to n-1 do
bool := (!bool)&&(m.(k) = 0.)
done ;
!bool;;

(*Cherche un produit scalaire inférieur à alpha_m. S'il y en a un, renvoie le numéro du vecteur correspondant (et non pas de la colonne dans T), sinon renvoie 0.*)

let contrainte_violée x alpha T =
let p = vect_length T in
let k = ref 0 in
  for i = 1 to p do
    if scal x T.(i-1) < alpha
    then match !k with
    |0 -> k := i
    |n -> if scal x T.(i-1) < scal x T.(n-1) then k := i ;
  done ;
!k ;;

(*Pivot de Gauss, résolution du système AX = B*)

(*Opération élémentaire du pivot de Gauss : soustraction de b fois la ligne i à la ligne j*)

let opération_élémentaire A X i j b =
let p = vect_length A in
  for k = 0 to p-1 do
    A.(k).(j) <- A.(k).(j) -. b*.A.(k).(i)
  done ;
X.(j) <- X.(j) -. b*.X.(i) ;
;;

(*Echange de deux lignes du système*)

let echange_elementaire A i j k =
let mem = A.(k).(i) in
A.(k).(i) <- A.(k).(j) ; A.(k).(j) <- mem ;;

let echange_lignes A X i j =
let mem = X.(i) in
X.(i) <- X.(j) ; X.(j) <- mem ;
let p = vect_length A in
  for k = 0 to p-1 do
   echange_elementaire A i j k
  done ;
;;

(*Résolution AX = B dans le cas où A est diagonale supérieure. Si pas de solution, false et un vecteur quelconque. Sinon, true et le vecteur solution.*)

let rec résolution_substitution A B = match vect_length A with
|1 ->
if A.(0).(0) = 0. then
  begin
    if B.(0) <> 0. then (false, B) else (true, [|0.|])
  end
else (true, [| B.(0) /. A.(0).(0)|])
|n when (A.(n-1).(n-1) = 0.)&&(B.(n-1) <> 0.) -> (false, B)
|n when A.(n-1).(n-1) = 0. ->
let B_1 = make_vect (n-1) 0. in
  for k = 0 to (n-2) do
    B_1.(k) <- B.(k)
  done ;
let A_1 = make_matrix (n-1) (n-1) 0. in
  for a = 0 to n-2 do
    for b = 0 to n-2 do
    A_1.(a).(b) <- A.(a).(b)
    done ;
  done ;
let (bool, X_1) = résolution_substitution A_1 B_1 in 
let X = make_vect n 0. in
  for k = 0 to n-2 do
    X.(k) <- X_1.(k)
  done ;
(bool, X)
|n ->
let y = B.(n-1) /. A.(n-1).(n-1) in
let B_1 = make_vect (n-1) 0. in
  for k = 0 to (n-2) do
    B_1.(k) <- B.(k) -. y *. A.(n-1).(k)
  done ;
let A_1 = make_matrix (n-1) (n-1) 0. in
  for a = 0 to n-2 do
    for b = 0 to n-2 do
      A_1.(a).(b) <- A.(a).(b)
    done ;
  done ;
let bool, X_1 = résolution_substitution A_1 B_1 in 
let X = make_vect n y in
  for k = 0 to n-2 do
    X.(k) <- X_1.(k)
  done ;
(bool, X)
;;

let résolution_pivot_gauss A B =
let rec aux A B n = match n with
|n when n = vect_length A -> résolution_substitution A B
|n ->
let p = vect_length A in
let pivot = ref p in
for i = p-1 downto n do
 if A.(n).(i) <> 0. then pivot := i
done ;
match !pivot with
|x when x=p -> aux A B (n+1)
|_ ->
  echange_lignes A B !pivot n ;
    for j = n+1 to p-1 do
      opération_élémentaire A B n j (A.(n).(j)/.A.(n).(n))
    done ;
aux A B (n+1)
in aux A B 0
;;

(*Algorithme de Pal Rujan.
On code A_m par un tableau de longueur p+1. La case 0 contient le cardinal de A_m et la case i contient 1 si c_i appartient à A_m et 0 sinon.*)

(*Correction de A.(0). Cet algorithme est en principe inutile mais se révèle nécessaire en pratique.*)

let comptage A =
let n = vect_length A in
let j = ref 0 in
  for k = 1 to n-1 do
    if A.(k) = 1 then j := !j + 1
  done ;
A.(0) <- !j ;;

(*Matrice contenant les c_i où i appartient à A_m, et matrice de l'équation à résoudre pour trouver les poids.*)

let matrices_intermédiaires A T Q =
let n = vect_length T in
let k = A.(0) in
let p = vect_length Q in
let Z = make_matrix (k+1) (k+1) 0. in
  for i = 1 to k do
    Z.(i).(0) <- 1. ; Z.(0).(i) <- -1.
  done ;
let x = ref 1 in
  for a = 1 to p do
    if A.(a) = 1 then
    begin
    let y = ref 1 in
      for b = 1 to p do
        if A.(b) = 1 then
        begin
        Z.(!x).(!y) <- Q.(a-1).(b-1) ; y := !y + 1
        end ;
      done ;
    x := !x + 1 ;
   end ;
  done ;
let S = make_matrix k n 0. in
let x = ref 0 in
  for a = 1 to p do
    if A.(a) = 1 then
    begin
      S.(!x) <- T.(a-1) ; x := !x + 1 ;
    end ;
  done ;
S,Z;;

(*Mise à jour de l'active set : A_m <- A_m+1*)

let active_set A v j =
let p = vect_length A in
let x = ref 0 in
 for i = 0 to vect_length v - 1 do
   if (v.(i) < v.(!x))&&(i <> j) then x := i
 done ;
if v.(!x) > 0.
then 0
else
let m = ref (!x + 1) in let k = ref 1 in
  while !m > 0 do
    if A.(!k) = 1 then m := !m - 1 ;
    k := !k + 1 ;
  done ;
A.(!k-1) <- 0 ; comptage A ;
!k-1 ;;


(*Retirer de A_m une contrainte différente de k*)

let retirer_contrainte A k =
let bool = ref true in
let x = ref 1 in
  while !bool do
    if (A.(!x) = 0)||(!x = k) then x := !x + 1 else bool := false
  done ;
A.(!x) <- 0 ; comptage A ;;

(*Renvoie les coefficients lambda_i >= 0 tels que sum lambda_i c_i = 0*)

let coefficients_convexes A u =
let p = vect_length A in
let v = make_vect (p-1) 0. in
let x = ref 0 in
  for k = 1 to p-1 do
    match A.(k) with
    | 1 -> v.(k-1) <- u.(!x) ; x := !x + 1
    | 0 -> ()
  done ;
v ;;

(*Marge d'un vecteur sur l'active set*)

let marge_active_set A u T =
let x = ref 1 in
  while A.(!x) = 0 do
    if A.(!x) = 0 then x := !x + 1
  done ;
let j = ref !x in
  for k = 1 to vect_length A - 1 do
    if (A.(k) = 1)&&(scal u T.(k-1) <= scal u T.(!j-1))
    then j := k
  done ;
scal u T.(!j-1) ;;

(*Renvoie si la solution de l'équation principale existe, le vecteur obtenu en supprimant le premier coefficient de la solution, ce premier coefficient, et S_m*)

let solution_intermédiaire T Q A =
let S,Z = matrices_intermédiaires A T Q in
let v = make_vect (A.(0)+1) 0. in
v.(0) <- 1. ;
let (bool, lambda) = résolution_pivot_gauss Z v in
let u = make_vect A.(0) 0. in
for k = 0 to A.(0)-1 do
  u.(k) <- lambda.(k+1)
done ;
bool, u, lambda.(0), S
;;

(*Met à jour l'active set jusqu'à ce que tous les lambda_i soient positifs strictement*)

let rec update_active_set T Q A j =
match solution_intermédiaire T Q A with
| bool, _, z, _ when (not bool)||(z < 0.) -> retirer_contrainte A j ; comptage A ; update_active_set T Q A j
| true, u, _, _ ->
match active_set A u j with
| 0 -> comptage A
| k -> comptage A ; update_active_set T Q A j ;;

(*Une étape de l'algorithme. Renvoie l'active set.*)

let rec aux_rujan T Q A l compteur =
if présence A l
then false, A, compteur else
let B = copie A in
let (bool, u, z, S) = solution_intermédiaire T Q A in
match z with
| 0. -> true, A, (compteur + 1)
| _ -> let w = multiplication S u in
normalisation w ;
let c = marge_active_set A w T in
match contrainte_violée w c T with
|0 -> true, A, (compteur + 1)
|k -> A.(k) <- 1 ; comptage A ; update_active_set T Q A k ; comptage A ; aux_rujan T Q A (B::l) (compteur+1)
;;

let algorithme_pal_rujan T =
  let rec aux T k acc = 
  let p = vect_length T in
  let Q = make_matrix p p 0. in
    for a = 0 to p-1 do
      for b = 0 to p-1 do
        Q.(a).(b) <- scal T.(a) T.(b)
      done ;
    done ;
  let A = make_vect (p+1) 0 in
  A.(0) <- 1 ; A.(k) <- 1 ;
  let (bool, B, n) = aux_rujan T Q A [] 0 in
  match bool, solution_intermédiaire T Q B with
  |false, (_, _, _, _) -> aux T (k+1) (acc+n)
  |_, (_, u, 0., _) -> (coefficients_convexes B u), 0., (acc+n)
  |_, (_, u, _, S) -> let w = multiplication S u in normalisation w ;
match marge_active_set B w T with
|c when c < 0. -> aux T (k+1) (acc+n)
|c -> (w, c, acc+n)
in aux T 1 0;;











(*Données prétraitées*)

let T = [|[|5.1;3.5;1.4;0.2;1.|];
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
[|6.3;3.3;6.0;2.5;-1.|];
[|5.8;2.7;5.1;1.9;-1.|];
[|7.1;3.0;5.9;2.1;-1.|];
[|6.3;2.9;5.6;1.8;-1.|];
[|6.5;3.0;5.8;2.2;-1.|];
[|7.6;3.0;6.6;2.1;-1.|];
[|4.9;2.5;4.5;1.7;-1.|];
[|7.3;2.9;6.3;1.8;-1.|];
[|6.7;2.5;5.8;1.8;-1.|];
[|7.2;3.6;6.1;2.5;-1.|];
[|6.5;3.2;5.1;2.0;-1.|];
[|6.4;2.7;5.3;1.9;-1.|];
[|6.8;3.0;5.5;2.1;-1.|];
[|5.7;2.5;5.0;2.0;-1.|];
[|5.8;2.8;5.1;2.4;-1.|];
[|6.4;3.2;5.3;2.3;-1.|];
[|6.5;3.0;5.5;1.8;-1.|];
[|7.7;3.8;6.7;2.2;-1.|];
[|7.7;2.6;6.9;2.3;-1.|];
[|6.0;2.2;5.0;1.5;-1.|];
[|6.9;3.2;5.7;2.3;-1.|];
[|5.6;2.8;4.9;2.0;-1.|];
[|7.7;2.8;6.7;2.0;-1.|];
[|6.3;2.7;4.9;1.8;-1.|];
[|6.7;3.3;5.7;2.1;-1.|];
[|7.2;3.2;6.0;1.8;-1.|];
[|6.2;2.8;4.8;1.8;-1.|];
[|6.1;3.0;4.9;1.8;-1.|];
[|6.4;2.8;5.6;2.1;-1.|];
[|7.2;3.0;5.8;1.6;-1.|];
[|7.4;2.8;6.1;1.9;-1.|];
[|7.9;3.8;6.4;2.0;-1.|];
[|6.4;2.8;5.6;2.2;-1.|];
[|6.3;2.8;5.1;1.5;-1.|];
[|6.1;2.6;5.6;1.4;-1.|];
[|7.7;3.0;6.1;2.3;-1.|];
[|6.3;3.4;5.6;2.4;-1.|];
[|6.4;3.1;5.5;1.8;-1.|];
[|6.0;3.0;4.8;1.8;-1.|];
[|6.9;3.1;5.4;2.1;-1.|];
[|6.7;3.1;5.6;2.4;-1.|];
[|6.9;3.1;5.1;2.3;-1.|];
[|5.8;2.7;5.1;1.9;-1.|];
[|6.8;3.2;5.9;2.3;-1.|];
[|6.7;3.3;5.7;2.5;-1.|];
[|6.7;3.0;5.2;2.3;-1.|];
[|6.3;2.5;5.0;1.9;-1.|];
[|6.5;3.0;5.2;2.0;-1.|];
[|6.2;3.4;5.4;2.3;-1.|];
[|5.9;3.0;5.1;1.8;-1.|]|] ;;

let f v =
let n = vect_length v in
let u = make_vect n 0. in
  for k= 0 to n-1 do
    u.(k) <- v.(k)
  done ;
let s = ref 0. in
  for k = 0 to n-1 do
    s := !s +. u.(k) *. u.(k)
  done ;
s := sqrt(!s) ;
for k = 0 to n-1 do
  u.(k) <- u.(k) /. !s
done ;
u ;;

let normalisation_données A =
let p = vect_length A in
let n = vect_length A.(0) in
let B = make_matrix p n 0. in
  for k = 0 to p-1 do
    B.(k) <- f A.(k)
  done ;
B;;

let g v = match v.(4) with
| 1. -> ()
| -1. -> for k = 0 to 3 do
           v.(k) <- -.v.(k)
         done ;;

let treat T =
let n = vect_length T in
  for k = 0 to n-1 do
    g T.(k)
  done ;;

treat T ;;

let A = normalisation_données T ;;

algorithme_pal_rujan A ;;


