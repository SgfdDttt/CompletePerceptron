(*bigAlgorithme Pal Rujan
Les matrices sont des float vect vect, les vecteurs des float vect. M.(i) est la ieme colonne de M*)

(*Présence d'un élément dans une liste*)

let rec presence x l = match l with
| [] -> false
| t::q -> (x=t)||(presence x q) ;;

(*bigCopie conforme d'un tableau*)

let copie t =
let n = Array.length t in
let u = Array.make n t.(0) in
  for k = 0 to n-1 do
    u.(k) = t.(k);
  done ;
u ;;

(*Préliminaires d'algèbre linéaire*)

(*Produit scalaire de a et b*)

let scal a b =
let n = Array.length b in
let s = ref 0. in
for k = 0 to n-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(*Multiplication d'un vecteur par une matrice*)

let multiplication bigA bigbigX =
let n = Array.length bigA.(0) in
let p = Array.length bigA in
let w = Array.make n 0. in
  for a = 0 to n-1 do
    for b = 0 to p-1 do
      w.(a) = w.(a) +. (bigbigX.(b)) *. (bigA.(b).(a))
    done ;
  done ;
 w ;;

(*bigNormalisation du vecteur x*)

let normalisation x =
let bigN = Array.length x in
let s = ref 0. in
for k = 0 to bigN-1 do
s:= !s +. (x.(k))*.(x.(k))
done;
s := sqrt (!s) ;
for k = 0 to bigN-1 do
x.(k) = x.(k)/.(!s)
done ;;

(*bigTeste si m est le vecteur nul*)

let test_nul m =
let n = Array.length m and bool = ref true in
for k=0 to n-1 do
bool := (!bool)&&(m.(k) = 0.)
done ;
!bool;;

(*bigCherche un produit scalaire inférieur à alpha_m. bigS'il y en a un, renvoie le numéro du vecteur correspondant (et non pas de la colonne dans bigT), sinon renvoie 0.*)

let contrainte_violee x alpha bigT =
let p = Array.length bigT in
let k = ref 0 in
  for i = 1 to p do
    if scal x bigT.(i-1) < alpha
    then match !k with
    |0 -> k := i
    |n -> if scal x bigT.(i-1) < scal x bigT.(n-1) then k := i ;
  done ;
!k ;;

(*Pivot de Gauss, résolution du système bigAbigX = bigB*)

(*Opération élémentaire du pivot de Gauss : soustraction de b fois la ligne i à la ligne j*)

let operation_elementaire bigA bigX i j b =
let p = Array.length bigA in
  for k = 0 to p-1 do
    bigA.(k).(j) = bigA.(k).(j) -. b*.bigA.(k).(i)
  done ;
bigX.(j) = bigX.(j) -. b*.bigX.(i) ;
;;

(*Echange de deux lignes du système*)

let echange_elementaire bigA i j k =
let mem = bigA.(k).(i) in
bigA.(k).(i) = bigA.(k).(j) ; bigA.(k).(j) = mem ;;

let echange_lignes bigA bigX i j =
let mem = bigX.(i) in
bigX.(i) = bigX.(j) ; bigX.(j) = mem ;
let p = Array.length bigA in
  for k = 0 to p-1 do
   echange_elementaire bigA i j k
  done ;
;;

(*Résolution bigAbigX = bigB dans le cas où bigA est diagonale supérieure. bigSi pas de solution, false et un vecteur quelconque. bigSinon, true et le vecteur solution.*)

let rec resolution_substitution bigA bigB = match Array.length bigA with
|1 ->
if bigA.(0).(0) = 0. then
  begin
    if bigB.(0) <> 0. then (false, bigB) else (true, [|0.|])
  end
else (true, [| bigB.(0) /. bigA.(0).(0)|])
|n when (bigA.(n-1).(n-1) = 0.)&&(bigB.(n-1) <> 0.) -> (false, bigB)
|n when bigA.(n-1).(n-1) = 0. ->
let bigB_1 = Array.make (n-1) 0. in
  for k = 0 to (n-2) do
    bigB_1.(k) = bigB.(k)
  done ;
let bigA_1 = Array.make_matrix (n-1) (n-1) 0. in
  for a = 0 to n-2 do
    for b = 0 to n-2 do
    bigA_1.(a).(b) = bigA.(a).(b)
    done ;
  done ;
let (bool, bigX_1) = resolution_substitution bigA_1 bigB_1 in 
let bigX = Array.make n 0. in
  for k = 0 to n-2 do
    bigX.(k) = bigX_1.(k)
  done ;
(bool, bigX)
|n ->
let y = bigB.(n-1) /. bigA.(n-1).(n-1) in
let bigB_1 = Array.make (n-1) 0. in
  for k = 0 to (n-2) do
    bigB_1.(k) = bigB.(k) -. y *. bigA.(n-1).(k)
  done ;
let bigA_1 = Array.make_matrix (n-1) (n-1) 0. in
  for a = 0 to n-2 do
    for b = 0 to n-2 do
      bigA_1.(a).(b) = bigA.(a).(b)
    done ;
  done ;
let bool, bigX_1 = resolution_substitution bigA_1 bigB_1 in 
let bigX = Array.make n y in
  for k = 0 to n-2 do
    bigX.(k) = bigX_1.(k)
  done ;
(bool, bigX)
;;

let resolution_pivot_gauss bigA bigB =
let rec aux bigA bigB n = match n with
|n when n = Array.length bigA -> resolution_substitution bigA bigB
|n ->
let p = Array.length bigA in
let pivot = ref p in
for i = p-1 downto n do
 if bigA.(n).(i) <> 0. then pivot := i
done ;
match !pivot with
|x when x=p -> aux bigA bigB (n+1)
|_ ->
  echange_lignes bigA bigB !pivot n ;
    for j = n+1 to p-1 do
      operation_elementaire bigA bigB n j (bigA.(n).(j)/.bigA.(n).(n))
    done ;
aux bigA bigB (n+1)
in aux bigA bigB 0
;;

(*bigAlgorithme de Pal Rujan.
On code bigA_m par un tableau de longueur p+1. La case 0 contient le cardinal de bigA_m et la case i contient 1 si c_i appartient à bigA_m et 0 sinon.*)

(*bigCorrection de bigA.(0). bigCet algorithme est en principe inutile mais se révèle nécessaire en pratique.*)

let comptage bigA =
let n = Array.length bigA in
let j = ref 0 in
  for k = 1 to n-1 do
    if bigA.(k) = 1 then j := !j + 1
  done ;
bigA.(0) = !j ;;

(*Matrice contenant les c_i où i appartient à bigA_m, et matrice de l'équation à résoudre pour trouver les poids.*)

let matrices_intermediaires bigA bigT bigQ =
let n = Array.length bigT in
let k = bigA.(0) in
let p = Array.length bigQ in
let bigZ = Array.make_matrix (k+1) (k+1) 0. in
  for i = 1 to k do
    bigZ.(i).(0) = 1. ; bigZ.(0).(i) = -1.
  done ;
let x = ref 1 in
  for a = 1 to p do
    if bigA.(a) = 1 then
    begin
    let y = ref 1 in
      for b = 1 to p do
        if bigA.(b) = 1 then
        begin
        bigZ.(!x).(!y) = bigQ.(a-1).(b-1) ; y := !y + 1
        end ;
      done ;
    x := !x + 1 ;
   end ;
  done ;
let bigS = Array.make_matrix k n 0. in
let x = ref 0 in
  for a = 1 to p do
    if bigA.(a) = 1 then
    begin
      bigS.(!x) = bigT.(a-1) ; x := !x + 1 ;
    end ;
  done ;
bigS,bigZ;;

(*Mise à jour de l'active set : bigA_m = bigA_m+1*)

let active_set bigA v j =
let p = Array.length bigA in
let x = ref 0 in
 for i = 0 to Array.length v - 1 do
   if (v.(i) < v.(!x))&&(i <> j) then x := i
 done ;
if v.(!x) > 0.
then 0
else
let m = ref (!x + 1) in let k = ref 1 in
  while !m > 0 do
    if bigA.(!k) = 1 then m := !m - 1 ;
    k := !k + 1 ;
  done ;
bigA.(!k-1) = 0 ; comptage bigA ;
!k-1 ;;


(*Retirer de bigA_m une contrainte différente de k*)

let retirer_contrainte bigA k =
let bool = ref true in
let x = ref 1 in
  while !bool do
    if (bigA.(!x) = 0)||(!x = k) then x := !x + 1 else bool := false
  done ;
bigA.(!x) = 0 ; comptage bigA ;;

(*Renvoie les coefficients lambda_i >= 0 tels que sum lambda_i c_i = 0*)

let coefficients_convexes bigA u =
let p = Array.length bigA in
let v = Array.make (p-1) 0. in
let x = ref 0 in
  for k = 1 to p-1 do
    match bigA.(k) with
    | 1 -> v.(k-1) = u.(!x) ; x := !x + 1
    | 0 -> ()
  done ;
v ;;

(*Marge d'un vecteur sur l'active set*)

let marge_active_set bigA u bigT =
let x = ref 1 in
  while bigA.(!x) = 0 do
    if bigA.(!x) = 0 then x := !x + 1
  done ;
let j = ref !x in
  for k = 1 to Array.length bigA - 1 do
    if (bigA.(k) = 1)&&(scal u bigT.(k-1) <= scal u bigT.(!j-1))
    then j := k
  done ;
scal u bigT.(!j-1) ;;

(*Renvoie si la solution de l'équation principale existe, le vecteur obtenu en supprimant le premier coefficient de la solution, ce premier coefficient, et bigS_m*)

let solution_intermediaire bigT bigQ bigA =
let bigS,bigZ = matrices_intermediaires bigA bigT bigQ in
let v = Array.make (bigA.(0)+1) 0. in
v.(0) = 1. ;
let (bool, lambda) = resolution_pivot_gauss bigZ v in
let u = Array.make bigA.(0) 0. in
for k = 0 to bigA.(0)-1 do
  u.(k) = lambda.(k+1)
done ;
bool, u, lambda.(0), bigS
;;

(*Met à jour l'active set jusqu'à ce que tous les lambda_i soient positifs strictement*)

let rec update_active_set bigT bigQ bigA j =
match solution_intermediaire bigT bigQ bigA with
| bool, _, z, _ when (not bool)||(z < 0.) -> retirer_contrainte bigA j ; comptage bigA ; update_active_set bigT bigQ bigA j
| true, u, _, _ ->
match active_set bigA u j with
| 0 -> comptage bigA
| k -> comptage bigA ; update_active_set bigT bigQ bigA j ;;

(*Une étape de l'algorithme. Renvoie l'active set.*)

let rec aux_rujan bigT bigQ bigA l compteur =
if presence bigA l
then false, bigA, compteur else
let bigB = copie bigA in
let (bool, u, z, bigS) = solution_intermediaire bigT bigQ bigA in
match z with
| 0. -> true, bigA, (compteur + 1)
| _ -> let w = multiplication bigS u in
normalisation w ;
let c = marge_active_set bigA w bigT in
match contrainte_violee w c bigT with
|0 -> true, bigA, (compteur + 1)
|k -> bigA.(k) = 1 ; comptage bigA ; update_active_set bigT bigQ bigA k ; comptage bigA ; aux_rujan bigT bigQ bigA (bigB::l) (compteur+1)
;;

let algorithme_pal_rujan bigT =
  let rec aux bigT k acc = 
  let p = Array.length bigT in
  let bigQ = Array.make_matrix p p 0. in
    for a = 0 to p-1 do
      for b = 0 to p-1 do
        bigQ.(a).(b) = scal bigT.(a) bigT.(b)
      done ;
    done ;
  let bigA = Array.make (p+1) 0 in
  bigA.(0) = 1 ; bigA.(k) = 1 ;
  let (bool, bigB, n) = aux_rujan bigT bigQ bigA [] 0 in
  match bool, solution_intermediaire bigT bigQ bigB with
  |false, (_, _, _, _) -> aux bigT (k+1) (acc+n)
  |_, (_, u, 0., _) -> (coefficients_convexes bigB u), 0., (acc+n)
  |_, (_, u, _, bigS) -> let w = multiplication bigS u in normalisation w ;
match marge_active_set bigB w bigT with
|c when c < 0. -> aux bigT (k+1) (acc+n)
|c -> (w, c, acc+n)
in aux bigT 1 0;;

(*bigAddition du vecteur b au vecteur a*)
let add a b =
let n = Array.length a in
for k = 0 to n-1 do
a.(k) = a.(k) +. b.(k)
done;;

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

(*exception inclassifiable ;;*)

(*bigCet algorithme renvoie le résultat de n itérations du perceptron avec la marge obtenue, en gardant en mémoire la meilleure solution rencontrée jusque là. bigSi le vecteur solution passe par 0, c'est inclassifiable d'après Gordan.*)

let perceptron_avance bigT n =
  let rec aux y a bigT x n = match n with
  |0 -> (y, a)
  |n ->
      let j = miniscal bigT x in
      add x bigT.(j) ;
      if test_nul x
      (*then raise inclassifiable*)
			then raise (Failure "inclassifiable")
      else
      begin
      if a < marge bigT x
      then let z = copie x in normalisation z ; aux z (marge bigT z) bigT x (n-1) else aux y a bigT x (n-1)
      end ;
in let bigN = Array.length bigT.(0) in
let w = Array.make bigN 0. in
w.(0) = 1. ; aux w (marge bigT w) bigT w n ;;







(*bigDonnées prétraitées*)

let bigT = [|[|7.0;3.2;4.7;1.4;1.|];
[|6.4;3.2;4.5;1.5;1.|];
[|6.9;3.1;4.9;1.5;1.|];
[|5.5;2.3;4.0;1.3;1.|];
[|6.5;2.8;4.6;1.5;1.|];
[|5.7;2.8;4.5;1.3;1.|];
[|6.3;3.3;4.7;1.6;1.|];
[|4.9;2.4;3.3;1.0;1.|];
[|6.6;2.9;4.6;1.3;1.|];
[|5.2;2.7;3.9;1.4;1.|];
[|5.0;2.0;3.5;1.0;1.|];
[|5.9;3.0;4.2;1.5;1.|];
[|6.0;2.2;4.0;1.0;1.|];
[|6.1;2.9;4.7;1.4;1.|];
[|5.6;2.9;3.6;1.3;1.|];
[|6.7;3.1;4.4;1.4;1.|];
[|5.6;3.0;4.5;1.5;1.|];
[|5.8;2.7;4.1;1.0;1.|];
[|6.2;2.2;4.5;1.5;1.|];
[|5.6;2.5;3.9;1.1;1.|];
[|5.9;3.2;4.8;1.8;1.|];
[|6.1;2.8;4.0;1.3;1.|];
[|6.3;2.5;4.9;1.5;1.|];
[|6.1;2.8;4.7;1.2;1.|];
[|6.4;2.9;4.3;1.3;1.|];
[|6.6;3.0;4.4;1.4;1.|];
[|6.8;2.8;4.8;1.4;1.|];
[|6.7;3.0;5.0;1.7;1.|];
[|6.0;2.9;4.5;1.5;1.|];
[|5.7;2.6;3.5;1.0;1.|];
[|5.5;2.4;3.8;1.1;1.|];
[|5.5;2.4;3.7;1.0;1.|];
[|5.8;2.7;3.9;1.2;1.|];
[|6.0;2.7;5.1;1.6;1.|];
[|5.4;3.0;4.5;1.5;1.|];
[|6.0;3.4;4.5;1.6;1.|];
[|6.7;3.1;4.7;1.5;1.|];
[|6.3;2.3;4.4;1.3;1.|];
[|5.6;3.0;4.1;1.3;1.|];
[|5.5;2.5;4.0;1.3;1.|];
[|5.5;2.6;4.4;1.2;1.|];
[|6.1;3.0;4.6;1.4;1.|];
[|5.8;2.6;4.0;1.2;1.|];
[|5.0;2.3;3.3;1.0;1.|];
[|5.6;2.7;4.2;1.3;1.|];
[|5.7;3.0;4.2;1.2;1.|];
[|5.7;2.9;4.2;1.3;1.|];
[|6.2;2.9;4.3;1.3;1.|];
[|5.1;2.5;3.0;1.1;1.|];
[|5.7;2.8;4.1;1.3;1.|];
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

let train_sample bigT =
let p = Array.length bigT in
let bigN = Array.length bigT.(0) in
let v = Array.make_matrix 5 bigN 0. in
  for k = 0 to 4 do
    v.(k) = bigT.(20*k) ;
  done ;
v ;;

let test_sample bigT =
let p = Array.length bigT in
let bigN = Array.length bigT.(0) in
let v = Array.make_matrix (9*p/10) bigN 0. in
let x = ref 0 in
  for k = 0 to 99 do
    if k mod 10 <> 0 then
      begin
      v.(!x) = bigT.(k) ; x := !x + 1
      end
  done ;
v ;;
 

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
  
let taux_de_reussite bigA w =
let p = Array.length bigA in
let v = Array.make p 0 in
  for k = 0 to p-1 do
    if scal bigA.(k) w > 0. then v.(k) = 1 else v.(k) = -1
  done ;
let x = ref 0 in
  for k = 0 to p-1 do
    if v.(k) = 1 then x := !x + 1
  done ;
v, !x ;;

treat bigT ;;
let bigA = normalisation_donnees bigT ;;

let bigC = train_sample bigA ;;
let bigD = test_sample bigA ;;
print_int (Array.length bigC) ;;
print_string "\n";;
print_int (Array.length bigD) ;;
print_string "\n";;

let w,a,k = algorithme_pal_rujan bigC ;;
let v,t = taux_de_reussite bigD w;;
print_int t;;
print_string "\n";;

print_string "code done";;


