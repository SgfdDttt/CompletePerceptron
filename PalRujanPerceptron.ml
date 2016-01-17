open LinearAlgebra;;

(*Pal Rujan algorithm
Matrices are float vect vect, vectors are float vect. M.(i) is the i-th column of M*)

(*Presence of an element in a list*)
let rec presence x l = match l with
| [] -> false
| t::q -> (x=t)||(presence x q) ;;

(*Within the columns of tt, looks for a vector whose scalar product with x is the smallest.*)
(*If that scalar product is smaller than alpha (the margin), it returns the number of that vector (and not the index of the column in tt), otherwise it returns 0*)
let smaller_scalar_product x alpha tt =
let p = Array.length tt in
let k = ref 0 in
  for i = 1 to p do
    if scal x tt.(i-1) < alpha
    then match !k with
    |0 -> k := i
    |n -> if scal x tt.(i-1) < scal x tt.(n-1) then k := i ;
  done ;
!k ;;

(*Pal Rujan algorithm.
The active set A_m is an array of length p+1. The slot of index 0 contains the cardinal of A_m. The slot i contains 1 if c_i belongs to A_m, and 0 otherwise.*)

(*Computing de A_m.(0). This algorithm is useful in practice.*)
let updating_cardinal aa =
let n = Array.length aa in
let j = ref 0 in
  for k = 1 to n-1 do
    if aa.(k) = 1 then j := !j + 1
  done ;
aa.(0) <- !j ;;

(*Matrice contenant les c_i où i appartient à bigA_m, et matrice de l'equation à resoudre pour trouver les poids.*)

let matrices_intermediaires aa tt qq =
let n = Array.length tt in
let k = aa.(0) in
let p = Array.length qq in
let zz = Array.make_matrix (k+1) (k+1) 0. in
  for i = 1 to k do
    zz.(i).(0) <- 1. ; zz.(0).(i) <- -1.
  done ;
let x = ref 1 in
  for a = 1 to p do
    if aa.(a) = 1 then
    begin
    let y = ref 1 in
      for b = 1 to p do
        if aa.(b) = 1 then
        begin
        zz.(!x).(!y) <- qq.(a-1).(b-1) ; y := !y + 1
        end ;
      done ;
    x := !x + 1 ;
   end ;
  done ;
let ss = Array.make_matrix k n 0. in
let x = ref 0 in
  for a = 1 to p do
    if aa.(a) = 1 then
    begin
      ss.(!x) <- tt.(a-1) ; x := !x + 1 ;
    end ;
  done ;
ss,zz;;

(*Mise à jour de l'active set : A_m = A_m+1*)
let active_set aa v j =
let p = Array.length aa in
let x = ref 0 in
 for i = 0 to Array.length v - 1 do
   if (v.(i) < v.(!x))&&(i <> j) then x := i
 done ;
if v.(!x) > 0.
then 0
else
let m = ref (!x + 1) in let k = ref 1 in
  while !m > 0 do
    if aa.(!k) = 1 then m := !m - 1 ;
    k := !k + 1 ;
  done ;
aa.(!k-1) <- 0 ; updating_cardinal aa ;
!k-1 ;;


(*Retirer de A_m une contrainte differente de k*)

let retirer_contrainte aa k =
let bool = ref true in
let x = ref 1 in
  while !bool do
    if (aa.(!x) = 0)||(!x = k) then x := !x + 1 else bool := false
  done ;
aa.(!x) <- 0 ; updating_cardinal aa ;;

(*Renvoie les coefficients lambda_i >= 0 tels que sum lambda_i c_i = 0*)
let coefficients_convexes aa u =
let p = Array.length aa in
let v = Array.make (p-1) 0. in
let x = ref 0 in
  for k = 1 to p-1 do
    match aa.(k) with
    | 1 -> v.(k-1) <- u.(!x) ; x := !x + 1
    | 0 -> ()
		| _ -> raise (Failure "coefficients_convexes: unexpected value")
  done ;
v ;;

(*Marge d'un vecteur sur l'active set*)
let marge_active_set aa u tt =
let x = ref 1 in
  while aa.(!x) = 0 do
    if aa.(!x) = 0 then x := !x + 1
  done ;
let j = ref !x in
  for k = 1 to Array.length aa - 1 do
    if (aa.(k) = 1)&&(scal u tt.(k-1) <= scal u tt.(!j-1))
    then j := k
  done ;
scal u tt.(!j-1) ;;

(*Renvoie si la solution de l'equation principale existe, le vecteur obtenu en supprimant le premier coefficient de la solution, ce premier coefficient, et bigS_m*)
let solution_intermediaire tt qq aa =
let ss,zz = matrices_intermediaires aa tt qq in
let v = Array.make (aa.(0)+1) 0. in
v.(0) <- 1. ;
let (bool, lambda) = solve_gaussian_pivot zz v in
let u = Array.make aa.(0) 0. in
for k = 0 to aa.(0)-1 do
  u.(k) <- lambda.(k+1)
done ;
bool, u, lambda.(0), ss
;;

(*Met à jour l'active set jusqu'à ce que tous les lambda_i soient positifs strictement*)
let rec update_active_set tt qq aa j =
match solution_intermediaire tt qq aa with
| bool, _, z, _ when (not bool)||(z < 0.) -> retirer_contrainte aa j ; updating_cardinal aa ; update_active_set tt qq aa j
| true, u, _, _ ->
match active_set aa u j with
| 0 -> updating_cardinal aa
| k -> updating_cardinal aa ; update_active_set tt qq aa j
;;

(*Une etape de l'algorithme. Renvoie l'active set.*)
let rec aux_rujan tt qq a_liste l compteur =
if presence a_liste l
then false, a_liste, compteur else
let bb = copy a_liste in
let (bool, u, z, ss) = solution_intermediaire tt qq a_liste in
match z with
| 0. -> true, a_liste, (compteur + 1)
| _ -> let w = matrix_mult ss u in
normalisation w ;
let c = marge_active_set a_liste w tt in
match smaller_scalar_product w c tt with
|0 -> true, a_liste, (compteur + 1)
|k -> a_liste.(k) <- 1 ; updating_cardinal a_liste ; update_active_set tt qq a_liste k ; updating_cardinal a_liste ; aux_rujan tt qq a_liste (bb::l) (compteur+1)
;;

let algorithme_pal_rujan tt =
  let rec aux tt k acc = 
  let p = Array.length tt in
  let qq = Array.make_matrix p p 0. in
    for a = 0 to p-1 do
      for b = 0 to p-1 do
        qq.(a).(b) <- scal tt.(a) tt.(b)
      done ;
    done ;
  let aa = Array.make (p+1) 0 in
  aa.(0) <- 1 ; aa.(k) <- 1 ;
  let (bool, bb, n) = aux_rujan tt qq aa [] 0 in
  match bool, solution_intermediaire tt qq bb with
  |false, (_, _, _, _) -> aux tt (k+1) (acc+n)
  |_, (_, u, 0., _) -> (coefficients_convexes bb u), 0., (acc+n)
  |_, (_, u, _, ss) -> let w = matrix_mult ss u in normalisation w ;
match marge_active_set bb w tt with
|c when c < 0. -> aux tt (k+1) (acc+n)
|c -> (w, c, acc+n)
in aux tt 1 0;;




