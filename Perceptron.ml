(*add vector b to vector a in place*)
let adda a b =
let n = Array.length a in
for k = 0 to n-1 do
a.(k) <- a.(k) +. b.(k)
done;;

(*scalar product of a and b*)
let scal a b =
let n = Array.length b in
let s = ref 0. in
for k = 0 to n-1 do
s:= !s +. (a.(k))*.(b.(k))
done;
!s;;

(* checks to see if all the coefficients in a vector are 0*)
let all_zero m =
let n = Array.length m and bool = ref true in
for k=0 to n-1 do
if m.(k) <> 0. then bool := false
done ;
!bool;;

(* If it converges, this algorithm returns a possible weight vector*)
let perceptron t =
let n = Array.length t.(0) and p = Array.length t in
let w = Array.make n 0. in
let finished = ref false in
while not (!finished) do
let delta = Array.make n 0. in
for k = 0 to p-1 do
if (scal t.(k) w) <= 0. then (adda delta t.(k))
done;
if all_zero delta then finished := true else adda w delta
done;
w;;

