(* Verification of the Archimedean generator-derivative representations in
   +gofcopula/+internal/+copulas/archimedean.m at high precision.
   Run: wolframscript -file generatorDerivatives.wls *)

Off[N::meprec]; (* residuals are exact zeros; N cannot tell them from 0 *)
chopped[expr_] := Abs[Chop[N[expr, 40], 10^-35]];

(* Frank: psi(x) = -(1/t) Log[1 + (E^-t - 1) E^-x];
   claim: psi^(k)(x) = (-1)^k Li_{-(k-1)}(qq)/t with qq = (1-E^-t) E^-x. *)
psiF[x_] := -(1/t) Log[1 + (Exp[-t] - 1) Exp[-x]];
frankErr[k_] := Max[Table[
  chopped[(D[psiF[x], {x, k}]
      - (-1)^k PolyLog[-(k - 1), (1 - Exp[-t]) Exp[-x]]/t)
      /. {t -> tv, x -> xv}],
  {tv, {1/2, 3, 8}}, {xv, {1/10, 1, 4}}]];
Do[Print["Frank psi^(", k, ") polylog identity max |error|: ",
  frankErr[k]], {k, 1, 4}];

(* AMH: psi(x) = (1-t)/(E^x - t);
   claim: psi^(k)(x) = ((1-t)/t) (-1)^k Li_{-k}(t E^-x). *)
psiA[x_] := (1 - t)/(Exp[x] - t);
amhErr[k_] := Max[Table[
  chopped[(D[psiA[x], {x, k}]
      - ((1 - t)/t) (-1)^k PolyLog[-k, t Exp[-x]])
      /. {t -> tv, x -> xv}],
  {tv, {1/4, 7/10}}, {xv, {1/10, 1, 4}}]];
Do[Print["AMH psi^(", k, ") polylog identity max |error|: ",
  amhErr[k]], {k, 1, 4}];

(* Gumbel: psi(x) = Exp[-x^(1/t)]. The implementation builds the polynomial
   P with psi^(k)(x) = psi(x) x^-k P(x^(1/t)) via the recurrence
   P_{n+1} = -n P_n - a z P_n + a z P_n', a = 1/t. Spot-check orders 2..4. *)
psiG[x_] := Exp[-x^(1/t)];
gumbelPoly[order_, a_] := Module[{p = {0, -a}, deriv, n},
  (* p holds coefficients of P in z, ascending; start P_1(z) = -a z *)
  Do[
    deriv = Table[i p[[i + 1]], {i, 1, Length[p] - 1}];
    p = -n PadRight[p, Length[p] + 1]
        - a PadRight[Prepend[p, 0], Length[p] + 1]
        + a PadRight[Prepend[deriv, 0], Length[p] + 1];
    , {n, 1, order - 1}];
  p];
gumbelErr[k_] := Max[Table[Module[{a = 1/tv, poly, claim},
  poly = gumbelPoly[k, a];
  claim = psiG[x] x^-k Sum[poly[[i + 1]] (x^(1/t))^i, {i, 0, Length[poly] - 1}]
      /. t -> tv;
  chopped[((D[psiG[x], {x, k}] /. t -> tv) - claim) /. x -> xv]],
  {tv, {3/2, 4}}, {xv, {1/2, 2, 10}}]];
Do[Print["Gumbel psi^(", k, ") recurrence max |error|: ",
  gumbelErr[k]], {k, 2, 4}];
