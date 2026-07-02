(* Verification of the closed-form conditional inverses used by the fast
   samplers in +gofcopula/+internal/+copulas/randomCore.m.
   h(v|u) is the conditional CDF dC/du; v(w) must satisfy h(v(w)|u) == w. *)

(* Clayton *)
hClayton = u^(-t - 1) (u^-t + v^-t - 1)^(-1/t - 1);
vClayton = (1 + u^-t (q^(-t/(1 + t)) - 1))^(-1/t);
(* FullSimplify cannot collapse the nested-radical residual; verify at
   40-digit precision over a (u, q, t) grid instead. *)
Off[N::meprec]; (* residuals are exact zeros; N cannot tell them from 0 *)
resC = Max[Table[
  Abs[Chop[N[((hClayton /. v -> vClayton) - q) /. {u -> uv, q -> qv, t -> tv}, 40],
    10^-35]],
  {uv, {1/20, 1/2, 19/20}}, {qv, {1/20, 1/2, 19/20}}, {tv, {1/2, 5/2, 8}}]];
On[N::meprec];
Print["Clayton inverse max |residual| over grid (40-digit, chopped at 1e-35): ", resC];

(* Frank *)
hFrank = Exp[-t u] (Exp[-t v] - 1)/((Exp[-t] - 1) + (Exp[-t u] - 1) (Exp[-t v] - 1));
vFrank = -(1/t) Log[1 + q (Exp[-t] - 1)/((1 - q) Exp[-t u] + q)];
resF = FullSimplify[(hFrank /. v -> vFrank) - q,
  Assumptions -> 0 < u < 1 && 0 < q < 1 && t > 0];
Print["Frank inverse residual simplifies to: ", resF];

(* FGM: h(v|u) = v + t(2u-1)v(v-1) = Av^2 + (1-A)v with A = t(2u-1) *)
hFGM = v + t (2 u - 1) v (v - 1);
vFGM = ((a - 1) + Sqrt[(1 - a)^2 + 4 a q])/(2 a) /. a -> t (2 u - 1);
resG = FullSimplify[(hFGM /. v -> vFGM) - q,
  Assumptions -> 0 < u < 1 && 0 < q < 1 && -1 < t < 1 && t != 0 && u != 1/2];
Print["FGM inverse residual simplifies to: ", resG];
