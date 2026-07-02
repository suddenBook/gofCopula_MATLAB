(* Symbolic verification of the analytic extreme-value Pickands derivatives
   implemented in +gofcopula/+internal/+copulas/evDependence.m.
   Run: wolframscript -file evDerivatives.wls  (output committed alongside). *)

(* Galambos: A = 1 - S^(-1/t), S = w^-t + (1-w)^-t *)
S = w^-t + (1 - w)^-t;
A = 1 - S^(-1/t);
A2claimed = (1 + t)*((w^(-t - 2) + (1 - w)^(-t - 2))*S^(-1/t - 1)
    - (w^(-t - 1) - (1 - w)^(-t - 1))^2*S^(-1/t - 2));
Print["Galambos A'' matches: ",
  FullSimplify[D[A, {w, 2}] - A2claimed == 0, Assumptions -> 0 < w < 1 && t > 0]];
A1claimed = -S^(-1/t - 1)*(w^(-t - 1) - (1 - w)^(-t - 1));
Print["Galambos A' matches: ",
  FullSimplify[D[A, w] - A1claimed == 0, Assumptions -> 0 < w < 1 && t > 0]];

(* Husler-Reiss: A = w Phi(z+) + (1-w) Phi(z-),
   z+- = 1/t +- (t/2) Log[w/(1-w)] *)
zp = 1/t + (t/2) Log[w/(1 - w)];
zm = 1/t - (t/2) Log[w/(1 - w)];
phi[x_] := Exp[-x^2/2]/Sqrt[2 Pi];
Print["HR identity w*phi(z+) = (1-w)*phi(z-): ",
  FullSimplify[w phi[zp] - (1 - w) phi[zm] == 0, Assumptions -> 0 < w < 1 && t > 0]];
Ahr = w CDF[NormalDistribution[], zp] + (1 - w) CDF[NormalDistribution[], zm];
Print["HR A' = Phi(z+) - Phi(z-): ",
  FullSimplify[D[Ahr, w] - (CDF[NormalDistribution[], zp] - CDF[NormalDistribution[], zm]) == 0,
    Assumptions -> 0 < w < 1 && t > 0]];
(* A'' claim, verified at high precision on a grid (FullSimplify is slow here). *)
Off[N::meprec]; (* residuals are exact zeros; N cannot tell them from 0 *)
A2hr = t/(2 w (1 - w)) (phi[zp] + phi[zm]);
errors = Table[
  Abs[Chop[N[(D[Ahr, {w, 2}] - A2hr) /. {t -> tv, w -> wv}, 40], 10^-35]],
  {tv, {3/10, 3/2, 5}}, {wv, {1/100, 1/4, 1/2, 3/4, 99/100}}];
Print["HR A'' max |error| over grid (40-digit): ", Max[errors]];

(* t-EV: A = w G(x) + (1-w) G(y), G = t CDF with nu+1 df; the implemented
   product-rule derivative is compared numerically at high precision. *)
r = (w/(1 - w))^(1/nu);
s = Sqrt[(nu + 1)/(1 - th^2)];
xx = (r - th) s; yy = (1/r - th) s;
Atev = w CDF[StudentTDistribution[nu + 1], xx] + (1 - w) CDF[StudentTDistribution[nu + 1], yy];
g[u_] := PDF[StudentTDistribution[nu + 1], u];
k = 1/(nu w (1 - w));
A1tev = CDF[StudentTDistribution[nu + 1], xx] - CDF[StudentTDistribution[nu + 1], yy]
    + s k (w r g[xx] - (1 - w) g[yy]/r);
errs = Table[
  Abs[Chop[N[(D[Atev, w] - A1tev) /. {nu -> 5, th -> thv, w -> wv}, 40], 10^-35]],
  {thv, {-3/10, 1/2}}, {wv, {1/20, 1/2, 19/20}}];
Print["tEV A' max |error| over grid (40-digit): ", Max[errs]];
