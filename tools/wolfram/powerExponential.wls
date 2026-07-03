(* Symbolic/numeric validation of the power-exponential (PE) copula math used in
   +gofcopula/+internal/+elliptical/. Reference: Gomez, Gomez-Villegas & Marin
   (1998). Run: wolframscript -file powerExponential.wls (output committed). *)

kd = n Gamma[n/2]/(Pi^(n/2) Gamma[1 + n/(2 beta)] 2^(1 + n/(2 beta)));

(* Density normalization and the beta=1 -> Gaussian reduction. *)
radial = Integrate[t^(n/2 - 1) Exp[-(1/2) t^beta], {t, 0, Infinity},
   Assumptions -> n > 0 && beta > 0];
Print["k normalizes the spherical PE density: ",
  FullSimplify[FunctionExpand[kd (2 Pi^(n/2)/Gamma[n/2]) (1/2) radial == 1],
    Assumptions -> n > 0 && beta > 0]];
Print["beta=1 gives k=(2Pi)^(-n/2) (Gaussian): ",
  FullSimplify[FunctionExpand[(kd /. beta -> 1) == (2 Pi)^(-n/2)], Assumptions -> n > 0]];

(* Sampler fact: T = R^(2beta) ~ Gamma(shape n/2beta, scale 2). *)
hR = n/(Gamma[1 + n/(2 beta)] 2^(n/(2 beta))) r^(n - 1) Exp[-(1/2) r^(2 beta)];
pdfT = Simplify[(hR /. r -> T^(1/(2 beta))) (1/(2 beta)) T^(1/(2 beta) - 1),
   Assumptions -> T > 0 && n > 0 && beta > 0];
Print["law of T=R^(2b) is Gamma(n/2b, scale 2): ",
  FullSimplify[FunctionExpand[pdfT == PDF[GammaDistribution[n/(2 beta), 2], T]],
    Assumptions -> T > 0 && n > 0 && beta > 0]];

(* Radial moment E[S^beta] = n/beta (S the squared Mahalanobis radius). *)
Print["E[(x'Sigma^-1 x)^beta] = n/beta: ", Simplify[(n/(2 beta)) 2 == n/beta]];

(* Copula entropy shape term: kappa(d,1)=0, i.e. Gaussian copula entropy = 1/2 log det R. *)
Aterm = FullSimplify[(1/2) (d/beta) - Log[kd /. n -> d]];
Print["A(d,1) = d/2 log(2 Pi E): ",
  FullSimplify[PowerExpand[(Aterm /. beta -> 1) - (d/2) Log[2 Pi E]] == 0, d > 0]];
Print["kappa(d,1) = A(d,1) - d*(1/2)log(2 Pi E) = 0: ",
  FullSimplify[PowerExpand[(Aterm /. beta -> 1) - d ((1/2) Log[2 Pi E])] == 0, d > 0]];

(* Non-closure: the 1-D margin of a 2-D PE is NOT univariate PE for beta != 1.
   If it were, -Log[margin] would be affine in |x|^(2 beta); show it is not. *)
b0 = 1/2; k2 = kd /. {n -> 2, beta -> b0};
marg[x_?NumericQ] := k2 NIntegrate[Exp[-(1/2) (x^2 + y^2)^b0], {y, -Infinity, Infinity}];
xs = {2/10, 6/10, 12/10, 2, 3};
nl = (-Log[marg[#]]) & /@ xs;
slope = (nl[[2]] - nl[[1]])/(xs[[2]] - xs[[1]]);
icpt = nl[[1]] - slope xs[[1]];
Print["margin of PE_2 is not PE_1 (nonzero residuals vs linear-in-|x|): ",
  N[Table[nl[[i]] - (icpt + slope xs[[i]]), {i, 3, 5}], 5]];
