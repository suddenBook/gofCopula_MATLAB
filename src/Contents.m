% gofCopula for MATLAB
% Version 0.5.0 (R2025b) 02-Jul-2026
%
% Models and results
%   gofcopula.CopulaModel       - Copula model specification.
%   gofcopula.GofResult         - Goodness-of-fit result value object.
%
% Coordinating interfaces
%   gofcopula.gof               - Run tests over copula families.
%   gofcopula.gofco             - Run tests for a CopulaModel.
%   gofcopula.gofCustomTest     - Bootstrap a user-defined statistic.
%   gofcopula.gofCheckTime      - Estimate bootstrap execution time.
%   gofcopula.gofGetHybrid      - Combine tests by Bonferroni correction.
%   gofcopula.gofOutputHybrid   - Format hybrid and individual results.
%
% Empirical and Kendall-process tests
%   gofcopula.gofCvM            - Empirical-copula Cramer-von Mises test.
%   gofcopula.gofKS             - Empirical-copula Kolmogorov-Smirnov test.
%   gofcopula.gofKendallCvM     - Kendall-process Cramer-von Mises test.
%   gofcopula.gofKendallKS      - Kendall-process Kolmogorov-Smirnov test.
%
% Rosenblatt-transform tests
%   gofcopula.gofRosenblattSnB  - Breymann aggregation test.
%   gofcopula.gofRosenblattSnC  - Cramer-von Mises aggregation test.
%   gofcopula.gofRosenblattGamma - Gamma aggregation test.
%   gofcopula.gofRosenblattChisq - Chi-square aggregation test.
%
% Other test families
%   gofcopula.gofKernel         - Kernel-based test.
%   gofcopula.gofWhite          - Information-matrix test.
%   gofcopula.gofPIOSTn         - PIOS jackknife statistic.
%   gofcopula.gofPIOSRn         - PIOS approximation statistic.
%   gofcopula.gofArchmSnB       - Archimedean Breymann test.
%   gofcopula.gofArchmSnC       - Archimedean Cramer-von Mises test.
%   gofcopula.gofArchmGamma     - Archimedean gamma test.
%   gofcopula.gofArchmChisq     - Archimedean chi-square test.
%
% Capability queries
%   gofcopula.CopulaTestTable   - Supported family/test dimensions.
%   gofcopula.gofCopula4Test    - Families supported by a test.
%   gofcopula.gofTest4Copula    - Tests supported by a family/dimension.
%
% Copula primitives
%   gofcopula.copulaCDF         - Copula distribution function.
%   gofcopula.copulaPDF         - Copula (log-)density.
%   gofcopula.copulaRandom      - Copula random number generation.
%   gofcopula.rosenblatt        - Sequential Rosenblatt transform.
%   gofcopula.runTest           - Run one bootstrap test (advanced).
%   gofcopula.runTestSerial     - Serial-dependence-robust bootstrap test.
%
% Documentation
%   ../docs/GettingStarted.m    - Executable introduction.
%   ../docs/MigrationGuide.md   - R-to-MATLAB API mapping.
%   ../docs/CapabilityMatrix.md - Family, test, and dimension support.
%   ../docs/Numerics.md         - Precision and NumericMode semantics.
%   ../docs/SerialDependence.md - Serial-dependence-robust bootstrap.
