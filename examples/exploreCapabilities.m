%% Explore supported copula/test combinations
% Query the capability table before selecting tests for a dataset.

%% Complete capability table

capabilities = gofcopula.CopulaTestTable();
disp(capabilities)

%% Tests for a five-dimensional Gumbel copula

tests = gofcopula.gofTest4Copula("gumbel",5);
disp(tests)

%% Families supported by the White test

families = gofcopula.gofCopula4Test("gofWhite");
disp(families)

