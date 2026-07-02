function plan = buildfile
%BUILDFILE Build automation for gofCopula for MATLAB.

import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.CleanTask
import matlab.buildtool.tasks.TestTask

plan = buildplan(localfunctions);

% Built-in tasks provide incremental execution and standard CI artifacts.
plan("clean") = CleanTask;
plan("check") = CodeIssuesTask(["src" "examples" "docs" "tools"], ...
    WarningThreshold=0, ...
    Results="results/code-issues.sarif");
plan("test") = TestTask("tests", ...
    SourceFiles="src", ...
    Dependencies="check", ...
    TestResults="results/test-results.xml") ...
    .addCodeCoverage(["results/coverage.xml" "results/coverage.mat"]);

% TestTask records coverage but does not summarize or assess it. The
% custom task ENFORCES the statement-coverage threshold.
plan("coverage").Dependencies = "test";
plan("coverage").Inputs = "results/coverage.mat";

plan.DefaultTasks = ["check" "test" "coverage"];
end

function coverageTask(context)
%COVERAGETASK Enforce the statement-coverage threshold.

coverageFile = fullfile(context.Plan.RootFolder,"results","coverage.mat");
if ~isfile(coverageFile)
    context.log("Coverage data not found; run the test task first.");
    return
end

data = load(coverageFile);
[summary,~] = coverageSummary(data.coverage,"statement");
executed = sum(summary(:,1));
total = sum(summary(:,2));
lineRate = executed / total;
context.log(sprintf("Statement coverage: %.1f%% (%d/%d)", ...
    100*lineRate,executed,total));

threshold = 0.80;
context.assertTrue(lineRate >= threshold, sprintf( ...
    "Statement coverage %.1f%% is below the required %.0f%%.", ...
    100*lineRate,100*threshold));
end
