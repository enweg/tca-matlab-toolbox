function results = testTCA(dynarePath)
    % `testTCA` Run the test suite of the toolbox.
    %
    %   `results = testTCA()` runs the tests of the TCA functions and of the
    %   `VAR`, `SVAR` and `LP` models. `results = testTCA(dynarePath)` also
    %   runs the `DSGE` tests using the Dynare installation in `dynarePath`.
    %
    %   ## Arguments
    %   - `dynarePath` (string, optional): Path to the Dynare `matlab` folder.
    %     If not provided, the `DSGE` tests are skipped.
    %
    %   ## Returns
    %   - `results` (matlab.unittest.TestResult): Results of all tests that ran.
    %
    %   ## Notes
    %   - Must be run from the root folder of the toolbox.

    addpath("./tests/");
    resultsUtils = run(utilsTest);
    resultsSimplifying = run(simplifyingTest);
    resultsSystemsForm = run(systemsFormTest);
    resultsTransmission = run(transmissionTest);

    addpath("./models/");
    resultsVAR = run(modelVARTest);
    resultsSVAR = run(modelSVARTest);
    resultsLP = run(modelLPTest);

    results = [resultsUtils, resultsSimplifying, resultsSystemsForm, resultsTransmission, resultsVAR, resultsSVAR, resultsLP];

    if nargin == 0
        warning("testTCA: Not testing Dynare functions because Dynare path was not provided.");
        return;
    end
    addpath(dynarePath)
    resultsDSGE = run(modelDSGETest);
    results = [results, resultsDSGE];
end
