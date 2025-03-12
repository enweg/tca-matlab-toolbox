function results = testTCA(dynarePath)
    addpath("./tests/");
    resultsUtils = run(utilsTest);
    resultsSimplifying = run(simplifyingTest);
    resultsSystemsForm = run(systemsFormTest);
    resultsTransmission = run(transmissionTest);
    results = [resultsUtils, resultsSimplifying, resultsSystemsForm, resultsTransmission];
    if nargin == 0
        warning("testTCA: Not testing Dynare functions because Dynare path was not provided.");
        return;
    end

    addpath("./dynare/")
    addpath(dynarePath)
    resultsDynare = run(dynareTest);
    results = [results, resultsDynare];
end
