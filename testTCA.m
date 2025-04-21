function results = testTCA(dynarePath)
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
