function results = testTCA(dynarePath)
    addpath("./tests/");
    resultsUtils = run(utilsTest);
    resultsSimplifying = run(simplifyingTest);
    resultsSystemsForm = run(systemsFormTest);
    resultsTransmission = run(transmissionTest);
    resultsDynare = run(dynareTest);
    results = [resultsUtils, resultsSimplifying, resultsSystemsForm, resultsTransmission, resultsDynare];
end
