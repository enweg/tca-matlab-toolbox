function results = testTCA()
    addpath("./tests/");
    resultsUtils = run(utilsTest);
    resultsSimplifying = run(simplifyingTest);
    resultsSystemsForm = run(systemsFormTest);
    resultsTransmission = run(transmissionTest);
    results = [resultsUtils, resultsSimplifying, resultsSystemsForm, resultsTransmission];
end
