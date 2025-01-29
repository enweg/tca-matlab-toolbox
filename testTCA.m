function results = testTCA()
    addpath("./tests/");
    resultsUtils = run(utilsTest);
    resultsSimplifying = run(simplifyingTest);
    resultsSystemsForm = run(systemsFormTest);
    results = [resultsUtils, resultsSimplifying, resultsSystemsForm];
end
