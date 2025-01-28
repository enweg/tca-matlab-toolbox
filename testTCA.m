function results = testTCA()
    addpath("./tests/");
    resultsUtils = run(utilsTest);
    resultsSimplifying = run(simplifyingTest);
    results = [resultsUtils, resultsSimplifying];
end
