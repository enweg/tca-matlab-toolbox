function tests = modelVARTest
    tests = functiontests(localfunctions);
end

function testBasicVARFunctions(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 10000;
    trendExponents = [0];
    B = 0.2 * randn(k, k*p + length(trendExponents));

    % Just testing if this version also runs
    errors = randn(k, T);
    Y = VAR.simulate(errors, B, 'trendExponents', trendExponents);
    % Now the more common version
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);

    assert(~model.isFitted())
    model.getDependent();
    model.getIndependent();
    model.getInputData();
    assert(model.nobs() == T - p);
    assert(~model.isStructural());

    verifyError(testCase, @() model.coeffs(), "VAR:NotFitted");
    verifyError(testCase, @() model.coeffs(true), "VAR:NotFitted");
    verifyError(testCase, @() model.fitted(), "VAR:NotFitted");
    verifyError(testCase, @() model.residuals(), "VAR:NotFitted");
    verifyError(testCase, @() model.ncoeffs(), "VAR:NotFitted");
    verifyError(testCase, @() model.makeCompanionMatrix(), "VAR:NotFitted");
    verifyError(testCase, @() model.spectralRadius(), "VAR:NotFitted");
    verifyError(testCase, @() model.isStable(), "VAR:NotFitted");
    verifyError(testCase, @() model.aic(), "VAR:NotFitted");
    verifyError(testCase, @() model.bic(), "VAR:NotFitted");
    verifyError(testCase, @() model.hqc(), "VAR:NotFitted");
    verifyError(testCase, @() model.sic(), "VAR:NotFitted");

    model.fit();
    assert(model.isFitted());
    % Below show no-longer throw any errors
    model.coeffs();
    model.coeffs(true);
    model.fitted();
    model.residuals();
    assert(model.ncoeffs() == k * (p*k + length(trendExponents)))
    model.makeCompanionMatrix();
    model.spectralRadius();
    model.isStable();
    model.aic();
    model.bic();
    model.sic();
    model.hqc();
end

function testVARCoefficientRecovery(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 100;
    trendExponents = [0];
    B = 0.2 * randn(k, k*p + length(trendExponents));

    errors = zeros(k, T);
    initial = 100 * ones(k*p, 1);
    Y = VAR.simulate(errors, B, 'trendExponents', trendExponents, 'initial', initial);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();

    Bhat = coeffs(model); 
    assert(all(max(abs(Bhat - B), [], 'all') < sqrt(eps())));
end

function testVARCovarianceRecovery(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 100000;
    trendExponents = 0:1;
    B = 0.2 * randn(k, k*p + length(trendExponents));
    SigmaU = [
        10 5 5
        5 10 5 
        5 5 10
    ];

    Y = VAR.simulate(T, B, 'SigmaU', SigmaU, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();

    testPctDiff = (model.SigmaU - SigmaU) ./ SigmaU;
    assert(all(max(abs(testPctDiff), [], 'all') < 1e-2));
end

function testVARInformationCriteria(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 10000;
    trendExponents = [0];
    B = 0.2 * randn(k, k*p + length(trendExponents));

    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);

    modelLarge = VAR(Y, p + 10, 'trendExponents', trendExponents);
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.aic_);
    assert(modelBest.p == p);
    icTable;
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.bic_);
    assert(modelBest.p == p);
    icTable;
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.sic_);
    assert(modelBest.p == p);
    icTable;
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.hqc_);
    assert(modelBest.p == p);
    icTable;


    
    SigmaU = [
        1 0.5 0.5;
        0.5 1 0.5; 
        0.5 0.5 1
    ];
    Y = VAR.simulate(T, B, 'SigmaU', SigmaU, 'trendExponents', trendExponents);

    modelLarge = VAR(Y, p + 10, 'trendExponents', trendExponents);
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.aic_);
    assert(modelBest.p == p);
    icTable;
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.bic_);
    assert(modelBest.p == p);
    icTable;
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.sic_);
    assert(modelBest.p == p);
    icTable;
    [modelBest, icTable] = modelLarge.fitAndSelect(@VAR.hqc_);
    assert(modelBest.p == p);
    icTable;
end

function testVARTrendExponents(testCase)
    % THESE ARE JUST IMPLEMENTATION TESTS TO AUTOMATICALLY CHECK IF ANY OF THE 
    % COMMONLY USED FUNCTIONS THROWS AN ERRORS. THE UNDERLYING FUNCTIONS 
    % ARE TESTED ELSEWHERE. 

    rng(6150533);
    k = 3;
    p = 2;
    T = 1000;

    % contstant and linear
    trendExponents = [0, 1];
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit()

    % polynomial trend
    trendExponents = 0:2;
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit()

    % no constant no trend
    trendExponents = [];
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit()

    % only trend
    trendExponents = [1];
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit()
end

function testVARLagLength(testCase)
    % THESE ARE JUST IMPLEMENTATION TESTS TO AUTOMATICALLY CHECK IF ANY OF THE 
    % COMMONLY USED FUNCTIONS THROWS AN ERRORS. THE UNDERLYING FUNCTIONS 
    % ARE TESTED ELSEWHERE. 

    rng(6150533);
    k = 3;
    trendExponents = 0:1;
    T = 10000;

    % no lags
    p = 0;
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();

    % one lag
    p = 1;
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();

    % many lags
    p = 10;
    B = 0.05 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();
end	

function testVARNumberVariables(testCase)
    % THESE ARE JUST IMPLEMENTATION TESTS TO AUTOMATICALLY CHECK IF ANY OF THE 
    % COMMONLY USED FUNCTIONS THROWS AN ERRORS. THE UNDERLYING FUNCTIONS 
    % ARE TESTED ELSEWHERE. 

    rng(6150533);
    p = 2;
    trendExponents = [0];
    T = 10000;

    % AR
    k = 1;
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();

    
    k = 2;
    B = 0.2 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();

    k = 20;
    B = 0.05 * randn(k, k*p + length(trendExponents));
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();
end

function testVARIRF(testCase)
    k = 3;
    p = 2;
    trendExponents = [0];
    m = length(trendExponents);
    B = 0.2 * randn(k, k*p + m);

    maxHorizon = 4;
    C = VAR.makeCompanionMatrix_(B, p, m);
    irfsCompanion = nan(k, k, maxHorizon + 1);
    for h = 0:maxHorizon
        tmp = C^h;
        irfsCompanion(:, :, h+1) = tmp(1:k, 1:k);
    end

    % Using exact coefficients
    irfs = VAR.IRF_(B(:, (m+1):end), p, maxHorizon);
    testDiff = irfs - irfsCompanion;
    assert(all(max(abs(testDiff), [], 'all') < sqrt(eps())));

    % THE FOLLOWING IS JUST AN IMPLEMENTATION TEST. VAR.IRF_ IS TESTED ABOVE.
    T = 10000;
    Y = VAR.simulate(T, B, 'trendExponents', trendExponents);
    model = VAR(Y, p, 'trendExponents', trendExponents);
    model.fit();
    irfsModel = model.IRF(maxHorizon);
end
