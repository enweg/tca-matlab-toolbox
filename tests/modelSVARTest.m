function tests = modelVARTest
    tests = functiontests(localfunctions);
end

function testSVARBasics(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 1000; 
    trendExponents = [0];
    m = length(trendExponents);

    A0 = randn(k, k);
    B = 0.2 * randn(k, k * p + m);
    APlus = A0 * B;

    Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);
    shocks = randn(k, T);
    Y = SVAR.simulate(shocks, A0, APlus, 'trendExponents', trendExponents);

    model = SVAR(Y, p, 'trendExponents', trendExponents);

    assert(~model.isFitted())
    model.getDependent();
    model.getIndependent();
    model.getInputData();
    model.getVariableNames();
    assert(model.nobs() == T - p);
    assert(model.isStructural());

    verifyError(testCase, @() model.coeffs(), "SVAR:NotFitted");
    verifyError(testCase, @() model.coeffs(true), "SVAR:NotFitted");
    verifyError(testCase, @() model.fitted(), "SVAR:NotFitted");
    verifyError(testCase, @() model.residuals(), "SVAR:NotFitted");
    verifyError(testCase, @() model.makeCompanionMatrix(), "SVAR:NotFitted");
    verifyError(testCase, @() model.spectralRadius(), "SVAR:NotFitted");
    verifyError(testCase, @() model.isStable(), "SVAR:NotFitted");
    verifyError(testCase, @() model.aic(), "SVAR:NotFitted");
    verifyError(testCase, @() model.bic(), "SVAR:NotFitted");
    verifyError(testCase, @() model.hqc(), "SVAR:NotFitted");
    verifyError(testCase, @() model.sic(), "SVAR:NotFitted");

    method = Recursive();
    model.fit(method);
    assert(model.isFitted());
    % Below show no-longer throw any errors
    model.coeffs();
    model.coeffs(true);
    model.fitted();
    model.residuals();
    model.makeCompanionMatrix();
    model.spectralRadius();
    model.isStable();
    model.aic();
    model.bic();
    model.sic();
    model.hqc();
end

function testSVARCoefficientsRecursive(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 1000; 
    trendExponents = [0];
    m = length(trendExponents);

    A0 = randn(k, k);
    A0 = tril(A0);
    S = diag(sign(diag(A0)));
    A0 = A0 * S;
    Phi0 = inv(A0);
    B = 0.2 * randn(k, k * p + m);
    APlus = A0 * B;
    SigmaU = Phi0 * Phi0';

    [A0Test, APlusTest] = Recursive.identifyVAR_(B, SigmaU);
    assert(all(max(abs(A0 - A0Test), [], 'all') < sqrt(eps())));
    assert(all(max(abs(APlus - APlusTest), [], 'all') < sqrt(eps())));

    % Checking implementations

    Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);
    modelVAR = VAR(Y, p, 'trendExponents', trendExponents);
    modelVAR.fit();
    method = Recursive();
    [A0Test, APlusTest] = method.identify(modelVAR);

    modelSVAR = SVAR(Y, p, 'trendExponents', trendExponents);
    modelSVAR.fit(method);
    [A0Test, APlusTest] = coeffs(modelSVAR);
end

function testSVARInformationCriteria(testCase)

    rng(6150533);
    k = 3;
    p = 2;
    T = 10000; 
    trendExponents = [0];
    m = length(trendExponents);

    A0 = randn(k, k);
    A0 = tril(A0);
    S = diag(sign(diag(A0)));
    A0 = A0 * S;
    B = 0.2 * randn(k, k * p + m);
    APlus = A0 * B;

    Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);
    modelLarge = SVAR(Y, p + 10, 'trendExponents', trendExponents);

    method = Recursive();
    [modelBest, icTable] = modelLarge.fitAndSelect(method);
    assert(modelBest.p == p);
    [modelBest, icTable] = modelLarge.fitAndSelect(method, @VAR.aic_);
    assert(modelBest.p == p);
    [modelBest, icTable] = modelLarge.fitAndSelect(method, @VAR.bic_);
    assert(modelBest.p == p);
    [modelBest, icTable] = modelLarge.fitAndSelect(method, @VAR.sic_);
    assert(modelBest.p == p);
    [modelBest, icTable] = modelLarge.fitAndSelect(method, @VAR.hqc_);
    assert(modelBest.p == p);
end

function testSVARIRFRecursive(testCase)
      
    rng(6150533);
    k = 3;
    p = 2;
    T = 10000; 
    trendExponents = [0];
    m = length(trendExponents);

    A0 = randn(k, k);
    A0 = tril(A0);
    S = diag(sign(diag(A0)));
    A0 = A0 * S;
    Phi0 = inv(A0);
    B = 0.2 * randn(k, k * p + m);
    APlus = A0 * B;
    SigmaU = Phi0 * Phi0';

    irfs = Recursive.identifyVARIrfs_(B(:, (m+1):end), SigmaU, p, 10);
    irfsSVAR = SVAR.IRF_(A0, APlus(:, (m+1):end), p, 10);

    assert(all(max(abs(irfs - irfsSVAR), [], 'all') < sqrt(eps())));

    % IMPLEMENTATION TESTS
    Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);
    modelVAR = VAR(Y, p, 'trendExponents', trendExponents);
    modelVAR.fit()
    method = Recursive();
    irfs = modelVAR.IRF(10, 'identificationMethod', method);

    model = SVAR(Y, p, 'trendExponents', trendExponents);
    method = Recursive();
    model.fit(method);
    irfs = model.IRF(10);
end

function testSVARInternalInstrument(testCase)
    k = 3; 
    p = 2;
    T = 100000;
    trendExponents = [0];
    m = length(trendExponents);

    A0 = randn(k, k);
    % A0 = tril(A0);
    % S = diag(sign(diag(A0)));
    % A0 = A0 * S;
    Phi0 = inv(A0);
    B = 0.2 * randn(k, k * p + m);
    APlus = A0 * B;
    SigmaU = Phi0 * Phi0';

    % introducing instrument into SVAR
    % first variable is instrument 
    % first "shock" is measurement error
    Phi0Tilde = zeros(k + 1, k+1);
    Phi0Tilde(2:end, 2:end) = Phi0;
    Phi0Tilde(1, 1:2) = [0.1, 1];

    BTilde = zeros(k+1, (k+1)*p + m);
    BTilde(2:end, 1:m) = B(:, 1:m);
    tmp = B(:, (m+1):end);
    for pp = 1:p
        Bi = tmp(:, (k*(pp-1)+1):(pp*k));
        BTilde(2:end, ((k+1)*(pp-1)+2+m):(pp*(k+1)+m)) = Bi;
    end
    SigmaUTilde = Phi0Tilde * Phi0Tilde';

    irfsTrue = SVAR.IRF_(A0, APlus(:, (m+1):end), p, 10);

    method = InternalInstrument(2);
    irfs = method.identifyVARIrfs_(BTilde(:, (m+1):end), SigmaUTilde, p, 10);
    testDiff = irfs(2:end, 1, :) - irfsTrue(:, 1, :) ./ irfsTrue(1, 1, 1);
    assert(all(max(abs(testDiff), [], 'all') < sqrt(eps())));

    method = InternalInstrument(2, 'normalisingHorizon', 1);
    irfs = method.identifyVARIrfs_(BTilde(:, (m+1):end), SigmaUTilde, p, 10);
    testDiff = irfs(2:end, 1, :) - irfsTrue(:, 1, :) ./ irfsTrue(1, 1, 2);
    assert(all(max(abs(testDiff), [], 'all') < sqrt(eps())));

    % THE REMAINING TESTS ARE IMPLEMENTATION TESTS. 

    Y = VAR.simulate(T, BTilde, 'SigmaU', SigmaUTilde, 'trendExponents', trendExponents);
    modelVAR = VAR(Y, p, 'trendExponents', trendExponents);
    modelVAR.fit();

    method = InternalInstrument(2);
    irfObj = modelVAR.IRF(10, 'identificationMethod', method);
    testDiff = irfObj.irfs(2:end, 1, :) - irfsTrue(:, 1, :) ./ irfsTrue(1, 1, 1);
    max(abs(testDiff), [], 'all');

    method = InternalInstrument('Y2');
    irfObj = modelVAR.IRF(10, 'identificationMethod', method);
    testDiff = irfObj.irfs(2:end, 1, :) - irfsTrue(:, 1, :) ./ irfsTrue(1, 1, 1);
    max(abs(testDiff), [], 'all');
    
    method = InternalInstrument('Y2', 'normalisingHorizon', 1);
    irfObj = modelVAR.IRF(10, 'identificationMethod', method);
    testDiff = irfObj.irfs(2:end, 1, :) - irfsTrue(:, 1, :) ./ irfsTrue(1, 1, 2);
    max(abs(testDiff), [], 'all');

    method = InternalInstrument('Y3');
    irfObj = modelVAR.IRF(10, 'identificationMethod', method);
    testDiff = irfObj.irfs(2:end, 1, :) - irfsTrue(:, 1, :) ./ irfsTrue(2, 1, 1);
    max(abs(testDiff), [], 'all');

    method = InternalInstrument('Y3', 'instrument', 2);
    irfObj = modelVAR.IRF(10, 'identificationMethod', method);
    irfObj;
end
