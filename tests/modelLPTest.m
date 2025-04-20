function tests = modelLPTest
    tests = functiontests(localfunctions);
end

function testLPConstruction(testCase)
    k = 4;
    T = 10;
    p = 2;

    Y = reshape(1:(T*k), T, k);
    data = array2table(Y);

    treatment = 1;
    horizons = 0:3;
    model = LP(data, treatment, p, horizons, 'includeConstant', true);

    % testing if constant was correctly defined
    assert(all(model.X(:, 1) == 1));

    % test the rest of the matrices
    for i = 1:length(horizons)
        h = horizons(i);

        % test treatment variables
        assert(isequal(model.X(1:end-h, 2:(treatment+1)), model.Y(1:end-h, 1:treatment, i) - h));

        % test lagged regressors
        for j = 1:p
            cols = (j-1)*k + 3 : j*k + 2;
            assert(isequal(model.X(1:end-h, cols), model.Y(1:end-h, :, i) - h - j));
        end
    end
end

function testLPBasicFunctions(testCase)
    k = 3; 
    p = 2;
    T = 1000;

    Y = randn(T, k);
    treatment = 1;
    horizons = 0:4;
    model = LP(Y, treatment, p, horizons, 'includeConstant', true);


    assert(~model.isFitted());
    assert(all(model.nobs() == T - p - horizons));
    model.getDependent();
    model.getIndependent();
    model.getInputData();
    model.getVariableNames();
    model.isStructural();

    % throws errors because not fitted
    verifyError(testCase, @() model.coeffs(), "LP:NotFitted");
    verifyError(testCase, @() model.fitted(), "LP:NotFitted");
    verifyError(testCase, @() model.residuals(), "LP:NotFitted");

    method = Recursive();
    model.fit(method);
    model.fit();
    model.coeffs();
    model.fitted();
    model.residuals();
end

function testLPRecursive(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 1000000;
    trendExponents = [0];
    m = length(trendExponents);

    A0 = tril(randn(k, k));
    S = diag(sign(diag(A0)));
    A0 = A0 * S;
    B = 0.2 * randn(k, k*p + m);
    APlus = A0 * B;

    maxHorizon = 4;
    Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);
    irfsTrue = SVAR.IRF_(A0, APlus(:, (m+1):end), p, maxHorizon);

    for treatment = 1:k
        model = LP(Y, treatment, p+10, 0:maxHorizon, 'includeConstant', true);
        method = Recursive();
        model.fit(method);

        irfObj = model.IRF(maxHorizon);
        irfsLP = irfObj.irfs(:, 1, :);

        testDiff = irfsLP - irfsTrue(:, treatment, :) ./ irfsTrue(treatment, treatment, 1);
        assert(all(max(abs(testDiff), [], 'all') < 1e-2));
    end
end

function testLPInformationCriteria(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 10000;
    trendExponents = [0];
    m = length(trendExponents);

    A0 = tril(randn(k, k));
    S = diag(sign(diag(A0)));
    A0 = A0 * S;
    B = 0.2 * randn(k, k*p + m);
    APlus = A0 * B;

    Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);

    modelLarge = LP(Y, 1, p+10, 0:4);
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

function testLPExternalInstrument(testCase)
    rng(6150533);
    k = 3;
    p = 2;
    T = 1000000;
    trendExponents = [0];
    m = length(trendExponents);

    A0 = tril(0.1 * randn(k, k));
    S = diag(sign(diag(A0)));
    A0 = A0 * S;
    B = 0.4 * randn(k, k*p + m);
    APlus = A0 * B;

    shocks = randn(k, T);
    Y = SVAR.simulate(shocks, A0, APlus, 'trendExponents', trendExponents);
    data = [shocks(1, :)' Y];
    data = array2table(data, 'VariableNames', {'instrument', 'Y1', 'Y2', 'Y3'});

    maxHorizon = 4;
    irfsTrue = SVAR.IRF_(A0, APlus(:, (m+1):end), p, maxHorizon);
    irfsTrue = irfsTrue(:, 1, :) ./ irfsTrue(1, 1, 1);

    % Defining treatment and instrument by number
    treatment = 2;
    model = LP(data, treatment, p, 0:maxHorizon);
    method = ExternalInstrument(2, 1);
    model.fit(method);
    irfObj = model.IRF(maxHorizon);
    irfsLP = irfObj.irfs(2:end, :, :);

    testDiff = irfsLP - irfsTrue;
    assert(all(max(abs(testDiff), [], 'all') < 1e-2));

    % Defining treatment and instrument by name
    treatment = 'Y1';
    model = LP(data, treatment, p, 0:maxHorizon);
    method = ExternalInstrument(treatment, {'instrument'});
    model.fit(method);
    irfObj = model.IRF(maxHorizon);
    irfsLPName = irfObj.irfs(2:end, :, :);

    testDiff = irfsLPName - irfsLP;
    assert(all(max(abs(testDiff), [], 'all') < sqrt(eps())));

    testDiff = irfsLPName - irfsTrue;
    assert(all(max(abs(testDiff), [], 'all') < 1e-2));

    % creating multiple instruments
    data.instrument2 = data.instrument + 0.1 * randn(T, 1);
    data.instrument3 = data.instrument + 0.1 * randn(T, 1);
    data = data(:, {'instrument', 'instrument2', 'instrument3', 'Y1', 'Y2', 'Y3'});

    treatment = 'Y1';
    model = LP(data, treatment, p, 0:maxHorizon);
    method = ExternalInstrument(treatment, {'instrument', 'instrument2', 'instrument3'});
    model.fit(method);
    irfObj = model.IRF(maxHorizon);
    irfsLPName = irfObj.irfs(4:end, :, :);

    testDiff = irfsLPName - irfsTrue;
    assert(all(max(abs(testDiff), [], 'all') < 1e-2));

    % Changing the normalising horizon
    treatment = 'Y2'; 
    % Changing shock because effect of shock1 on Y1 close to zero for h=1
    data.instrument = shocks(2, :)';
    data.instrument2 = data.instrument + 0.1 * randn(T, 1);
    data.instrument3 = data.instrument + 0.1 * randn(T, 1);
    model = LP(data, treatment, p, 0:maxHorizon);
    method = ExternalInstrument(treatment, {'instrument', 'instrument2', 'instrument3'}, 'normalisingHorizon', 1);
    model.fit(method);
    irfObj = model.IRF(maxHorizon);
    irfsLPName = irfObj.irfs(4:end, :, :);
    % Need to adjust the true IRFs too
    irfsTrue = SVAR.IRF_(A0, APlus(:, (m+1):end), p, maxHorizon);
    irfsTrue = irfsTrue(:, 2, :) ./ irfsTrue(2, 2, 2);

    testDiff = irfsLPName - irfsTrue;
    assert(all(max(abs(testDiff), [], 'all') < 1e-2));
end

