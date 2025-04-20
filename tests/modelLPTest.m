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

% TODO: implement the remaining test
