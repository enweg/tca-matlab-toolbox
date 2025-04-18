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

% TODO: implement the remaining test
