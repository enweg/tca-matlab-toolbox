function tests = simplifyingTest
    tests = functiontests(localfunctions);
end

function testSystemsForm(testCase)
    % Structural model
    A0 = randn(3, 3);
    A0inv = inv(A0);
    As = arrayfun(@(x) randn(3, 3), 1:2, 'UniformOutput', false);
    Psis = arrayfun(@(x) randn(3, 3), 1:2, 'UniformOutput', false);
    Sigma = A0inv * A0inv';

    [L, D] = makeLD(Sigma);
    O = zeros(3, 3);
    Q = A0 * inv(L);

    B = [
        (eye(3) - D*L), O, O, O;
        D*Q'*As{1}, (eye(3) - D*L), O, O;
        D*Q'*As{2}, D*Q'*As{1}, (eye(3) - D*L), O;
        O, D*Q'*As{2}, D*Q'*As{1}, (eye(3) - D*L)
    ];

    Omega = [
        D*Q', O, O, O;
        D*Q'*Psis{1}, D*Q', O, O;
        D*Q'*Psis{2}, D*Q'*Psis{1}, D*Q', O;
        O, D*Q'*Psis{2}, D*Q'*Psis{1}, D*Q'
    ];

    % Reduced form model
    As = cellfun(@(A) A0inv * A, As, 'UniformOutput', false);
    Psis = cellfun(@(Psi) A0inv * Psi * A0, Psis, 'UniformOutput', false);

    % testing
    BTest = makeB(As, Sigma, 1:3, 3);
    assert(all(max(abs(B - BTest), [], 'all') < sqrt(eps())));

    BTest = makeB(As, Sigma, 1:3, 1);
    assert(all(max(abs(BTest - B(1:6, 1:6)), [], 'all') < sqrt(eps())));

    OmegaTest = makeOmega(A0inv, Psis, Sigma, 1:3, 3);
    assert(all(max(abs(OmegaTest - Omega), [], 'all') < sqrt(eps())));

    OmegaTest = makeOmega(A0inv, Psis, Sigma, 1:3, 1);
    assert(all(max(abs(OmegaTest - Omega(1:6, 1:6)), [], 'all') < sqrt(eps())));
end

