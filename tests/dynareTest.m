function tests = dynareTests
    tests = functiontests(localfunctions);
end

function testDynareVarmaSW2007(testCase)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETTING UP THE TEST
    % Assumes that the model was run once so that the output folder exists
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load ./tests/SW2007/SW2007/Output/SW2007_results.mat
    vars = struct2cell(oo_.irfs);
    names = fieldnames(oo_.irfs);
    % Dynamically create variables
    for i = 1:numel(names)
        % Assign in the current function workspace
        eval([names{i} ' = transpose(vars{i});']); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END SETTING UP THE TEST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Phi0, As, Psis, p, q] = dynareToVarma(M_, oo_, options_);
    irfsVarma = varmaIrfs(Phi0, As, Psis, 19);

    tol = 1e-9;  % accepted tolerance throughout all tests
    n = size(Phi0, 1);

    irfsDSGEea = [dy_ea dc_ea dinve_ea pinfobs_ea dw_ea robs_ea labobs_ea];
    irfsDSGEea = reshape(irfsDSGEea', n, 1, []);
    [shockSize, idx] = getShockSize(M_, "ea");  % only need this for idx 
    testea = max(vec(abs(irfsDSGEea - irfsVarma(:, idx, :))));
    assert(testea < tol);

    irfsDSGEeb = [dy_eb dc_eb dinve_eb pinfobs_eb dw_eb robs_eb labobs_eb];
    irfsDSGEeb = reshape(irfsDSGEeb', n, 1, []);
    [shock_size, idx] = getShockSize(M_, "eb");
    testeb = max(vec(abs(irfsDSGEeb - irfsVarma(:, idx, :))));
    assert(testeb < tol);

    irfsDSGEeg = [dy_eg dc_eg dinve_eg pinfobs_eg dw_eg robs_eg labobs_eg];
    irfsDSGEeg = reshape(irfsDSGEeg', n, 1, []);
    [shock_size, idx] = getShockSize(M_, "eg");
    testeg = max(vec(abs(irfsDSGEeg - irfsVarma(:, idx, :))));
    assert(testeg < tol);

    irfsDSGEeqs = [dy_eqs dc_eqs dinve_eqs pinfobs_eqs dw_eqs robs_eqs labobs_eqs];
    irfsDSGEeqs = reshape(irfsDSGEeqs', n, 1, []);
    [shock_size, idx] = getShockSize(M_, "eqs");
    testeqs = max(vec(abs(irfsDSGEeqs - irfsVarma(:, idx, :))));
    assert(testeqs < tol);

    irfsDSGEem = [dy_em dc_em dinve_em pinfobs_em dw_em robs_em labobs_em];
    irfsDSGEem = reshape(irfsDSGEem', n, 1, []);
    [shock_size, idx] = getShockSize(M_, "em");
    testem = max(vec(abs(irfsDSGEem - irfsVarma(:, idx, :))));
    assert(testem < tol);

    irfsDSGEepinf = [dy_epinf dc_epinf dinve_epinf pinfobs_epinf dw_epinf robs_epinf labobs_epinf];
    irfsDSGEepinf = reshape(irfsDSGEepinf', n, 1, []);
    [shock_size, idx] = getShockSize(M_, "epinf");
    testepinf = max(vec(abs(irfsDSGEepinf - irfsVarma(:, idx, :))));
    assert(testepinf < tol);

    irfsDSGEew = [dy_ew dc_ew dinve_ew pinfobs_ew dw_ew robs_ew labobs_ew];
    irfsDSGEew = reshape(irfsDSGEew', n, 1, []);
    [shock_size, idx] = getShockSize(M_, "ew");
    testew = max(vec(abs(irfsDSGEew - irfsVarma(:, idx, :))));
    assert(testew < tol);
end

function testDynareVarmaGali2015(testCase)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETTING UP THE TEST
    % Assumes that the model was run once so that the output folder exists
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load ./tests/Gali2015/Gali2015Chapter3/Output/Gali2015Chapter3_results.mat
    vars = struct2cell(oo_.irfs);
    names = fieldnames(oo_.irfs);
    % Dynamically create variables
    for i = 1:numel(names)
        % Assign in the current function workspace
        eval([names{i} ' = transpose(vars{i});']); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END SETTING UP THE TEST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Phi0, As, Psis, p, q] = dynareToVarma(M_, oo_, options_);
    irfsVarma = varmaIrfs(Phi0, As, Psis, 19);

    tol = 1e-9;  % accepted tolerance throughout all tests
    n = size(Phi0, 1);

    irfsDSGEeps_nu = [pi_eps_nu, i_eps_nu, y_gap_eps_nu];
    irfsDSGEeps_nu = reshape(irfsDSGEeps_nu', n, 1, []);
    [shockSize, idx] = getShockSize(M_, "eps_nu");  % only need this for idx 
    testeps_nu = max(vec(abs(irfsDSGEeps_nu - irfsVarma(:, idx, :))));
    assert(testeps_nu < tol);

    irfsDSGEeps_z = [pi_eps_z, i_eps_z, y_gap_eps_z];
    irfsDSGEeps_z = reshape(irfsDSGEeps_z', n, 1, []);
    [shockSize, idx] = getShockSize(M_, "eps_z");  % only need this for idx 
    testeps_z = max(vec(abs(irfsDSGEeps_z - irfsVarma(:, idx, :))));
    assert(testeps_z < tol);

    irfsDSGEeps_a = [pi_eps_a, i_eps_a, y_gap_eps_a];
    irfsDSGEeps_a = reshape(irfsDSGEeps_a', n, 1, []);
    [shockSize, idx] = getShockSize(M_, "eps_a");  % only need this for idx 
    testeps_a = max(vec(abs(irfsDSGEeps_a - irfsVarma(:, idx, :))));
    assert(testeps_a < tol);
end
