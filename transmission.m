function effects = transmission(from, arr1, arr2, q, method)
    [varAnd, varNot, multiplier] = getVarNumsAndMultiplier(q);
    if isequal(method, "BOmega")
        disp('INFO: using method BOmega to compute transmission effect.');
        % apply the BOmega method
        % arr1 = B, arr2 = Omega
        effects = transmissionBOmega(from, arr1, arr2, varAnd, varNot, multiplier);
    elseif isequal(method, "irf")
        disp('INFO: using method irfs to compute transmission effect.');
        % apply the irf method 
        % arr1 = irfs, arr2 = irfsOrtho
        effectsCell = cell(1, length(varAnd)); 
        parfor ii = 1:length(varAnd) % Use parfor for parallel execution
            vAnd = varAnd{ii};
            vNot = varNot{ii};
            m = multiplier(ii);
            effectsCell{ii} = transmissionIrfs(from, arr1, arr2, vAnd, vNot, m);
        end
        effects = effectsCell{1};
        for ii = 2:length(effectsCell)
            effects = effects + effectsCell{ii};
        end
    else
        error("Unrecognised method. Must be either BOmega or irf");
    end
end
