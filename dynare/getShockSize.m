function [shockSize, idx] = getShockSize(M_, shockName)
    shocks = dynareCellArrayToVec(M_.exo_names);
    idx = find(shocks == shockName);
    shockSize = sqrt(M_.Sigma_e(idx, idx));
end
