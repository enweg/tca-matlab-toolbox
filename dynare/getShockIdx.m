function idxShock = getShockIdx(shockName, M_)
    shocks = dynareCellArrayToVec(M_.exo_names);
    idxShock = find(shocks == shockName);
end
