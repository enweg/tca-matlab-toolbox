function idx = getVariableIdx(varname, options_)
    varnames = dynareCellArrayToVec(options_.varobs);
    idx = find(varnames == varname);
end
