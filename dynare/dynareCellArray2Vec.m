function v=dynareCellArray2Vec(ca)
    % `dynareCellArray2Vec` Convert a Dynare cell array into a standard vector.
    %
    %   `v = dynareCellArray2Vec(ca)` converts a cell array of strings, typically 
    %   used by Dynare for storing names or labels, into a standard string vector.
    %
    %   ## Arguments
    %   - `ca` (Cell array): A cell array where each cell contains a string.
    %
    %   ## Return
    %   - `v` (Vector): A vector of strings, converted from the input cell array.
    v = repelem("", length(ca));
    for i=1:length(ca)
        v(i) = string(cell2mat(ca(i)));
    end
    v = vec(v);
end
