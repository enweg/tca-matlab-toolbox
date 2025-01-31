function [B, Omega] = applyAndToB(B, Omega, from, var)
    % `applyAndToB` Manipulate `B` and `Omega` to ensure var lies on all paths.
    %
    %   `applyAndToB(B, Omega, from, var)` modifies the matrices `B` and `Omega` such that 
    %   the variable indexed by `var` is present on all paths. This is achieved by:
    %   - Zeroing out all edges going directly from the shock (indexed by from) 
    %     to any variables ordered after `var`.
    %   - Zeroing out any edges going from variables ordered before `var` to any 
    %     variables ordered after `var`.
    %
    %   ## Arguments
    %   - `B` (matrix): Part of the system's representation.
    %   - `Omega` (matrix): Part of the system's representation.
    %   - `from` (integer): The shock index.
    %   - `var` (integer): The variable index that must lie on all paths.
    %
    %   ## Notes
    %   - This function is intended for internal use only.
    %
    %   See also `applyNotToB`, `makeB`, `makeOmega`, and `makeSystemsForm`.

    Omega((var+1):end, from) = 0;
    B((var+1):end, 1:(var-1)) = 0;
end
