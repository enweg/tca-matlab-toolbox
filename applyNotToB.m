function [BTilde, OmegaTilde] = applyNotToB(B, Omega, from, var)
    % `applyNotToB` Manipulate B and Omega to ensure var lies on no paths.
    %
    %   `applyNotToB(B, Omega, from, var)` modifies the matrices `B` and `Omega` such that 
    %   the variable indexed by `var` is absent from all paths. This is achieved by:
    %   - Zeroing out the edge from the shock (indexed by from) to `var`.
    %   - Zeroing out all edges from variables ordered before var to `var`.
    %
    %   ## Arguments
    %   - `B` (matrix): Part of the systems form.
    %   - `Omega` (matrix): Part of the systems form.
    %   - `from` (integer): The shock index.
    %   - `var` (integer): The variable index that cannot lie on any paths.
    %
    %   ## Notes
    %   - This function is intended for internal use only.
    %   
    %   See also `applyAndToB`, `makeB`, `makeOmega`, and `makeSystemsForms`.

    OmegaTilde = Omega; 
    BTilde = B;
    OmegaTilde(var, from) = 0;
    BTilde(var, 1:(var-1)) = 0;
end
