function P = permmatrix(order)
    % `permmatrix` Create a permutation matrix.
    %
    %   `P = permmatrix(order)` returns the permutation matrix whose rows are
    %   the rows of the identity matrix in the order given by `order`, so
    %   that `P * x` equals `x(order)`.
    %
    %   ## Arguments
    %   - `order` (vector): Ordering of variables given by the transmission matrix.
    %
    %   ## Returns
    %   - `P` (matrix): Permutation matrix of size `length(order)`.
    %
    %   ## Notes
    %   - This function is intended for internal use only.
    %
    %   See also `makeB`, `makeOmega`

    P = eye(length(order));
    P = P(order, :);
end
