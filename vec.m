function v = vec(m)
    % `vec` Stack the columns of an array into a column vector.
    %
    %   `v = vec(m)` returns `m(:)`.
    %
    %   ## Arguments
    %   - `m` (array): Any array.
    %
    %   ## Returns
    %   - `v` (vector): Column vector with all elements of `m`.

    v = m(:);
end
