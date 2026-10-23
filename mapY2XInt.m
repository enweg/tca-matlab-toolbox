function xIdx = mapY2XInt(i, t, K, order)
    % `mapY2XInt` Map a variable and horizon of the dynamic system to the index of a systems form variable `x`.
    %
    %   `xIdx = mapY2XInt(i, t, K, order)` returns the index of the systems
    %   form variable corresponding to `y_{i, t}`.
    %
    %   ## Arguments
    %   - `i` (integer): Index of the variable in the original ordering.
    %   - `t` (integer): Horizon of the variable. The contemporaneous period is
    %     horizon `0`.
    %   - `K` (integer): Number of variables in the system.
    %   - `order` (vector): Ordering of variables given by the transmission matrix.
    %
    %   ## Returns
    %   - `xIdx` (integer): Index of the corresponding variable in the systems form.
    %
    %   See also `mapX2YInt`, `mapY2X`

    xIdx = K * t + find(order == i, 1);
end
