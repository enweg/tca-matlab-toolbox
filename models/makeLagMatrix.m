function lagMatrix = makeLagMatrix(Y, p)
    % `makeLagMatrix` Create lag matrix from input data.
    %
    %   `lagMatrix = makeLagMatrix(Y, p)` returns a matrix where
    %   columns contain 1 to `p` lags of `Y`.
    %
    %   ## Arguments
    %   - `Y` (matrix): Input data matrix (T x k).
    %   - `p` (integer): Number of lags.
    %
    %   ## Returns
    %   - `lagMatrix` (matrix): Matrix of lagged values. First `k`
    %     columns are lag 1, next `k` are lag 2, and so on.
    %
    %   ## Notes
    %   - Missing values due to lagging are filled with NaN.
    [T, k] = size(Y);
    lagMatrix = nan(T, k * p);
    for l = 1:p
        lagMatrix((l+1):end, ((l-1)*k+1):(l*k)) = Y(1:(end-l), :);
    end
end
