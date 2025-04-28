function leadMatrix = makeLeadMatrix(Y, p)
    % `makeLeadMatrix` Create lead matrix from input data.
    %
    %   `leadMatrix = makeLeadMatrix(Y, p)` returns a matrix where
    %   columns contain 1 to `p` leads of `Y`.
    %
    %   ## Arguments
    %   - `Y` (matrix): Input data matrix (T x k).
    %   - `p` (integer): Number of leads.
    %
    %   ## Returns
    %   - `leadMatrix` (matrix): Matrix of lead values. First `k`
    %     columns are lead 1, next `k` are lead 2, and so on.
    %
    %   ## Notes
    %   - Missing values due to leading are filled with NaN.
    [T, k] = size(Y);
    leadMatrix = nan(T, k * p);
    for l = 1:p
        leadMatrix(1:(end-l), ((l-1)*k+1):(l*k)) = Y((l+1):end, :);
    end
end
