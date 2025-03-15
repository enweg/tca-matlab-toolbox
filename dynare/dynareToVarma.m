function [Phi0, As, Psis, p, q] = dynareToVarma(M_, oo_, options_, maxKappa)
    % `dynareToVarma` Transform a DSGE model into a VARMA representation.
    %
    %   `[Phi0, As, Psis, p, q] = dynareToVarma(M_, oo_, options_, maxKappa)`
    %   converts a linearized DSGE model estimated using Dynare into a VARMA form,
    %   following the method of Morris (2016).
    %
    %   ## Arguments
    %   - `M_` (struct): Returned by Dynare. 
    %   - `oo_` (struct): Returned by Dynare. 
    %   - `options_` (struct): Returned by Dynare. 
    %   - `maxKappa` (integer, optional): Tuning parameter in Morris (2016). 
    %      Related to the maximum AR order via `maxArOrder=maxKappa+1`.
    %      Defaults to 20. 
    %
    %   ## Returns
    %   - `Phi0` (matrix): Impact matrix linking structural shocks to reduced-form errors.
    %   - `As` (cell array): AR coefficient matrices `{A_1, A_2, ..., A_p}`.
    %   - `Psis` (cell array): MA coefficient matrices `{Psi_1, Psi_2, ..., Psi_q}`.
    %   - `p` (integer): The determined autoregressive order of the VARMA representation.
    %   - `q` (integer): The determined moving average order of the VARMA representation.
    %
    %   ## Methodology
    %   The function follows the approach outlined in Morris (2016) and returns
    %   a VARMA of the form: 
    %   $$
    %   y_t = \sum_{i=1}^{p} A_i y_{t-i} + \sum_{j=1}^{q} \Psi_j u_{t-j} + u_t,
    %   $$
    %   where:
    %   - $u_t = \Phi_0 \varepsilon_t$, with $\varepsilon_t$ being structural shocks.
    %
    %   ## Reference
    %   - Morris, S. D. (2016). "VARMA representation of DSGE models." *Economics Letters*, 138, 30–33.
    %     [https://doi.org/10.1016/j.econlet.2015.11.027](https://doi.org/10.1016/j.econlet.2015.11.027)
    %
    %   See also `getABCD`, `varmaIrfs`.

    if ~isfield(options_, "varobs")
        error("dynareToVarma: No observed variables were defined in the mod file.")
    end

    % Default choice for maximum VAR order following notation in Morris 2016. 
    if nargin==3
        maxKappa = 20;
    end

    [A, B, C, D] = getABCD(M_, oo_, options_);

    % Basic assumption is that D is invertible. 
    condTolerance = 1e-20;
    if rcond(D) < condTolerance
        error("Matrix D must not be singular.")
    end

    n = size(C, 1); 
    m = size(A, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Case I: C is invertible. In that case, p=q=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if n==m && rcond(C) > condTolerance
        p = 1;
        q = 1;

        % Finding AR, MA matrices
        CInv = inv(C);
        DInv = inv(D);
        Phi0 = D;
        As = {C*A*CInv};
        Psis = {C*(B - A*CInv*D)*DInv};
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Case II: Follow the general proposition of Morris 2016
    % but adjusted for our notation.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    kappa = 0;
    F = [C];
    while kappa < maxKappa 
        kappa = kappa + 1;
        F = [C*A^kappa; F];
        if rank(F) == size(F, 2)
            FPlus = pinv(F);  % Checked rank condition above.
            if rank(FPlus(:, 1:n)) == n
                break  % Both rank conditions are satisfied.
            end
        end
    end
    if rank(F) ~= size(F, 2)
        error("dynareToVarma: Rank condition for F is not satisfied. Could not find VARMA representation of DSGE.")
    end
    FPlus = pinv(F);  % Checked rank condition above
    if rank(FPlus(:, 1:n)) ~= n
        error("dynareToVarma: Rank condition for FPlus is not satisfied. Could not find VARMA representation of DSGE.")
    end
    p = kappa + 1;
    q = kappa + 1;

    % Finding VARMA Coefficients
    Theta = FPlus(:, 1:n);
    ThetaPlus = pinv(Theta);  % Checked rank conditions above.

    G = cell(1, kappa+1);
    for col=1:(kappa+1)
        Gcol = zeros(n*(kappa+1), n);
        for row=1:(col-1)
            rowStart = (row-1)*n+1;
            rowEnd = row*n;
            Gcol(rowStart:rowEnd, :) = C*A^(col-row-1)*B;
        end
        rowStart = (col-1)*n+1;
        rowEnd = col*n;
        Gcol(rowStart:rowEnd, :) = D;
        G{col} = Gcol;
    end

    As = cell(1, kappa+1);
    Phi0 = ThetaPlus * FPlus * G{1};
    A0 = inv(Phi0);  
    for k=1:kappa
        colStart = (k-1)*n + 1;
        colEnd = k*n;
        As{k} = ThetaPlus * (A*FPlus(:, colStart:colEnd) - FPlus(:, (colStart+n):(colEnd+n)));
    end
    colStart = (kappa+1-1)*n + 1;
    colEnd = (kappa+1)*n;
    As{end} = ThetaPlus * A * FPlus(:, colStart:colEnd);

    Psis = cell(1, kappa+1);
    for k=1:kappa 
        Psis{k} = ThetaPlus * (FPlus * G{k+1} - A * FPlus * G{k}) * A0;
    end
    Psis{end} = ThetaPlus * (B - A * FPlus * G{kappa+1}) * A0;
end
