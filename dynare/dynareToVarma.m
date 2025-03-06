function [Phi0, As, Psis, p, q] = dynareToVarma(M_, oo_, options_, maxKappa)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x_t = Ax_{t-1} + B\varepsilon_t
    % y_t = Cx_{t-1} + D\varepsilon_t
    % x_t is the state vector, y_t is the observed vector (defined under the 
    % varobs block in dynare mod file)
    %
    % Output model is a SVARMA
    % A_0y_t = A_1y_{t-1} + ... + A_py_{t-p} + \varepsilon_t 
    %        + Psi_1\varepsilon_{t-1} + ... + Psi_{t-q}\varepsilon_{t-q}
    % As = {A_0, A_1, ...}
    % Psis = {Psi_0, Psi_1, ...}
    %
    %
    % IMPORTANT: original version retuned SVARMA from above. This version returns 
    % VARMA version plus Phi0 -- the contemporaneous impact matrix. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~isfield(options_, "varobs")
        error("dynareToVarma: No observed variables were defined in the mod file.")
    end

    % Dynare handles shock sizes seperately. Let's get these first. 
    S = sqrt(M_.Sigma_e);
    % S = eye(size(M_.Sigma_e, 1));
    SInv = inv(S);

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
        Phi0 = D * S;
        As = {C*A*CInv};
        Psis = {C*(B - A*CInv*D)*SInv*DInv};
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Case II: Follow the general proposition of Morris 2016
    % but adjusted for our notation.
    %
    % Morris, S. D. (2016). VARMA representation of DSGE models. 
    % Economics Letters, 138, 30–33. https://doi.org/10.1016/j.econlet.2015.11.027
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

    % Premultiplication of As{1} below is because we form the structural VARMA 
    % representation and thus need to pre-multiply all non-structural VARMA 
    % matrices by As{1} = A_0.
    As = cell(1, kappa+1);
    Phi0 = ThetaPlus * FPlus * G{1};
    Phi0 = Phi0 * S;  % because of how dynare handles shock sizes
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
        Psis{k} = ThetaPlus * (FPlus * G{k+1} - A * FPlus * G{k}) * S * A0;
    end
    Psis{end} = ThetaPlus * (B - A * FPlus * G{kappa+1}) * S * A0;
end
