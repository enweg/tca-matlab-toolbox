function [p, q]=determineVarmaOrder(M_, options_, oo_, maxKappa, obsVar)
    if nargin==4
        [Ax, Ay, B, Cx, Cy, D, vars] = getAbcdVarma(M_, options_, oo_);
    elseif nargin==5
        [Ax, Ay, B, Cx, Cy, D, vars] = getAbcdVarma(M_, options_, oo_, obsVar);
    end

    % Main assumption
    if rank(D) < size(D, 1)
        error("D is not invertible.");
    end

    % Corollary 1 and 2
    m = size(Ax, 1);
    n = size(Cx, 1);
    if n == m && rank(Cx) == size(Cx, 1)
        p = 2;
        q = 1;
        if all(abs(Ay) < 1e-10) & all(abs(Cy) < 1e-10)
            p = 1;
            q = 1;
        end
        return; 
    end

    % Proposition 1
    p = Inf;
    q = Inf; 
    C = [Cx Cy];
    A = [[Ax Ay]; [Cx Cy]];
    FKappa = Cx;
    for kappa=1:maxKappa
        FKappa = [(C*A^(kappa-1)*[Ax' Cx']')' FKappa'];
        FKappa = FKappa';
        FKappaPlus = pinv(FKappa);
        PsiKappa = FKappaPlus(:, 1:n);
        if rank(PsiKappa) == size(PsiKappa, 2)
            p = kappa + 2;
            q = kappa + 1;
            return;
        end
    end

    error("Could not find a VARMA representation with maxKappa="+maxKappa)
end




