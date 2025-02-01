function [A0, Phis, Psis, vars]=getVarmaCoeffs(M_, options_, oo_, p, q, obsVar)
    if nargin==5
        [Ax, Ay, Bx, Cx, Cy, D, vars] = getAbcdVarma(M_, options_, oo_);
    elseif nargin==6
        [Ax, Ay, Bx, Cx, Cy, D, vars] = getAbcdVarma(M_, options_, oo_, obsVar);
    end
    
    n = size(Cx, 1);
    m = size(Ax, 1);
    Phis = zeros(n, n, p);
    Psis = zeros(n, n, q);

    vars = vars((end - n + 1):end);
    vars = getVariableOrdering(M_, options_, oo_, vars);

    % Main assumption
    if rank(D) < size(D, 1)
        error("D is not invertible.");
    end

    if p<=2 && q<=1
        % case of corollary 1 or 2
        % double checking rank of Cx
        if rank(Cx) ~= size(Cx, 1)
            error("Cannot have p<=2 and q<=1 if Cx is not full rank.")
        end
        Phis(:, :, 1) = D \ ((Cx * Ax) / Cx + Cy);
        if p==2
            Phis(:, :, 2) = D \ (Ay - ((Cx * Ax) / Cx) * Cy);
        end
        Psis(:, :, 1) = D \ ((Cx * Bx / D - Cx * Ax / Cx)*D);
        A0 = D \ eye(n);
        return;
    end

    % case of general proposition 1
    if p - 1 ~= q
        error("Proposition states that p-1=q. This is not satisfied.");
    end
    kappa = p-2;

    C = [Cx Cy];
    A = [Ax Ay; Cx Cy];
    B = [Bx; D];

    F = C; 
    for k=1:kappa
        F = [C*A^k; F];
    end
    Fx = F(:, 1:m);
    FxPlus = pinv(Fx);
    Fy = F(:, (m+1):end);
    
    Gs = zeros(n*(kappa+1), n, kappa+1);
    Gi = D;
    Gs(:, :, 1) = [Gi; zeros(n*kappa, n)];
    Gi = [C*B; Gi];
    Gs(:, :, 2) = [Gi; zeros(n*(kappa-1), n)];
    for i=3:(kappa+1)
        Gi = [C*A^(i-2)*B; Gi];
        Gs(:, :, i) = [Gi, zeros(n*(kappa+1-i), n)];
    end
      
    coefEt = (FxPlus * Gs(:, :, 1));
    pinvCoefEt = pinv(coefEt);
    A0 =  pinvCoefEt * FxPlus(:, 1:n); 
    Phis(:, :, end) = pinvCoefEt *  (Ay - Ax * FxPlus * Fy);
    Phis(:, :, end-1) = pinvCoefEt * (Ax * FxPlus(:, (end-n+1):end) + FxPlus * Fy);
    for k=1:kappa
        Phis(:, :, kappa) = pinvCoefEt * (Ax * FxPlus(:, ((k-1)*n+1):(k*n)) - FxPlus(:, (k*n+1):((k+1)*n)));
    end
    
    Psis(:, :, 1) = pinvCoefEt * (Bx - Ax * FxPlus*Gs(:, :, end));
    for k=1:kappa
        Psis(:, :, k) = pinvCoefEt * (FxPlus * Gs(:, :, k+1) - Ax * FxPlus * Gs(:, :, k));
    end
end


