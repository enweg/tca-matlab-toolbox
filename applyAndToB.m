function [B, Omega] = applyAndToB(B, Omega, from, var)
    Omega((var+1):end, from) = 0;
    B((var+1):end, 1:(var-1)) = 0;
end
