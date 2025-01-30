function [BTilde, OmegaTilde] = applyNotToB(B, Omega, from, var)
    OmegaTilde = Omega; 
    BTilde = B;
    OmegaTilde(var, from) = 0;
    BTilde(var, 1:(var-1)) = 0;
end
