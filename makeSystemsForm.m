function [B, Omega] = makeSystemsForm(Phi0, As, Psis, Sigma, order, maxHorizon)
    % 1. Create B and Omega using helper functions
    B = makeB(As, Sigma, order, maxHorizon);
    Omega = makeOmega(Phi0, Psis, Sigma, order, maxHorizon);
end
