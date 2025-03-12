function [B, Omega] = makeSystemsForm(Phi0, As, Psis, order, maxHorizon)
    % `makeSystemsForm` Transform an SVARMA dynamic model into the system representation $x = Bx + \Omega\varepsilon$.
    %
    %   ## Arguments
    %   - `Phi0` (matrix): The matrix of contemporaneous structural impulse responses.
    %   - `As` (cell array of matrices): A vector of reduced-form autoregressive (AR) 
    %     matrices, where the first entry corresponds to the AR matrix for the first lag, etc.
    %   - `Psis` (cell array of matrices): A vector of reduced-form moving average (MA) 
    %     matrices, where the first entry corresponds to the MA matrix for the first lag, etc.
    %   - `order` (vector): The vector of intergers indicating the order of 
    %     variables, typically determined by the transmission matrix.
    %   - `maxHorizon` (integer): The maximum time horizon to consider for the systems form, 
    %     with `0` representing the contemporaneous period.
    %
    %   ## Returns
    %   - `B` (matrix)
    %   - `Omega` (matrix)
    %
    %   See also `makeB`, `makeOmega`.

    % 1. Create B and Omega using helper functions
    Sigma = Phi0 * Phi0';
    B = makeB(As, Sigma, order, maxHorizon);
    Omega = makeOmega(Phi0, Psis, order, maxHorizon);
end
