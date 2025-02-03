function irfs=irf_varma(A0inv, Bs, Psis, horizon)
    % `varmaIrf` Compute the Impulse Response Functions (IRFs) for a structural VARMA model.
    %
    %   `varmaIrf(A0, Bs, Psis, horizon)` computes the impulse response functions 
    %   of a structural Vector Autoregressive Moving Average (VARMA) model over a 
    %   specified time horizon.
    %
    %   ## Arguments
    %   - `A0inv` (matrix, n × n): Structural impact effect matrix, where `n` is the number of variables.
    %   - `Bs` (Cell array of length p): Array of reduced form AR matrices (n × n) where p is the AR order.
    %   - `Psis` (Cell array of length q): Array of reduced form MA matrices (n × n) where q is the MA order.
    %   - `horizon` (integer): Number of periods for which to compute the IRFs.
    %
    %   ## Returns
    %   - `irfs` (3D array, n × n × (horizon+1)): Impulse response functions. The dimensions 
    %     correspond to the number of variables, the number of shocks, and the number of periods, 
    %     respectively.
    %
    
    p = size(Bs, 3);
    q = size(Psis, 3);
    n = size(A0inv, 1);

    % calculating irfs
    irfs = zeros(n, n, horizon+1);
    irfs(:, :, 1) = A0inv;
    for h=1:horizon
        for i=1:min(p, h)
            irfs(:, :, h+1) = irfs(:, :, h+1) + Bs{i}*irfs(:, :, h-i+1);
        end
        if h <= q
            irfs(:, :, h+1) = irfs(:, :, h+1) + Psis{h} * A0inv;
        end
    end
end
