function irfs=varmaIrfs(Phi0, As, Psis, horizon)
    % `varmaIrfs` Compute structural impulse response functions (IRFs) for a VARMA model.
    %
    %   `irfs = varmaIrfs(Phi0, As, Psis, horizon)` computes the structural impulse
    %   response functions (IRFs) of a VARMA model, given the structural shock impact
    %   matrix, autoregressive (AR) coefficients, and moving average (MA) coefficients.
    %
    %   ## Model Specification
    %   The VARMA model is defined as:
    %   $$
    %   y_t = \sum_{i=1}^{p} A_i y_{t-i} + \sum_{j=1}^{q} \Psi_j u_{t-j} + u_t,
    %   $$
    %   where:
    %   - $u_t = \Phi_0 \varepsilon_t$, with $\varepsilon_t$ being structural shocks.
    %
    %   ## Arguments
    %   - `Phi0` (matrix): Impact matrix linking structural shocks to reduced-form errors.
    %   - `As` (cell array): AR coefficient matrices `{A_1, A_2, ..., A_p}`.
    %   - `Psis` (cell array): MA coefficient matrices `{Psi_1, Psi_2, ..., Psi_q}`.
    %   - `horizon` (integer): Number of periods for which IRFs are computed.
    %      `horizon=0` means only contemporaneous impulses are computed which are 
    %      the same as `Phi0`.
    %
    %   ## Returns
    %   - `irfs` (3D array): Structural IRFs of size `(n, m, horizon+1)`, where `n`
    %     is the number of endogenous variables, `m` is the number of shocks, 
    %     assumed to satisfy `m=n`. The IRFs capture the dynamic response
    %     of each variable to a unit shock over the specified horizon.
    %

    p = length(As);
    q = length(Psis);
    n = size(Phi0, 1);

    % calculating irfs
    irfs = zeros(n, n, horizon+1);
    irfs(:, :, 1) = Phi0;
    for h=1:horizon
        for i=1:min(p, h)
            irfs(:, :, h+1) = irfs(:, :, h+1) + As{i}*irfs(:, :, h-i+1);
        end
        if h <= q
            irfs(:, :, h+1) = irfs(:, :, h+1) + Psis{h} * Phi0;
        end
    end
end
