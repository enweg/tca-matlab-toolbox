classdef DSGE < handle & Model
    % `DSGE` Dynamic Stochastic General Equilibrium (DSGE) model.
    %
    %   This class specifies a DSGE model structure. The model must
    %   be previously computed using Dynare. It provides access to
    %   the Dynare output structures.
    %
    %   ## Properties
    %   - `M_` (struct): Model structure returned by Dynare.
    %   - `options_` (struct): Options structure from Dynare.
    %   - `oo_` (struct): Output structure with results from Dynare.
    %
    %   ## Notes
    %   - The model must have already been solved in Dynare.
    %   - This class serves as a wrapper to interface with Dynare output.
    properties
        % The following are all structs returned by Dynare
        M_
        options_
        oo_
    end

    methods (Static)

        function checkDynare_()
            % `checkDynare_` Ensure Dynare is properly loaded and available.

            % First check if Dynare is loaded
            if exist("dynare", "file") == 2
                pathDynare = fileparts(which("dynare.m"));
                if isempty(pathDynare)
                    error("Could not determine Dynare installation path.");
                end

                if exist("kalman_transition_matrix", "file") == 2
                    % Older versions of Dynare have the function in the root directory.
                    % The function is therefore loaded as soon as base Dynare is loaded.
                    return;
                end

                % Newer versions of Dynare include the file as part of the
                % stochastic_solver subdirectory. The function is thus no-longer loaded
                % as soon as Dynare itself is added to the path. We must add
                % stochastic_solver manually to the path.

                pathStochasticSolver = fullfile(pathDynare, 'stochastic_solver');

                % Check if the folder exists
                if exist(pathStochasticSolver, 'dir')
                    addpath(pathStochasticSolver);
                    disp('[INFO]: Dynare exists and is ready.');
                else
                    error("Could not find 'kalman_transition_matrix' in Dynare functions.");
                end
            else
                error("Dynare does not exist in path. Please add Dynare to your path first.");
            end
        end

        function order = defineOrder_(vars, options_)
            % `defineOrder_` Determine the ordering of observed variables.
            %
            %   `defineOrder_(vars, options_)` returns an ordering vector `order`
            %   that maps the variables in `vars` to their corresponding positions
            %   in the list of observed variables of a DSGE model estimated using
            %   Dynare. It defines the transmission matrix and can be used in
            %   `makeB`, `makeOmega`, `makeSystemsForm`, `makeConditionY`,
            %   `notThrough`, `through`, and `transmission`.
            %
            %   ## Arguments
            %   - `vars` (vector): A list of observed variable names.
            %   - `options_` (struct): Options structure returned by Dynare.
            %
            %   ## Returns
            %   - `order` (vector): Indices of `vars` in the original observed
            %     variable list.
            %
            %   See also `transmission`, `through`, `notThrough`, `makeSystemsForm`.
            varsOriginal = DSGE.dynareCellArrayToVec_(options_.varobs);
            order = zeros(length(vars), 1);
            for ii = 1:length(vars)
                order(ii) = find(varsOriginal == vars(ii));
            end
        end

        function v = dynareCellArrayToVec_(ca)
          v = repelem("", length(ca));
          for i=1:length(ca)
            v(i) = string(cell2mat(ca(i)));
          end
          v = v(:);
        end

        function [Phi0, As, Psis, p, q] = dynareToVarma_(M_, oo_, options_, pMax)
            % `dynareToVarma_` Transform a DSGE model into VARMA representation.
            %
            %   `[Phi0, As, Psis, p, q] = dynareToVarma_(M_, oo_, options_, pMax)`
            %   converts a linearized DSGE model estimated using Dynare into a
            %   VARMA form, following the method of Morris (2016).
            %
            %   ## Arguments
            %   - `M_` (struct): Model structure returned by Dynare.
            %   - `oo_` (struct): Output structure returned by Dynare.
            %   - `options_` (struct): Options structure returned by Dynare.
            %   - `pMax` (integer, optional): Maximum AR order to try. 
            %
            %   ## Returns
            %   - `Phi0` (matrix): Impact matrix linking shocks to reduced-form errors.
            %   - `As` (cell array): AR coefficient matrices `{A_1, ..., A_p}`.
            %   - `Psis` (cell array): MA coefficient matrices `{Psi_1, ..., Psi_q}`.
            %   - `p` (integer): Determined autoregressive order.
            %   - `q` (integer): Determined moving average order.
            %
            %   ## Methodology
            %   The function follows the approach outlined in Morris (2016)
            %   and returns a VARMA of the form:
            %   $$
            %   y_t = \sum_{i=1}^{p} A_i y_{t-i} + \sum_{j=1}^{q} \Psi_j u_{t-j} + u_t,
            %   $$
            %   where:
            %   - $u_t = \Phi_0 \varepsilon_t$, with $\varepsilon_t$
            %     being structural shocks.
            %
            %   ## Reference
            %   - Morris, S. D. (2016). "VARMA representation of DSGE models."
            %     *Economics Letters*, 138, 30–33.
            %     [https://doi.org/10.1016/j.econlet.2015.11.027](https://doi.org/10.1016/j.econlet.2015.11.027)
            %
            %   See also `getABCD_`, `varmaIrfs_`.

            % Checking if Dynare is setup
            DSGE.checkDynare_();

            if ~isfield(options_, "varobs")
                error("dynareToVarma: No observed variables were defined in the mod file.")
            end

            % Default choice for maximum VAR order
            if nargin==3
                pMax = 20;
            end

            [A, B, C, D] = DSGE.getABCD_(M_, oo_, options_);
            % Note that this state-space form still allows for shocks that have
            % non-unity variance. Since we work with shocks that have unity
            % variance, we change below the shock related matrices B and D. 
            S = sqrt(DSGE.getShockVariances_(M_));
            % B = B*S;
            % D = D*S;

            % Basic assumption is that D is invertible.
            if size(D, 1) ~= size(D, 2)
                error("Matrix D must be square")
            end
            if rank(D) ~= size(D, 1)
                error("Matrix D must not be singular")
            end
            Phi0 = D*S;

            n = size(C, 1);
            m = size(A, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Case I: C is invertible. In that case, p=q=1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if n==m && rank(C) == n
                p = 1;
                q = 1;

                % Finding AR, MA matrices
                CInv = inv(C);
                DInv = inv(D);
                As = {C*A*CInv};
                Psis = {C*(B - A*CInv*D)*DInv};
                return;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Case II: Follow the general proposition of Morris 2016.
            % Algorithm was slightly adjusted for better implementation.  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Finding smallest M for which F matrix is full column rank. 
            % If F is full column rank for M, then VARMA has AR order M+1. 
            % Thus, our search cannot exceed pMax-1.
            B = B * inv(D);
            F = zeros(0, n);
            p = missing;
            for M = 0:(pMax-1)
                F = [C*A^M; F];
                if rank(F) == m
                    p = M + 1;
                    break;
                end
            end

            if ismissing(p)
                error("Could not find VARMA with p <= %i", pMax)
            end
            q = p;

            FPlus = pinv(F);

            % Recovering AR coefficients
            Bcoeffs = F * A * FPlus;
            Bcoeffs = Bcoeffs(1:n, :);
            As = VAR.coeffsToCellArray_(Bcoeffs);

            % Constructing G matrix to recover MA coefficients
            GFirstRow = eye(n);
            for M = 0:(q-2)
                tmp = C * A^M * B;
                GFirstRow = [GFirstRow tmp];
            end
            G = zeros(n * q, n * q);
            for i = 1:q
                skip = (i - 1) * n;
                r = (i - 1) * n + 1;
                G(r:(r + n - 1), (skip + 1):end) = GFirstRow(:, 1:(end-skip));
            end

            % Recovering MA coefficients
            Psicoeffs = zeros(n, n * q);
            Psicoeffs(:, 1:(end - n)) = G(1:n, (n + 1):end);
            FAFPlusG = F * A * FPlus * G;
            Psicoeffs = Psicoeffs - FAFPlusG(1:n, :);
            FB = F * B;
            Psicoeffs(:, (end - n + 1):end) = Psicoeffs(:, (end - n + 1):end) + FB(1:n, :);
            Psis = VAR.coeffsToCellArray_(Psicoeffs);
        end

        function [A, B, C, D] = getABCD_(M_, oo_, options_)
            % `getABCD_` Obtain the ABCD state-space representation of a DSGE model.
            %
            %   `[A, B, C, D] = getABCD_(M_, oo_, options_)` computes the state-space
            %   representation
            %   $$
            %   \begin{split}
            %   x_t &= Ax_{t-1} + B\varepsilon_t \\
            %   y_t &= Cx_{t-1} + D\varepsilon_t
            %   \end{split}
            %   of a DSGE model estimated using Dynare. Only the minimal state
            %   representation is returned.
            %
            %   ## Arguments
            %   - `M_` (struct): Returned by Dynare.
            %   - `oo_` (struct): Returned by Dynare.
            %   - `options_` (struct): Returned by Dynare.
            %
            %   ## Returns
            %   - `A` (matrix): State transition matrix. See above equation.
            %   - `B` (matrix): Control input matrix capturing exogenous shocks.
            %     See above equation.
            %   - `C` (matrix): Observation matrix mapping state variables to
            %     observed variables. See above equation.
            %   - `D` (matrix): Observation noise matrix. See above equation.
            %
            %   ## Notes
            %   - Requires MATLAB's Control System Toolbox.
            %

            if ~isfield(options_,'varobs_id')
                warning('getABCD: No observables have been defined using the varobs-command.')
                return;
            end

            % Dynare re-orders variables into the order static, backward, mixed, forward.
            % The state variables are the backward and mixed variables.
            % Thus, in the DR (internal order) ordering, the state variables are given by
            % the following indices, where nspred is the number of state variables.
            ipred = M_.nstatic+(1:M_.nspred)';
            % options_.varobs_id is in declaration order. Need to change this to internal DR
            % order for ABCD matrices.
            obs_var=oo_.dr.inv_order_var(options_.varobs_id);

            % get state transition matrices
            [A,B] = kalman_transition_matrix(oo_.dr,ipred,1:M_.nspred);
            % get observation equation matrices
            [C,D] = kalman_transition_matrix(oo_.dr,obs_var,1:M_.nspred);

            % We need the minimum state representation
            if user_has_matlab_license('control_toolbox')
                [A,B,C,D]=minreal(A,B,C,D); % Matlab control toolbox; TODO: find alternative
            else
                error('Control System Toolbox is missing')
            end
        end

        function idxShock = getShockIdx_(shockName, M_)
            shocks = DSGE.dynareCellArrayToVec_(M_.exo_names);
            idxShock = find(shocks == shockName);
        end

        function shockSize = getShockSize_(shockName, M_)
            % `getShockSize_` Obtain the standard deviation of a specified shock.
            %
            %   `shockSize = getShockSize_(shockName, M_)` computes the standard
            %   deviation (size) of a specified shock in a DSGE model
            %   estimated using Dynare.
            %
            %   ## Arguments
            %   - `shockName` (string): The name of the shock whose size and
            %     index are required.
            %   - `M_` (struct): Returned by Dynare.
            %
            %   ## Returns
            %   - `shockSize` (double): The standard deviation of the specified shock.
            %

            idx = DSGE.getShockIdx_(shockName, M_);
            shockSize = sqrt(M_.Sigma_e(idx, idx));
        end

        function SigmaShock = getShockVariances_(M_)
            % `getShockVariances_` Obtain the covariance matrix
            % of structural shocks.
            %
            %   `SigmaShock = getShockVariances_(M_)` returns the
            %   covariance matrix of the structural shocks in a
            %   DSGE model estimated using Dynare. This function
            %   additionally enforces that the shock covariance
            %   matrix is diagonal. Models with correlated shocks
            %   (non-diagonal covariance) are not supported.
            %
            %   ## Arguments
            %   - `M_` (struct): Model structure returned by
            %     Dynare. Must contain the field `Sigma_e`,
            %     representing the covariance matrix of shocks.
            %
            %   ## Returns
            %   - `SigmaShock` (matrix): Covariance matrix of the
            %     model’s structural shocks as given in
            %     `M_.Sigma_e`, provided it is diagonal.
            %

            SigmaShock = M_.Sigma_e;

            % Ensure shock covariance matrix is diagonal
            if ~isequal(SigmaShock, diag(diag(SigmaShock)))
                error(['Non-diagonal shock covariance matrices are not currently supported.']);
            end
        end

        function idx = getVariableIdx_(varname, options_)
            varnames = DSGE.dynareCellArrayToVec_(options_.varobs);
            idx = find(varnames == varname);
        end

        function irfs = varmaIrfs_(Phi0, As, Psis, horizon)
            % `varmaIrfs_` Compute structural impulse response functions (IRFs)
            % for a VARMA model.
            %
            %   `irfs = varmaIrfs_(Phi0, As, Psis, horizon)` computes the
            %   structural impulse response functions (IRFs) of a VARMA model,
            %   given the structural shock impact matrix, autoregressive (AR)
            %   coefficients, and moving average (MA) coefficients.
            %
            %   ## Model Specification
            %   The VARMA model is defined as:
            %   $$
            %   y_t = \sum_{i=1}^{p} A_i y_{t-i} + \sum_{j=1}^{q} \Psi_j u_{t-j} + u_t,
            %   $$
            %   where:
            %   - $u_t = \Phi_0 \varepsilon_t$, with $\varepsilon_t$ being
            %     structural shocks.
            %
            %   ## Arguments
            %   - `Phi0` (matrix): Impact matrix linking structural shocks to
            %     reduced-form errors.
            %   - `As` (cell array): AR coefficient matrices
            %     `{A_1, A_2, ..., A_p}`.
            %   - `Psis` (cell array): MA coefficient matrices
            %     `{Psi_1, Psi_2, ..., Psi_q}`.
            %   - `horizon` (integer): Number of periods for which IRFs are
            %     computed. `horizon=0` means only contemporaneous impulses are
            %     computed which are the same as `Phi0`.
            %
            %   ## Returns
            %   - `irfs` (3D array): Structural IRFs of size `(n, m, horizon+1)`,
            %     where `n` is the number of endogenous variables, `m` is the
            %     number of shocks, assumed to satisfy `m=n`. The IRFs capture
            %     the dynamic response of each variable to a unit shock over
            %     the specified horizon.
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

        function irfs = stateSpaceIrfs_(A, B, C, D, S, horizon)
            % `stateSpaceIrfs_` Compute structural IRFs from a DSGE state-space system.
            %
            %   `irfs = stateSpaceIrfs_(A, B, C, D, S, horizon)` computes
            %   impulse response functions from the state-space representation
            %   returned by `getABCD_`:
            %   $$
            %   \begin{split}
            %   x_t &= Ax_{t-1} + B\varepsilon_t \\
            %   y_t &= Cx_{t-1} + D\varepsilon_t
            %   \end{split}
            %   $$
            %   where the shocks are scaled by `S` so that the responses match
            %   Dynare's convention of one-standard-deviation structural shocks.
            %
            %   ## Arguments
            %   - `A` (matrix): State transition matrix.
            %   - `B` (matrix): Shock loading matrix in the state equation.
            %   - `C` (matrix): Observation matrix.
            %   - `D` (matrix): Contemporaneous shock impact matrix.
            %   - `S` (matrix): Diagonal matrix of shock standard deviations.
            %   - `horizon` (integer): Maximum IRF horizon. `horizon=0`
            %     returns only the contemporaneous impact.
            %
            %   ## Returns
            %   - `irfs` (3D array): Structural IRFs of size
            %     `(n_observables, n_shocks, horizon+1)`.
            %
            %   ## Notes
            %   - Horizon 0 equals `D*S`.
            %   - For `h >= 1`, the response is `C*A^(h-1)*B*S`.
            %
            %   See also `getABCD_`, `varmaIrfs_`

            nObs = size(C, 1);
            nShocks = size(D, 2);
            irfs = zeros(nObs, nShocks, horizon+1);
            irfs(:, :, 1) = D * S;

            stateResponse = B * S;
            for h=1:horizon
                irfs(:, :, h+1) = C * stateResponse;
                stateResponse = A * stateResponse;
            end
        end

    end

    methods

        function obj = DSGE(M_, options_, oo_)
            obj.M_ = M_;
            obj.options_ = options_;
            obj.oo_ = oo_;
        end

        function idx = getVariableIdx(obj, varname)
            % `getVariableIdx` Get index of an observed variable in DSGE model.
            %
            %   `idx = getVariableIdx(obj, varname)` returns the index of
            %   `varname` in the list of observed variables of a DSGE model.
            %
            %   ## Arguments
            %   - `obj` (DSGE): DSGE model object.
            %   - `varname` (char): Name of the observed variable.
            %
            %   ## Returns
            %   - `idx` (integer): Index of the observed variable.
            %
            %   See also `getShockIdx`, `getShockSize`
            idx = DSGE.getVariableIdx_(varname, obj.options_);
        end
        function idx = getShockIdx(obj, shockname)
            % `getShockIdx` Get index of a structural shock in DSGE model.
            %
            %   `idx = getShockIdx(obj, shockname)` returns the index of
            %   `shockname` in the list of shocks of a DSGE model.
            %
            %   ## Arguments
            %   - `obj` (DSGE): DSGE model object.
            %   - `shockname` (char): Name of the structural shock.
            %
            %   ## Returns
            %   - `idx` (integer): Index of the structural shock.
            %
            %   See also `getVariableIdx`, `getShockSize`
            idx = DSGE.getShockIdx_(shockname, obj.M_);
        end
        function shockSize = getShockSize(obj, shockname)
            % `getShockSize` Get size (standard deviation) of a structural shock.
            %
            %   `shockSize = getShockSize(obj, shockname)` returns the standard
            %   deviation of the specified structural shock.
            %
            %   ## Arguments
            %   - `obj` (DSGE): DSGE model object.
            %   - `shockname` (char): Name of the structural shock.
            %
            %   ## Returns
            %   - `shockSize` (number): Standard deviation of the shock.
            %
            %   See also `getVariableIdx`, `getShockIdx`
            shockSize = DSGE.getShockSize_(shockname, obj.M_);
        end

        function [Phi0, As, Psis] = coeffs(obj)
            % `coeffs` Return VARMA coefficients from the DSGE model.
            %
            %   `[Phi0, As, Psis] = coeffs(obj)` computes and returns the
            %   VARMA coefficients of the DSGE model.
            %
            %   ## Returns
            %   - `Phi0` (matrix): Contemporaneous impact effects of
            %     structural shocks.
            %   - `As` (cell array): Reduced-form AR coefficient matrices.
            %   - `Psis` (cell array): Reduced-form MA coefficient matrices.
            %
            %   ## Notes
            %   - Internally calls `dynareToVarma_` to extract VARMA form.
            %   - The method follows Morris (2016) for VARMA approximation.
            %
            %   See also `DSGE.dynareToVarma_`
            [Phi0, As, Psis, p, q] = DSGE.dynareToVarma_(obj.M_, obj.oo_, obj.options_);
        end

        function varnames = getVariableNames(obj)
            varnames = DSGE.dynareCellArrayToVec_(obj.options_.varobs);
        end

        function shocks = getShockNames(obj)
            % `getShockNames` returns a vector of shock names.
            shocks = DSGE.dynareCellArrayToVec_(obj.M_.exo_names);
        end

        function flag = isStructural(obj)
            flag = true;
        end

        function flag = isFitted(obj)
            flag = true;
        end
        function fit(obj)
            error("DSGE must be estimated using Dynare.");
        end
        function fitAndSelect(obj)
            error("DSGE must be estimated using Dynare.");
        end

        function getInputData(obj)
            error("No input data provided.")
        end
        function getIndependent(obj)
            error("No input data provided.")
        end
        function getDependent(obj)
            error("No input data provided.")
        end
        function nobs(obj)
            error("No input data provided.")
        end
        function residuals(obj)
            error("No input data provided.")
        end
        function fitted(obj)
            error("No input data provided.")
        end

        function irfObj = IRF(obj, maxHorizon)
            % `IRF` Compute impulse response functions for DSGE model.
            %
            %   `irfObj = IRF(obj, maxHorizon)` computes IRFs of the DSGE
            %   model up to horizon `maxHorizon`. Uses the VARMA representation
            %   of the DSGE.
            %
            %   ## Arguments
            %   - `obj` (DSGE): DSGE model object.
            %   - `maxHorizon` (integer): Maximum forecast horizon.
            %
            %   ## Returns
            %   - `irfObj` (IRFContainer): Container with computed IRFs.
            %
            %   ## Notes
            %   - Uses VARMA representation for IRF computation.
            %
            %   See also `coeffs`, `dynareToVarma_`, `varmaIrfs_`
            [Phi0, As, Psis] = obj.coeffs();
            irfs = DSGE.varmaIrfs_(Phi0, As, Psis, maxHorizon);
            varnames = obj.getVariableNames();
            irfObj = IRFContainer(irfs, varnames, obj);
        end

        function irfObj = IRFStateSpace(obj, maxHorizon)
            % `IRFStateSpace` Compute impulse response functions from the DSGE state-space form.
            %
            %   `irfObj = IRFStateSpace(obj, maxHorizon)` computes IRFs of the
            %   DSGE model directly from the state-space representation
            %   obtained via `getABCD_`.
            %
            %   ## Arguments
            %   - `obj` (DSGE): DSGE model object.
            %   - `maxHorizon` (integer): Maximum forecast horizon.
            %
            %   ## Returns
            %   - `irfObj` (IRFContainer): Container with computed IRFs.
            %
            %   ## Notes
            %   - Structural shocks are scaled by their standard deviations so
            %     that the returned IRFs match Dynare output.
            %   - The resulting IRFs are equivalent to those obtained from the
            %     VARMA representation used by `IRF`.
            %
            %   See also `IRF`, `getABCD_`, `varmaIrfs_`
            [A, B, C, D] = DSGE.getABCD_(obj.M_, obj.oo_, obj.options_);
            S = sqrt(DSGE.getShockVariances_(obj.M_));
            irfs = DSGE.stateSpaceIrfs_(A, B, C, D, S, maxHorizon);
            varnames = obj.getVariableNames();
            irfObj = IRFContainer(irfs, varnames, obj);
        end

        function effects = transmission(obj, shock, condition, order, maxHorizon)
            % `transmission` Compute transmission effects in a DSGE model.
            %
            %   `effects = transmission(obj, shock, condition, order, maxHorizon)`
            %   computes the transmission effects for a `shock` under a
            %   `condition` based on `order`, up to `maxHorizon`.
            %
            %   ## Arguments
            %   - `obj` (DSGE): DSGE model object.
            %   - `shock` (char or int): Shock name or index.
            %   - `condition` (Q): Transmission condition object.
            %   - `order` (cell array or vector): Variable ordering.
            %   - `maxHorizon` (integer): Maximum forecast horizon.
            %
            %   ## Returns
            %   - `effects` (3D array): Transmission effects over horizons:
            %       - First dimension: Endogenous variables.
            %       - Second dimension: Selected shock.
            %       - Third dimension: Horizon.
            %
            %   ## Notes
            %   - `shock` and `order` can be provided by name or index.
            %
            %   See also `DSGE.through`, `DSGE.notThrough`
            if ~ischar(shock) && ~isnumeric(shock)
                error("The shock must either be given as integer or using the shock's name.");
            end
            if ~isa(condition, 'Q')
                error("The provided transmission condition is not valid.")
            end

            shockIdx = shock;
            if ischar(shock)
                shockNames = obj.getShockNames();
                shockIdx = find(cellfun(@(c) isequal(c, shock), shockNames), 1, 'first');
            end

            orderIdx = obj.vars2idx_(order);
            [Phi0, As, Psis, p, q] = DSGE.dynareToVarma_(obj.M_, obj.oo_, obj.options_);
            [B, Omega] = makeSystemsForm(Phi0, As, Psis, orderIdx, maxHorizon);
            effects = transmission(shockIdx, B, Omega, condition, "BOmega", orderIdx);
        end

    end
end
