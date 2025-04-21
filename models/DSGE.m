classdef DSGE < handle & Model
    properties
        % The following are all structs returned by Dynare
        M_
        options_
        oo_
    end

    methods (Static)

        function checkDynare_()
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
            % `defineOrder` Determine the ordering (the transmission matrix) of
            % observed variables in a DSGE model. 
            %
            %   `defineOrder(vars, options_)` returns an ordering vector `order` that
            %   maps the variables in `vars` to their corresponding positions in the
            %   list of observed variables of a DSGE model estimated using Dynare. It 
            %   therefore defines the transmission matrix and can be used in `makeB`, 
            %   `makeOmega`, `makeSystemsForm`, `makeConditionY`, `notThrough`, `through`, 
            %   and `transmission`. 
            %
            %   ## Arguments
            %   - `vars` (vector): A list of observed variable names.
            %   - `options_` (struct): Returned by Dynare.
            %
            %   ## Returns
            %   - `order` (vector): A vector containing the indices of `vars` in the
            %     original observed variable list.
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

        function [Phi0, As, Psis, p, q] = dynareToVarma_(M_, oo_, options_, maxKappa)
            % `dynareToVarma` Transform a DSGE model into a VARMA representation.
            %
            %   `[Phi0, As, Psis, p, q] = dynareToVarma(M_, oo_, options_, maxKappa)`
            %   converts a linearized DSGE model estimated using Dynare into a VARMA form,
            %   following the method of Morris (2016).
            %
            %   ## Arguments
            %   - `M_` (struct): Returned by Dynare. 
            %   - `oo_` (struct): Returned by Dynare. 
            %   - `options_` (struct): Returned by Dynare. 
            %   - `maxKappa` (integer, optional): Tuning parameter in Morris (2016). 
            %      Related to the maximum AR order via `maxArOrder=maxKappa+1`.
            %      Defaults to 20. 
            %
            %   ## Returns
            %   - `Phi0` (matrix): Impact matrix linking structural shocks to reduced-form errors.
            %   - `As` (cell array): AR coefficient matrices `{A_1, A_2, ..., A_p}`.
            %   - `Psis` (cell array): MA coefficient matrices `{Psi_1, Psi_2, ..., Psi_q}`.
            %   - `p` (integer): The determined autoregressive order of the VARMA representation.
            %   - `q` (integer): The determined moving average order of the VARMA representation.
            %
            %   ## Methodology
            %   The function follows the approach outlined in Morris (2016) and returns
            %   a VARMA of the form: 
            %   $$
            %   y_t = \sum_{i=1}^{p} A_i y_{t-i} + \sum_{j=1}^{q} \Psi_j u_{t-j} + u_t,
            %   $$
            %   where:
            %   - $u_t = \Phi_0 \varepsilon_t$, with $\varepsilon_t$ being structural shocks.
            %
            %   ## Reference
            %   - Morris, S. D. (2016). "VARMA representation of DSGE models." *Economics Letters*, 138, 30–33.
            %     [https://doi.org/10.1016/j.econlet.2015.11.027](https://doi.org/10.1016/j.econlet.2015.11.027)
            %
            %   See also `getABCD`, `varmaIrfs`.

            % Checking if Dynare is setup
            DSGE.checkDynare_();

            if ~isfield(options_, "varobs")
                error("dynareToVarma: No observed variables were defined in the mod file.")
            end

            % Default choice for maximum VAR order following notation in Morris 2016. 
            if nargin==3
                maxKappa = 20;
            end

            [A, B, C, D] = DSGE.getABCD_(M_, oo_, options_);

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
                Phi0 = D;
                As = {C*A*CInv};
                Psis = {C*(B - A*CInv*D)*DInv};
                return;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Case II: Follow the general proposition of Morris 2016
            % but adjusted for our notation.
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

            As = cell(1, kappa+1);
            Phi0 = ThetaPlus * FPlus * G{1};
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
                Psis{k} = ThetaPlus * (FPlus * G{k+1} - A * FPlus * G{k}) * A0;
            end
            Psis{end} = ThetaPlus * (B - A * FPlus * G{kappa+1}) * A0;
        end

        function [A, B, C, D] = getABCD_(M_, oo_, options_)
            % `getABCD` Obtain the ABCD state-space representation of a DSGE model.
            %
            %   `[A, B, C, D] = getABCD(M_, oo_, options_)` computes the state-space
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
            %   - `B` (matrix): Control input matrix capturing exogenous shocks. See above equation.
            %   - `C` (matrix): Observation matrix mapping state variables to observed variables. See above equation.
            %   - `D` (matrix): Observation noise matrix. See above equation.
            %
            %   ## Notes
            %   - Requires MATLAB's Control Toolbox. 
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
            % `getShockSize` Obtain the standard deviation of a specified shock.
            %
            %   `shockSize = getShockSize(shockName, M_)` computes the standard
            %   deviation (size) of a specified shock in a DSGE model estimated using Dynare.
            %
            %   ## Arguments
            %   - `shockName` (string): The name of the shock whose size and index are required.
            %   - `M_` (struct): Returned by Dynare.
            %
            %   ## Returns
            %   - `shockSize` (double): The standard deviation of the specified shock.
            %

            idx = DSGE.getShockIdx_(shockName, M_);
            shockSize = sqrt(M_.Sigma_e(idx, idx));
        end

        function idx = getVariableIdx_(varname, options_)
            varnames = DSGE.dynareCellArrayToVec_(options_.varobs);
            idx = find(varnames == varname);
        end

        function irfs = varmaIrfs_(Phi0, As, Psis, horizon)
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

    end

    methods

        function obj = DSGE(M_, options_, oo_)
            obj.M_ = M_; 
            obj.options_ = options_; 
            obj.oo_ = oo_;
        end

        function idx = getVariableIdx(obj, varname)
            idx = DSGE.getVariableIdx_(varname, obj.options_);
        end
        function idx = getShockIdx(obj, shockname)
            idx = DSGE.getShockIdx_(shockname, obj.M_);
        end
        function shockSize = getShockSize(obj, shockname)
            shockSize = DSGE.getShockSize_(shockname, obj.M_);
        end

        function [Phi0, As, Psis] = coeffs(obj)
            % returns the VARMA coeffs
            [Phi0, As, Psis, p, q] = DSGE.dynareToVarma_(obj.M_, obj.oo_, obj.options_);
        end

        function varnames = getVariableNames(obj)
            varnames = DSGE.dynareCellArrayToVec_(obj.options_.varobs);
        end

        function shocks = getShockNames(obj)
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
            [Phi0, As, Psis] = obj.coeffs();
            irfs = DSGE.varmaIrfs_(Phi0, As, Psis, maxHorizon);
            varnames = obj.getVariableNames();
            irfObj = IRFContainer(irfs, varnames, obj);
        end

        % TODO: test
        function effects = transmission(obj, shock, condition, order, maxHorizon)
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
            shockSize = obj.getShockSize(shock);
            [Phi0, As, Psis, p, q] = DSGE.dynareToVarma_(obj.M_, obj.oo_, obj.options_);
            [B, Omega] = makeSystemsForm(Phi0, As, Psis, orderIdx, maxHorizon);
            effects = transmission(shockIdx, B, Omega, condition, "BOmega", orderIdx) * shockSize;
        end

    end
end
