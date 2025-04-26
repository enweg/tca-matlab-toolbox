classdef VAR < handle & Model
    % `VAR` Vector Autoregressive (VAR) model in matrix form.
    %
    %   A VAR of lag order `p` is specified as:
    %
    %   ```
    %   y_t = C e_t + B_1 y_{t-1} + ... + B_p y_{t-p} + u_t
    %   ```
    %
    %   where:
    %   - `e_t` is a vector of deterministic components (constant, trends, etc).
    %   - `C, B_i` are conformable coefficient matrices.
    %   - `u_t` is vector white noise.
    %
    %   Compact form:
    %
    %   ```
    %   y_t' = z_t' B_+' + u_t'
    %   ```
    %
    %   with:
    %   - `z_t = [e_t; y_{t-1}; ...; y_{t-p}]`
    %   - `B_+ = [C, B_1, ..., B_p]`
    %
    %   Stacking from `t = p+1` to `T`:
    %
    %   ```
    %   Y = X B_+' + U
    %   ```
    %
    %   ## Properties
    %   - `B` (matrix): Coefficient matrix `[C B_1 ... B_p].
    %   - `SigmaU` (matrix): Covariance matrix of the error term.
    %   - `p` (integer): Lag order of the VAR.
    %   - `trendExponents` (vector): Time trend exponents (e.g., `[0, 1]`
    %     implies constant and linear trend).
    %   - `inputData` (table or matrix): Original data used to estimate
    %     the VAR.
    %   - `Y` (matrix): Left-hand side outcomes `y_t`, size `(T-p) x k`.
    %   - `X` (matrix): Right-hand side regressors `z_t`, size `(T-p) x (k*p + m)`
    %     where `m` is the number of deterministic domponents.
    %   - `U` (matrix): Residuals `u_t`, size `(T-p) x k`.
    %   - `Yhat` (matrix): Fitted values `X * B_+'`, size `(T-p) x k`.
    properties
        B                   % Coefficient matrices [C, B_1, ..., B_p]
        SigmaU              % Covariance of error term
        p                   % Lag-length
        trendExponents      % Deterministic trend exponents for time trend

        inputData           % Data used to fit the model
        Y                   % LHS matrix
        X                   % RHS matrix
        U                   % Residual Matrix
        Yhat                % Fitted values
    end

    methods(Static)
        function val = ic_(SigmaU, nCoeffs, ct)
            % `ic_` Compute a generic information criterion value.
            %
            %   `val = ic_(SigmaU, nCoeffs, ct)` returns the value of an
            %   information criterion based on the log determinant of the
            %   residual covariance matrix and a complexity penalty term.
            %
            %   ## Arguments
            %   - `SigmaU` (matrix): Covariance matrix of the VAR residuals.
            %   - `nCoeffs` (integer): Total number of estimated coefficients.
            %   - `ct` (number): Complexity term adjusting for sample size.
            %
            %   ## Returns
            %   - `val` (number): Computed information criterion value.
            val = logdet(SigmaU) + ct * nCoeffs; 
        end
        function val = aic_(SigmaU, nCoeffs, T)
            % `aic_` Compute Akaike Information Criterion (AIC).
            %
            %   `val = aic_(SigmaU, nCoeffs, T)` returns the AIC value given the
            %   residual covariance matrix, number of coefficients, and sample size.
            %
            %   ## Arguments
            %   - `SigmaU` (matrix): Covariance matrix of the VAR residuals.
            %   - `nCoeffs` (integer): Total number of estimated coefficients.
            %   - `T` (integer): Number of effective observations.
            %
            %   ## Returns
            %   - `val` (number): Computed AIC value.
            ct = 2 / T;
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = hqc_(SigmaU, nCoeffs, T)
            % `hqc_` Compute Hannan-Quinn Information Criterion (HQC).
            %
            %   `val = hqc_(SigmaU, nCoeffs, T)` returns the HQC value given the
            %   residual covariance matrix, number of coefficients, and sample size.
            %
            %   ## Arguments
            %   - `SigmaU` (matrix): Covariance matrix of the VAR residuals.
            %   - `nCoeffs` (integer): Total number of estimated coefficients.
            %   - `T` (integer): Number of effective observations.
            %
            %   ## Returns
            %   - `val` (number): Computed HQC value.
            ct = 2 * log(log(T)) / T;
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = sic_(SigmaU, nCoeffs, T)
            % `sic_` Compute Schwarz Information Criterion (SIC/BIC).
            %
            %   `val = sic_(SigmaU, nCoeffs, T)` returns the SIC value given the
            %   residual covariance matrix, number of coefficients, and sample size.
            %
            %   ## Arguments
            %   - `SigmaU` (matrix): Covariance matrix of the VAR residuals.
            %   - `nCoeffs` (integer): Total number of estimated coefficients.
            %   - `T` (integer): Number of effective observations.
            %
            %   ## Returns
            %   - `val` (number): Computed SIC value.
            ct = log(T) / T;
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = bic_(SigmaU, nCoeffs, T)
            % `bic_` Compute Bayesian Information Criterion (SIC/BIC).
            %
            %   `val = sic_(SigmaU, nCoeffs, T)` returns the BIC value given the
            %   residual covariance matrix, number of coefficients, and sample size.
            %
            %   ## Arguments
            %   - `SigmaU` (matrix): Covariance matrix of the VAR residuals.
            %   - `nCoeffs` (integer): Total number of estimated coefficients.
            %   - `T` (integer): Number of effective observations.
            %
            %   ## Returns
            %   - `val` (number): Computed BIC value.
            %
            %   ## Notes 
            %   - BIC is the same as SIC.
            val = VAR.sic_(SigmaU, nCoeffs, T);
        end

        function BCellArray = coeffsToCellArray_(B)
            % `coeffsToCellArray_` Transform coefficient matrix into cell array
            % of lag matrices.
            %
            %   `BCellArray = coeffsToCellArray_(B)` converts the
            %   coefficient matrix `[B_1 B_2 ... B_p]` into a cell array where
            %   each element corresponds to one lag matrix `B_i`.
            %
            %   ## Arguments
            %   - `B` (matrix): Stacked coefficient matrix excluding
            %     deterministic components (i.e., the matrix does not include
            %     the constant or trend coefficients `C`). Size is `(k, k*p)`
            %     where `k` is the number of variables and `p` is the lag order.
            %
            %   ## Returns
            %   - `BCellArray` (cell array): A 1-by-`p` cell array where each
            %     cell contains the `(k x k)` lag coefficient matrix for one lag.
            %
            %   ## Notes
            %   - Assumes that `B` has already been stripped of coefficients on
            %     deterministic components.
            %
            [k, kp] = size(B);
            p = kp / k;
            BCellArray = cell(1, p);
            for i = 1:p
                BCellArray{i} = B(:, ((i-1)*k + 1):(i*k));
            end
        end

        function Y = simulate(errorsOrT, B, varargin)
            % `simulate` Simulate a VAR process given errors or time periods.
            %
            %   `Y = simulate(errorsOrT, B, varargin)` simulates a VAR model using
            %   either provided error terms or by generating new errors from 
            %   a Normal distribution.
            %
            %   ## Arguments
            %   - `errorsOrT` (matrix or integer): Either a `(k x T)` matrix of
            %     error terms or an integer specifying the number of periods `T`
            %     to simulate.
            %   - `B` (matrix): Coefficient matrix `[C B_1 ... B_p]` where `p`
            %     is the lag order.
            %   - `varargin`: Name-value pairs for optional arguments:
            %     - `trendExponents` (vector): Exponents for deterministic
            %       trends. Default is `[0]` (constant term).
            %     - `initial` (vector): Initial values for lags, size `(p*k, 1)`.
            %       Default is zeros.
            %     - `SigmaU` (matrix): Covariance matrix for error generation if
            %       simulating errors. Default is identity matrix.
            %
            %   ## Returns
            %   - `Y` (matrix): Simulated data matrix, size `(T x k)`.

            % Default values
            opts.trendExponents = [0];
            opts.initial = [];  % will become zeros as default
            opts.SigmaU = [];   % will become identity matrix by default
            % Parsing user options (name-value pairs)
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error("VAR.simulate: " + varargin{i} + " is not a valid simulation option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if isscalar(errorsOrT)
                errors = randn(size(B, 1), errorsOrT);
                if ~isempty(opts.SigmaU)
                    errors = chol(opts.SigmaU)' * errors;
                end
            elseif ismatrix(errorsOrT)
                errors = errorsOrT;
            else
                error("VAR.simulate: First argument must either be a k*T matrix of errors terms or an integer denoting the number of periods to simulate.")
            end

            [k, T] = size(errors);
            m = length(opts.trendExponents);
            kp = size(B, 2) - m;

            if mod(kp, k) ~= 0
                error("VAR.simulate: Dimensions of B are wrong.")
            end

            if isempty(opts.initial)
                opts.initial = zeros(kp, 1);
            end

            Y = zeros(size(errors));
            Zt = ones(m + kp, 1);
            Zt((m+1):end) = opts.initial;
            for t = 1:T
                % trend terms
                Zt(1:m) = arrayfun(@(x) t^x, opts.trendExponents);
                Y(:, t) = B * Zt + errors(:,  t);
                if kp > 0
                    % Rotating in the new data point into Zt
                    Zt((m+1+k):end) = Zt((m+1):(end-k));
                    Zt((m+1):(m+k)) = Y(:, t);
                end
            end

            Y = Y';
        end

        function irfs = IRF_(B, p, maxHorizon)
            % `IRF_` Compute impulse response functions for a VAR model.
            %
            %   `irfs = IRF_(B, p, maxHorizon)` computes impulse response
            %   functions (IRFs) for horizons from 0 to `maxHorizon`, given a
            %   coefficient matrix `B`.
            %
            %   ## Arguments
            %   - `B` (matrix): Stacked coefficient matrix `[B_1 ... B_p]`,
            %     excluding deterministic components.
            %   - `p` (integer): Lag order of the VAR.
            %   - `maxHorizon` (integer): Maximum horizon for IRFs.
            %
            %   ## Returns
            %   - `irfs` (3D array): Impulse responses of size `(k x k x
            %     (maxHorizon+1))`, where:
            %       - First dimension: Endogenous variables (responses).
            %       - Second dimension: Shocks (impulses).
            %       - Third dimension: Horizon.
            %
            %   ## Notes
            %   - `B` must not include coefficients on deterministic components
            %     (such as constants or trends).
            k = size(B, 1);
            irfs = zeros(k, k, maxHorizon + 1);
            irfs(:, :, 1) = eye(k);

            for h = 1:maxHorizon
                for j = 1:min(h, p)
                    Bj = B(:, ((j-1)*k+1):(j*k));
                    irfs(:, :, h+1) = irfs(:, :, h+1) + Bj * irfs(:, :, h+1-j);
                end
            end
        end

        function C = makeCompanionMatrix_(B, p, m)
            % `makeCompanionMatrix_` Form the companion matrix of a VAR model.
            %
            %   `C = makeCompanionMatrix_(B, p, m)` constructs the companion
            %   matrix for a VAR(p) model given its coefficient matrix `B`.
            %
            %   ## Arguments
            %   - `B` (matrix): Coefficient matrix `[C B_1 ... B_p]`, where
            %     `C` are coefficients on deterministic components and `B_i`
            %     are lag matrices.
            %   - `p` (integer): Lag order of the VAR.
            %   - `m` (integer): Number of deterministic components.
            %
            %   ## Returns
            %   - `C` (matrix): Companion matrix of the VAR(p) system.
            %
            %   ## Notes
            %   - The companion matrix has the structure:
            %
            %     $$
            %     \begin{bmatrix}
            %     B_1 & B_2 & \cdots & B_p \\
            %     I   & 0   & \cdots & 0 \\
            %     0   & I   & \cdots & 0 \\
            %     \vdots & \vdots & \ddots & \vdots \\
            %     0   & 0   & \cdots & I
            %     \end{bmatrix}
            %     $$
            B = B(:, (m+1):end);
            k = size(B, 1);
            C = diag(ones(k * p - k, 1), -k);
            C(1:k, :) = B;
        end

        function rho = spectralRadius_(C)
            % `spectralRadius_` Compute the spectral radius of a matrix.
            %
            %   `rho = spectralRadius_(C)` returns the spectral radius of the
            %   companion matrix `C`, defined as the maximum absolute value of
            %   its eigenvalues.
            %
            %   ## Arguments
            %   - `C` (matrix): Companion matrix of the VAR model.
            %
            %   ## Returns
            %   - `rho` (number): Spectral radius of the companion matrix.
            %
            %   See also `makeCompanionMatrix_`
            rho = max(abs(eig(C)));
        end
    end

    methods

        function obj = VAR(data, p, varargin)
            % `VAR` Construct a VAR(p) model.
            %
            %   `obj = VAR(data, p, varargin)` creates a VAR with lag
            %   length `p` based on the provided dataset.
            %
            %   ## Arguments
            %   - `data` (table or matrix): Input dataset for the VAR model.
            %   - `p` (integer): Lag order of the VAR.
            %   - `varargin`: Name-value pairs for optional arguments:
            %     - `trendExponents` (vector): Exponents for deterministic
            %       trends. Defaults to `[0]` (constant term).
            %     - `B` (matrix): Coefficient matrix. Default is empty (must be
            %       estimated).
            %     - `SigmaU` (matrix): Covariance matrix of residuals. Default
            %       is empty (must be estimated).
            %
            %   ## Returns
            %   - `obj` (VAR): A VAR model.
            %
            %   See also `fit`, `simulate`

            % Default values
            opts.trendExponents = [0];
            opts.B = [];
            opts.SigmaU = [];

            % Parsing options (name-value pairs)
            for i=1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error("VAR: " + varargin{i} + " is not a valid option.")
                end

                opts.(varargin{i}) = varargin{i+1};
            end

            if istable(data)
                dataMatrix = table2array(data);
            else
                dataMatrix = data;
                varnames = arrayfun(@(i) "Y" + i, 1:size(data, 2));
                data = array2table(data, 'VariableNames', varnames);
            end

            U = [];
            Yhat = [];
            Y = dataMatrix((p+1):end, :);
            X = makeLagMatrix(dataMatrix, p);
            X = X((p+1):end, :);
            for i = 1:length(opts.trendExponents)
                te = opts.trendExponents(i);
                trend = arrayfun(@(x) x^te, (p+1):size(dataMatrix, 1))';
                X = [trend X];
            end

            obj.B = opts.B;
            obj.SigmaU = opts.SigmaU;
            obj.p = p;
            obj.trendExponents = opts.trendExponents;
            obj.inputData = data;
            obj.Y = Y;
            obj.X = X;
            obj.U = U;
            obj.Yhat = Yhat;
        end

        function flag = isFitted(obj)
            if size(obj.Yhat) > 0
                flag = true;
                return;
            end
            flag = false;
        end

        function B = coeffs(obj, excludeDeterministic)
            % `coeffs` Return the VAR coefficient matrix.
            %
            %   `B = coeffs(obj, excludeDeterministic)` returns the VAR
            %   coefficient matrix `[C, B_1, ..., B_p]`. If
            %   `excludeDeterministic` is true, returns `[B_1, ..., B_p]`
            %   instead, excluding deterministic components.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %   - `excludeDeterministic` (logical, optional): If true,
            %     exclude coefficients on deterministic components.
            %     Defaults to false.
            %
            %   ## Returns
            %   - `B` (matrix): VAR coefficient matrix.
            requireFitted(obj);
            if nargin < 2
                excludeDeterministic = false;
            end
            if ~excludeDeterministic
                B = obj.B;
                return;
            end

            m = length(obj.trendExponents);
            B = obj.B(:, (m+1):end);
        end
        
        function Yhat = fitted(obj)
            % `fitted` Return the fitted values of the VAR model.
            %
            %   `Yhat = fitted(obj)` returns the matrix of fitted values
            %   with size `(T-p) x k`, where `T` is the number of observations
            %   and `k` is the number of variables.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `Yhat` (matrix): Matrix of fitted values.
            requireFitted(obj);
            Yhat = obj.Yhat;
        end

        function U = residuals(obj)
            % `residuals` Return the residuals of the VAR model.
            %
            %   `U = residuals(obj)` returns the matrix of VAR residuals
            %   with size `(T-p) x k`, where `T` is the number of observations
            %   and `k` is the number of variables.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `U` (matrix): Matrix of residuals.
            requireFitted(obj);
            U = obj.U;
        end

        function n = nobs(obj)
            n = size(obj.Y, 1);
        end

        function n = ncoeffs(obj)
            requireFitted(obj);
            n = prod(size(coeffs(obj)));
        end

        function Y = getDependent(obj)
            Y = obj.Y;
        end

        function X = getIndependent(obj)
            X = obj.X;
        end

        function data = getInputData(obj)
            data = obj.inputData;
        end

        function varnames = getVariableNames(obj)
            data = getInputData(obj);
            varnames = data.Properties.VariableNames;
        end

        function flag = isStructural(obj)
            flag = false;
        end

        function C = makeCompanionMatrix(obj)
            % `makeCompanionMatrix` Form the companion matrix of the VAR model.
            %
            %   `C = makeCompanionMatrix(obj)` constructs the companion matrix
            %   for the fitted VAR(p) model.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `C` (matrix): Companion matrix of the VAR(p)
            requireFitted(obj);
            m = length(obj.trendExponents);
            C = VAR.makeCompanionMatrix_(coeffs(obj), obj.p, m); 
        end

        function rho = spectralRadius(obj)
            % `spectralRadius` Compute the spectral radius of the VAR model.
            %
            %   `rho = spectralRadius(obj)` returns the spectral radius of the
            %   companion matrix associated with the fitted VAR model.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `rho` (number): Spectral radius of the companion matrix.
            requireFitted(obj);
            rho = VAR.spectralRadius_(makeCompanionMatrix(obj));
        end

        function flag = isStable(obj)
            requireFitted(obj);
            flag = (spectralRadius(obj) < 1);
        end
        function val = aic(obj)
            % `aic` Compute Akaike Information Criterion (AIC) for VAR model.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `val` (number): AIC value.
            %
            %   See also `hqc`, `sic`, `bic`, `fit`
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.aic_(SigmaU, nCoeffs, T);
        end
        function val = hqc(obj)
            % `hqc` Compute Hannan-Quinn Criterion (HQC) for VAR model.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `val` (number): HQC value.
            %
            %   See also `aic`, `sic`, `bic`, `fit`
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.hqc_(SigmaU, nCoeffs, T);
        end
        function val = sic(obj)
            % `sic` Compute Schwarz Information Criterion (SIC) for VAR model.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `val` (number): SIC value.
            %
            %   See also `aic`, `hqc`, `bic`, `fit`
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.sic_(SigmaU, nCoeffs, T);
        end
        function val = bic(obj)
            % `bic` Compute Bayesian Information Criterion (BIC) for VAR model.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %
            %   ## Returns
            %   - `val` (number): BIC value.
            %
            %   See also `aic`, `hqc`, `sic`, `fit`
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.bic_(SigmaU, nCoeffs, T);
        end

        function fit(obj)
            % `fit` Estimate the VAR model using ordinary least squares (OLS).
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            X = obj.X;
            Y = obj.Y;

            obj.B = (Y' * X) / (X' * X);
            obj.Yhat = X * obj.B';
            obj.U = Y - obj.Yhat;
            obj.SigmaU = obj.U' * obj.U / nobs(obj);
        end

        function [modelBest, icTable] = fitAndSelect(obj, icFunction)
            % `fitAndSelect` Estimate and select the best VAR model by IC.
            %
            %   `[modelBest, icTable] = fitAndSelect(obj, icFunction)` fits
            %   the VAR model for different lag lengths and selects the one
            %   minimizing the information criterion. Maximum lag length is 
            %   given by the lag length of the provided model. 
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %   - `icFunction` (function handle, optional): Information
            %     criterion function to minimize. Defaults to `aic`.
            %
            %   ## Returns
            %   - `modelBest` (VAR): Best fitting model.
            %   - `icTable` (table): Table of lag lengths and IC values.
            %
            %   See also `fit`, `aic`, `bic`, `hqc`, `sic`
            if nargin < 2
                icFunction = @VAR.aic_;
            end

            pMax = obj.p;
            ps = 0:pMax;
            ics = zeros(length(ps), 1);

            obj.fit();
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            T = nobs(obj);
            icBest = icFunction(SigmaU, nCoeffs, T);
            ics(pMax+1) = icBest;
            modelBest = obj;

            for p = (pMax-1):-1:0
                modelTmp = VAR(getInputData(obj), p, 'trendExponents', obj.trendExponents);
                modelTmp.fit();
                nCoeffsTmp = ncoeffs(modelTmp);
                UTmp = residuals(modelTmp);
                UTmp = UTmp((pMax-p+1):end, :);
                SigmaUTmp = UTmp' * UTmp / T;

                icTmp = icFunction(SigmaUTmp, nCoeffsTmp, T);
                ics(p+1) = icTmp;
                if icTmp < icBest
                    icBest = icTmp;
                    modelBest = modelTmp;
                end
            end

            icTable = table(ps', ics, 'VariableNames', {'p', 'IC'});
        end

        function irfObj = IRF(obj, maxHorizon, varargin)
            % `IRF` Compute impulse response functions for the VAR model.
            %
            %   `irfObj = IRF(obj, maxHorizon, varargin)` computes IRFs up to
            %   horizon `maxHorizon`. If an `identificationMethod` is provided,
            %   structural IRFs are computed.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %   - `maxHorizon` (integer): Maximum horizon for IRFs.
            %   - `varargin`: Name-value pairs for options:
            %     - `identificationMethod` (an `IdentificationMethod`): Optional 
            %       method to compute structural IRFs.
            %
            %   ## Returns
            %   - `irfObj` (IRFContainer): Object containing computed IRFs.
            %
            %   ## Notes
            %   - Without an identification method, reduced-form IRFs are
            %     computed.
            %   - With an identification method, structural IRFs are computed.
            %
            %   See also `IRF_`, `IRFContainer`, `fit`, `IdentificationMethod`
            opts.identificationMethod = missing;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error("VAR.IRF: " + varargin{i} + " is not a valid option.")
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            requireFitted(obj);
            if ismissing(opts.identificationMethod)
                B = obj.coeffs(true);
                irfs = VAR.IRF_(B, obj.p, maxHorizon);
            else
                % This way users can easily implement new identification methods.
                irfs = opts.identificationMethod.identifyIrfs(obj, maxHorizon);
            end
            varnames = getVariableNames(obj);
            irfObj = IRFContainer(irfs, varnames, obj, opts.identificationMethod);
        end

        function effects = transmission(obj, shock, condition, order, maxHorizon, varargin)
            % `transmission` Compute transmission effects in a VAR model.
            %
            %   `effects = transmission(obj, shock, condition, order, maxHorizon, varargin)`
            %   computes the transmission effects of a `shock` under a
            %   `condition`, using the transmission matrix defined by `order`,
            %   up to `maxHorizon`.
            %
            %   ## Arguments
            %   - `obj` (VAR): VAR model object.
            %   - `shock` (integer): Index of the shock variable.
            %   - `condition` (Q): Transmission condition object.
            %   - `order` (cell array of char): Variable transmission ordering.
            %   - `maxHorizon` (integer): Maximum horizon.
            %   - `varargin`: Name-value pairs for options:
            %     - `identificationMethod` (`IdentificationMethod`): Required 
            %       method to compute structural IRFs.
            %
            %   ## Returns
            %   - `effects` (3D array): Transmission effects over horizons, where:
            %       - First dimension: Endogenous variables (responses).
            %       - Second dimension: Shocks (of size one for the single
            %         selected shock).
            %       - Third dimension: Horizon.
            %
            %   See also `VAR.through`, `VAR.notThrough`

            opts.identificationMethod = missing;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end
            if ismissing(opts.identificationMethod)
                error("To compute transmission effects from a VAR an identificationMethod must be provided.")
            end

            effects = opts.identificationMethod.identifyTransmission(obj, shock, condition, order, maxHorizon);
        end
    end
end
