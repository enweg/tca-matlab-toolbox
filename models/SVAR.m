classdef SVAR < handle & Model 
    % `SVAR` Structural Vector Autoregressive (SVAR) model.
    %
    %   An SVAR of lag order `p` is specified as:
    %
    %   ```
    %   A0 y_t = C e_t + A1 y_{t-1} + ... + Ap y_{t-p} + \varepsilon_t
    %   ```
    %
    %   where:
    %   - `e_t` is a vector of deterministic components (constant, trends).
    %   - `C`, `A0`, `Ai` are conformable matrices.
    %   - `\varepsilon_t` are structural shocks.
    %
    %   Assuming `A0` is invertible, the model can be rewritten as a
    %   reduced-form VAR:
    %
    %   ```
    %   y_t'A0' = z_t' A_+' + u_t'
    %   ```
    %
    %   where:
    %   - `z_t = [e_t; y_{t-1}; ...; y_{t-p}]`
    %   - `A_+ = [C, A_1, ..., A_p]`
    %
    %   Assuming `A0` is invertible, the reduced-form VAR can be obtained as 
    %
    %   ```
    %   y_t' = z_t' A_+'(A_0')^{-1} + u_t'(A_0')^{-1}
    %   ```
    %   which can be represented using a `VAR` object.
    %
    %   ## Properties
    %   - `A0` (matrix): Contemporaneous coefficient matrix.
    %   - `APlus` (matrix): Stacked coefficient matrix `[C A1 ... Ap]`.
    %   - `p` (integer): Lag order of the (S)VAR model.
    %   - `trendExponents` (vector): Time trend exponents (e.g., `[0, 1]`).
    %   - `VARModel` (VAR): Reduced-form VAR representation.
    %
    %   See also `VAR`
    properties
        A0
        APlus
        p
        trendExponents
        VARModel
    end

    methods (Static)
        function Y = simulate(shocksOrT, A0, APlus, varargin)
            % `simulate` Simulate a Structural VAR (SVAR) process.
            %
            %   `Y = simulate(shocksOrT, A0, APlus, varargin)` simulates an SVAR
            %   model using either provided structural shocks or by generating
            %   random shocks.
            %
            %   ## Arguments
            %   - `shocksOrT` (matrix or integer): Either a `(k x T)` matrix of
            %     structural shocks or an integer specifying the number of
            %     periods to simulate.
            %   - `A0` (matrix): Contemporaneous coefficient matrix.
            %   - `APlus` (matrix): Coefficient matrix `[C A1 ... Ap]`.
            %   - `varargin`: Name-value pairs for options:
            %     - `trendExponents` (vector): Exponents for deterministic
            %       trends. Defaults to `[0]` (constant).
            %     - `initial` (vector): Initial lag values, default is zeros.
            %
            %   ## Returns
            %   - `Y` (matrix): Simulated data matrix, size `(T x k)`.
            %
            %   ## Notes
            %   - If `shocksOrT` is a scalar, shocks are drawn from a standard
            %     normal distribution with identity covariance.
            %
            %   See also `VAR.simulate`
            opts.trendExponents = [0];
            opts.initial = [];
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if isscalar(shocksOrT)
                shocks = randn(size(A0, 1), shocksOrT);
            elseif ismatrix(shocksOrT)
                shocks = shocksOrT;
            else
                error("First argument must either be a k*T matrix of shocks or an integer denoting the number of periods to simulate.")
            end

            Phi0 = inv(A0);
            errors = Phi0 * shocks;
            B = Phi0 * APlus;
            
            Y = VAR.simulate(errors, B, 'trendExponents', opts.trendExponents, 'initial', opts.initial);
        end

        function irfs = IRF_(A0, APlus, p, maxHorizon)
            % `IRF_` Compute impulse response functions for an SVAR model.
            %
            %   `irfs = IRF_(A0, APlus, p, maxHorizon)` computes the structural
            %   impulse response functions up to `maxHorizon` given the
            %   contemporaneous matrix `A0` and stacked lag matrices `APlus`.
            %
            %   ## Arguments
            %   - `A0` (matrix): Contemporaneous coefficient matrix.
            %   - `APlus` (matrix): Stacked lag coefficient matrix `[A1 ... Ap]`,
            %     excluding deterministic components.
            %   - `p` (integer): Lag order of the SVAR.
            %   - `maxHorizon` (integer): Maximum horizon for IRFs.
            %
            %   ## Returns
            %   - `irfs` (3D array): Impulse response functions of size
            %     `(k x k x (maxHorizon+1))`, where:
            %       - First dimension: Endogenous variables (responses).
            %       - Second dimension: Structural shocks.
            %       - Third dimension: Horizon.
            %
            %   ## Notes
            %   - `APlus` must not include coefficients on deterministic terms.
            %
            %   See also `VAR.IRF_`
            Phi0 = inv(A0);
            B = Phi0 * APlus; 

            irfs = VAR.IRF_(B, p, maxHorizon);
            for h = 0:maxHorizon
                irfs(:, :, h+1) = irfs(:, :, h+1) * Phi0;
            end
        end
    end

    methods
        function obj = SVAR(data, p, varargin)
            % `SVAR` Construct a Structural VAR (SVAR) model object.
            %
            %   `obj = SVAR(data, p, varargin)` creates an SVAR object with lag
            %   length `p`, based on the provided dataset and structural
            %   specification.
            %
            %   ## Arguments
            %   - `data` (table or matrix): Input dataset for the SVAR model.
            %   - `p` (integer): Lag order of the SVAR.
            %   - `varargin`: Name-value pairs for optional arguments:
            %     - `trendExponents` (vector): Exponents for deterministic
            %       trends. Defaults to `[0]` (constant term).
            %     - `A0` (matrix): Contemporaneous coefficient matrix.
            %     - `APlus` (matrix): Stacked coefficient matrix `[A1, ..., Ap]`.
            %     - `VARModel` (VAR): Precomputed reduced-form VAR model.
            %
            %   ## Returns
            %   - `obj` (SVAR): Constructed SVAR model object.
            %
            %   ## Notes
            %   - If `VARModel`, `A0`, `APlus` are not provided, they are 
            %     estimated from the data.
            %
            %   See also `VAR`

            opts.trendExponents = [0];
            opts.A0 = [];
            opts.APlus = [];
            opts.VARModel = missing;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if ismissing(opts.VARModel)
                opts.VARModel = VAR(data, p, 'trendExponents', opts.trendExponents);
            end

            obj.A0 = opts.A0; 
            obj.APlus = opts.APlus;
            obj.p = p;
            obj.trendExponents = opts.trendExponents;
            obj.VARModel = opts.VARModel;
        end

        function flag = isFitted(obj)
            flag = obj.VARModel.isFitted();
        end

        function [A0, APlus] = coeffs(obj, excludeDeterministic)
            % `coeffs` Return the SVAR coefficient matrices.
            %
            %   `[A0, APlus] = coeffs(obj, excludeDeterministic)` returns the
            %   contemporaneous matrix `A0` and the lag coefficient matrix
            %   `APlus` for the SVAR model.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %   - `excludeDeterministic` (logical, optional): If true,
            %     exclude coefficients on deterministic components from `APlus`.
            %     Defaults to false.
            %
            %   ## Returns
            %   - `A0` (matrix): Contemporaneous coefficient matrix.
            %   - `APlus` (matrix): Coefficient matrix. If
            %     `excludeDeterministic` is true, returns only lag matrices
            %     `[A1 ... Ap]`. Otherwise `APlus = [C A_1 ... A_p]` is returned.

            requireFitted(obj);
            if nargin < 2
                excludeDeterministic = false;
            end

            A0 = obj.A0; 
            if ~excludeDeterministic
                APlus = obj.APlus;
                return;
            end

            m = length(obj.trendExponents);
            APlus = obj.APlus(:, (m+1):end);
        end

        function Yhat = fitted(obj)
            % `fitted` Return the fitted values of the SVAR model.
            %
            %   `Yhat = fitted(obj)` returns the fitted values from the
            %   reduced-form VAR model associated with the SVAR.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `Yhat` (matrix): Matrix of fitted values `(T-p) x k`.
            requireFitted(obj);
            Yhat = obj.VARModel.fitted();
        end

        function U = residuals(obj)
            % `residuals` Return the residuals of the SVAR model.
            %
            %   `U = residuals(obj)` returns the residuals from the
            %   reduced-form VAR model associated with the SVAR.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `U` (matrix): Matrix of residuals `(T-p) x k`.
            requireFitted(obj);
            U = obj.VARModel.residuals();
        end

        function n = nobs(obj)
            n = obj.VARModel.nobs();
        end

        function Y = getDependent(obj)
            Y = obj.VARModel.getDependent();
        end

        function X = getIndependent(obj)
            X = obj.VARModel.getIndependent();
        end

        function data = getInputData(obj)
            data = obj.VARModel.getInputData();
        end

        function varnames = getVariableNames(obj)
            varnames = obj.VARModel.getVariableNames();
        end

        function flag = isStructural(obj)
            flag = true;
        end

        function C = makeCompanionMatrix(obj)
            % `makeCompanionMatrix` Form the companion matrix of the SVAR model.
            %
            %   `C = makeCompanionMatrix(obj)` returns the companion matrix
            %   associated with the reduced-form VAR model of the SVAR.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `C` (matrix): Companion matrix.
            requireFitted(obj);
            C = obj.VARModel.makeCompanionMatrix();
        end

        function rho = spectralRadius(obj)
            % `spectralRadius` Compute the spectral radius of the SVAR model.
            %
            %   `rho = spectralRadius(obj)` returns the spectral radius of the
            %   companion matrix associated with the SVAR.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `rho` (number): Spectral radius of the companion matrix.
            requireFitted(obj);
            rho = obj.VARModel.spectralRadius();
        end

        function flag = isStable(obj)
            % `isStable` Check if the SVAR model is stable.
            %
            %   `flag = isStable(obj)` returns true if the spectral radius of
            %   the companion matrix is less than 1.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `flag` (logical): True if the model is stable, false otherwise.
            requireFitted(obj);
            flag = obj.VARModel.isStable();
        end

        function val = aic(obj)
            % `aic` Compute Akaike Information Criterion (AIC) for SVAR model.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `val` (number): AIC value.
            %
            %   See also `bic`, `hqc`, `sic`, `VAR.aic`
            requireFitted(obj);
            val = obj.VARModel.aic();
        end
        function val = bic(obj)
            % `bic` Compute Bayesian Information Criterion (BIC) for SVAR model.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `val` (number): BIC value.
            %
            %   See also `aic`, `hqc`, `sic`, `VAR.bic`
            requireFitted(obj);
            val = obj.VARModel.bic();
        end
        function val = sic(obj)
            % `sic` Compute Schwarz Information Criterion (SIC) for SVAR model.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `val` (number): SIC value.
            %
            %   See also `aic`, `hqc`, `bic`, `VAR.sic`
            requireFitted(obj);
            val = obj.VARModel.sic();
        end
        function val = hqc(obj)
            % `hqc` Compute Hannan-Quinn Criterion (HQC) for SVAR model.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %
            %   ## Returns
            %   - `val` (number): HQC value.
            %
            %   See also `aic`, `bic`, `sic`, `VAR.hqc`
            requireFitted(obj);
            val = obj.VARModel.hqc();
        end

        function fit(obj, identificationMethod)
            % `fit` Estimate the SVAR model using an identification method.
            %
            %   `fit(obj, identificationMethod)` first fits the reduced-form
            %   VAR model using ordinary least squares (OLS), then identifies
            %   the structural matrices using the provided identification method.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %   - `identificationMethod` (object): An object of type
            %     `IdentificationMethod` used to identify `A0` and `APlus`.
            %
            %   ## Notes
            %   - `identificationMethod` must implement an `identify` method
            %     taking a VAR object and returning `A0` and `APlus`.
            %
            %   See also `VAR.fit`, `IdentificationMethod`
            obj.VARModel.fit();
            [A0, APlus] = identificationMethod.identify(obj.VARModel);
            obj.A0 = A0;
            obj.APlus = APlus;
        end

        function [modelBest, icTable] = fitAndSelect(obj, identificationMethod, icFunction)
            % `fitAndSelect` Estimate and select the best SVAR model by IC.
            %
            %   `[modelBest, icTable] = fitAndSelect(obj, identificationMethod, icFunction)`
            %   fits the SVAR model for different lag lengths and selects the one
            %   minimizing the information criterion. The maximuml lag length 
            %   is determined by the lag length of the given model. 
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %   - `identificationMethod` (`IdentificationMethod`): Identification 
            %     method for SVAR.
            %   - `icFunction` (function handle, optional): Information criterion
            %     function to minimize. Defaults to `aic`.
            %
            %   ## Returns
            %   - `modelBest` (SVAR): Best fitting SVAR model.
            %   - `icTable` (table): Table of lag lengths and IC values.
            %
            %   ## Notes
            %   - Maximum lag length is given by the lag length of the provided model.
            %
            %   See also `SVAr.fit`, `aic`, `bic`, `hqc`, `sic`, `VAR.fitAndSelect`
            if nargin < 3
                icFunction = @VAR.aic_;
            end

            [modelVARBest, icTable] = fitAndSelect(obj.VARModel, icFunction);
            [A0, APlus] = identificationMethod.identify(modelVARBest);

            data = obj.getInputData();
            modelBest = SVAR(data, modelVARBest.p, 'trendExponents', obj.trendExponents, 'A0', A0, 'APlus', APlus, 'VARModel', modelVARBest);
        end

        function irfObj = IRF(obj, maxHorizon)
            % `IRF` Compute structural impulse response functions for SVAR.
            %
            %   `irfObj = IRF(obj, maxHorizon)` computes structural IRFs up to
            %   horizon `maxHorizon` from an estimated SVAR model.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %   - `maxHorizon` (integer): Maximum horizon for IRFs.
            %
            %   ## Returns
            %   - `irfObj` (IRFContainer): Container with computed IRFs.
            %
            %   ## Notes
            %   - The IRFs have dimensions `(k x k x (maxHorizon+1))`:
            %       - First dimension: Responding variables.
            %       - Second dimension: Structural shocks.
            %       - Third dimension: Horizons.
            %
            %   See also `VAR.IRF`, ~IRFContainer`
            requireFitted(obj);
            B = obj.VARModel.coeffs(true);
            irfs = VAR.IRF_(B, obj.p, maxHorizon);

            Phi0 = inv(obj.A0);
            for h = 0:maxHorizon
                irfs(:, :, h+1) = irfs(:, :, h+1)  * Phi0;
            end

            varnames = getVariableNames(obj);
            irfObj = IRFContainer(irfs, varnames, obj);
        end

        function effects = transmission(obj, shock, condition, order, maxHorizon, varargin)
            % `transmission` Compute transmission effects in an SVAR model.
            %
            %   `effects = transmission(obj, shock, condition, order, maxHorizon, varargin)`
            %   computes the transmission effects of a `shock` under a `condition`,
            %   using the transmission matrix defined by `order`, up to `maxHorizon`.
            %
            %   ## Arguments
            %   - `obj` (SVAR): SVAR model object.
            %   - `shock` (integer): Index of the shock variable.
            %   - `condition` (Q): Transmission condition object.
            %   - `order` (cell array of char): Variable transmission ordering.
            %   - `maxHorizon` (integer): Maximum horizon.
            %   - `varargin`: Name-value pairs for options:
            %     - `identificationMethod` (IdentificationMethod, optional):
            %       Method to identify the SVAR if not yet fitted.
            %
            %   ## Returns
            %   - `effects` (3D array): Transmission effects over horizons, where:
            %       - First dimension: Endogenous variables (responses).
            %       - Second dimension: Shocks (only the selected shock).
            %       - Third dimension: Horizon.
            %
            %   ## Notes
            %   - If `identificationMethod` is provided, the model is fitted first.
            %
            %   See also `SVAR.through`, `SVAR.notThrough`, `VAR.transmission`
            opts.identificationMethod = missing;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end
            if ~ismissing(opts.identificationMethod)
                obj.fit(opts.identificationMethod);
            end

            requireFitted(obj);
            if ~isnumeric(shock)
                error("Shock must be provided as integer for SVAR models.")
            end
            if ~isa(condition, 'Q')
                error("The provided transmission condition is not valid.")
            end

            shockIdx = shock;
            orderIdx = obj.vars2idx_(order);

            B = obj.VARModel.coeffs(true);
            Bs = VAR.coeffsToCellArray_(B);
            Phi0 = inv(obj.A0);
            Psis = cell(0, 0);
            [B, Omega] = makeSystemsForm(Phi0, Bs, Psis, orderIdx, maxHorizon);

            effects = transmission(shockIdx, B, Omega, condition, "BOmega", orderIdx); 
        end

    end
end
