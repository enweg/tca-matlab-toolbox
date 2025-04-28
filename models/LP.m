classdef LP < handle & Model
    % `LP` Local Projection (LP) model for estimating IRFs.
    %
    %   Local Projection (LP) model for estimating impulse response functions (IRFs)
    %   in a flexible and semi-parametric manner.
    %
    %   Each LP regression estimates the dynamic response of an outcome variable at
    %   future horizon `h` to a one-period change in a treatment variable at time `t`,
    %   controlling for contemporaneous and lagged covariates.
    %
    %   The regression model is specified as:
    %
    %   ```
    %   w_{i,t+h} = \mu_{i,h} + \theta_{i,h} x_t + \gamma_{i,h}' r_t +
    %               \sum_{l=1}^p \delta_{i,h,l} w_{t-l} + \xi_{i,h,t}
    %   ```
    %
    %   where `w_t = (r_t', x_t, q_t')` and:
    %   - `x_t` is the treatment variable
    %   - `r_t` contains contemporaneous controls (all variables before `x_t`)
    %   - `p` is the number of lags included
    %   - `\theta_{i,h}` is the relative IRF of `x_t` on the `i`-th variable at 
    %     horizon `h`.
    %
    %   The treatment variable may be endogenous. Structural interpretation of IRFs
    %   can be achieved using valid instruments—see `ExternalInstrument` for one such
    %   method. If the treatment satisfies a conditional ignorability assumption 
    %   (a recursive assumption in macro), then the coefficient has a structural 
    %   interpretation even without the use of instruments. For this to hold, 
    %   `x_t - E(x_t|r_t, w_{t-1}, ..., w_{t-p})` must be equal to the structural shock.
    %
    %   ## Properties
    %   - `data` (table or matrix): Input time series dataset.
    %   - `treatment` (char or integer): Treatment variable.
    %   - `p` (integer): Number of lags.
    %   - `horizons` (vector): Forecast horizons for projections.
    %   - `includeConstant` (logical): Whether to include an intercept.
    %   - `B` (array): Coefficient estimates per horizon.
    %   - `Y` (array): Dependent variables per horizon.
    %   - `X` (matrix): Common regressor matrix.
    %   - `U` (array): Residuals per horizon.
    %   - `Yhat` (array): Fitted values per horizon.
    %
    properties
        data 
        treatment
        p
        horizons
        includeConstant

        B
        Y
        X
        U
        Yhat
    end

    methods (Static)
        function [X, Y] = createXY_(data, treatment, p, horizons, varargin)
            % `createXY_` Construct regressor matrix `X` and response array `Y`.
            %
            %   `[X, Y] = createXY_(data, treatment, p, horizons, varargin)`
            %   prepares the design matrices for local projection estimation.
            %
            %   ## Arguments
            %   - `data` (matrix or table): Input time series dataset.
            %   - `treatment` (char or int): Treatment variable.
            %   - `p` (int): Lag length.
            %   - `horizons` (vector): Forecast horizons.
            %   - `varargin`: Name-value pairs for options:
            %     - `includeConstant` (logical): Include constant column in `X`
            %       (Default is true).
            %
            %   ## Returns
            %   - `X` (matrix): Common regressor matrix for all horizons.
            %   - `Y` (3D array): Outcome variables stacked over horizons (along
            %     the third dimension).
            %
            %   ## Notes
            %   - `X` is structured as [deterministic contemporaneous treatment lagged]
            %   - `Y` stacks future outcomes across selected horizons along the
            %     third dimension.
            %   - Shape of `Y`: (observations, variables, horizons).
            %
            %   ## Details
            %   - All variables ordered before the treatment in `data` are 
            %     included as contemporaneous controls. 
            %

            opts.includeConstant = true;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if ~ismatrix(data)
                error("Data should be provided in the form of a matrix or table.");
            end

            dataMatrix = data;
            if istable(data)
                dataMatrix = table2array(data);
            end

            [T, k] = size(dataMatrix);
            idxTreatment = findVariableIndex(data, treatment);

            dataMatrixLag = makeLagMatrix(dataMatrix, p);
            X = [dataMatrix((p+1):end, 1:idxTreatment) dataMatrixLag((p+1):end, :)];
            if opts.includeConstant
                X = [ones(T-p, 1) X];
            end

            maxHorizon = max(horizons);
            Y = makeLeadMatrix(dataMatrix, maxHorizon);
            Y = [dataMatrix Y];
            Y = Y((p+1):end, :);
            Y = reshape(Y, T-p, k, []);
            Y = Y(:, :, horizons + 1);
        end
    end

    methods
        function obj = LP(data, treatment, p, horizons, varargin)
            % `LP` Construct a Local Projection (LP) model object.
            %
            %   `obj = LP(data, treatment, p, horizons, varargin)` initializes
            %   an LP object for estimating impulse response functions.
            %
            %   ## Arguments
            %   - `data` (matrix or table): Input time series dataset.
            %   - `treatment` (char or int): Treatment variable.
            %   - `p` (integer): Lag length.
            %   - `horizons` (vector): Forecast horizons.
            %   - `varargin`: Name-value pairs for options:
            %     - `includeConstant` (logical): Include constant in regressors
            %       (Defaults to true).

            opts.includeConstant = true; 
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            [X, Y] = LP.createXY_(data, treatment, p, horizons, 'includeConstant', opts.includeConstant);

            if ~istable(data)
                data = array2table(data, 'VariableNames', arrayfun(@(i) "Y" + i, 1:size(data, 2)));
            end
            obj.data = data; 
            obj.treatment = treatment; 
            obj.p = p;
            obj.horizons = horizons; 
            obj.includeConstant = opts.includeConstant;
            obj.B = [];
            obj.Y = Y;
            obj.X = X;
            obj.U = [];
            obj.Yhat = [];
        end

        function flag = isFitted(obj)
            flag = size(obj.Yhat, 1) > 0;
        end

        function B = coeffs(obj, excludeDeterministic)
            % `coeffs` Return coefficient estimates from LP model.
            %
            %   `B = coeffs(obj, excludeDeterministic)` returns the estimated
            %   coefficients from the fitted LP model.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %   - `excludeDeterministic` (logical, optional): If true, excludes
            %     constant terms from the coefficients. Defaults to false.
            %
            %   ## Returns
            %   - `B` (3D array): Coefficients array with dimensions:
            %       - First dimension: Outcome variable.
            %       - Second dimension: Regressors.
            %       - Third dimension: Horizon.
            %
            %   See also `fit`
            requireFitted(obj);
            if nargin < 2
                excludeDeterministic = false;
            end

            if ~excludeDeterministic
                B = obj.B;
                return;
            end

            B = obj.B(:, (obj.includeConstant+1):end, :);
        end

        function Yhat = fitted(obj)
            % `fitted` Return the fitted values from the LP model.
            %
            %   `Yhat = fitted(obj)` returns the fitted values obtained from
            %   the local projection regressions.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %
            %   ## Returns
            %   - `Yhat` (3D array): Fitted values with dimensions:
            %       - First dimension: Time.
            %       - Second dimension: Outcome variable.
            %       - Third dimension: Horizon.
            %
            %   See also `residuals`, `coeffs`, `fit`
            requireFitted(obj);
            Yhat = obj.Yhat;
        end

        function U = residuals(obj)
            % `residuals` Return residuals from the LP model.
            %
            %   `U = residuals(obj)` returns the residuals from the local
            %   projection regressions.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %
            %   ## Returns
            %   - `U` (3D array): Residuals with dimensions:
            %       - First dimension: Time.
            %       - Second dimension: Outcome variable.
            %       - Third dimension: Horizon.
            %
            %   See also `fitted`, `coeffs`, `fit`
            requireFitted(obj);
            U = obj.U;
        end

        function n = nobs(obj)
            n = size(obj.data, 1) - obj.p - obj.horizons;
        end

        function Y = getDependent(obj)
            Y = obj.Y;
        end

        function X = getIndependent(obj)
            X = obj.X;
        end

        function data = getInputData(obj)
            data = obj.data;
        end

        function varnames = getVariableNames(obj)
            data = getInputData(obj);
            varnames = data.Properties.VariableNames;
        end

        function flag = isStructural(obj)
            flag = true;
        end

        function fit(obj, identificationMethod)
            % `fit` Estimate the LP model with an identification method.
            %
            %   `fit(obj, identificationMethod)` estimates the LP model,
            %   identifying causal effects with respect to the treatment.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %   - `identificationMethod` (object, optional): Identification
            %     method. Must be of type `IdentificationMethod`. Defaults to
            %     `Recursive`.
            %
            %   See also `coeffs`, `fitted`, `residuals`, `Recursive`
            if nargin < 2
                identificationMethod = Recursive();
            end
            B = identificationMethod.identify(obj);
            
            Yhat = nan(size(obj.Y));
            U = nan(size(obj.Y));

            for i = 1:length(obj.horizons)
                Yhat(:, :, i) = obj.X * B(:, :, i)'; 
                U(:, :, i) = obj.Y(:, :, i) - Yhat(:, :, i);
            end

            obj.B = B;
            obj.Yhat = Yhat; 
            obj.U = U;
        end

        function [modelBest, icTable] = fitAndSelect(obj, identificationMethod, icFunction)
            % `fitAndSelect` Select optimal lag length for LP model.
            %
            %   `[modelBest, icTable] = fitAndSelect(obj, identificationMethod, icFunction)`
            %   selects the optimal lag length based on an equivalent VAR model.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %   - `identificationMethod` (`IdentificationMethod`, optional): 
            %     Identification method. Defaults to `Recursive`.
            %   - `icFunction` (function handle, optional): Information criterion
            %     function to minimize. Defaults to `aic`.
            %
            %   ## Returns
            %   - `modelBest` (LP): Best fitting LP model.
            %   - `icTable` (table): Table of lag lengths and IC values.
            %
            %   ## Notes
            %   - Maximum lag length considered is the lag length of `obj`.
            %
            %   See also `fit`, `VAR.fitAndSelect`, `Recursive`
            if nargin == 1
                identificationMethod = Recursive();
            end
            if nargin < 3
                icFunction = @VAR.aic_;
            end

            % best model will be selected by finding the best VAR alternative
            data = obj.getInputData();
            if obj.includeConstant
                trendExponents = [0];
            else
                trendExponents = [];
            end
            modelVAR = VAR(data, obj.p, 'trendExponents', trendExponents);
            [modelVARBest, icTable] = modelVAR.fitAndSelect(icFunction);

            modelBest = LP(data, obj.treatment, modelVARBest.p, obj.horizons, 'includeConstant', obj.includeConstant);
            modelBest.fit(identificationMethod);
        end

        function irfObj = IRF(obj, maxHorizon, varargin)
            % `IRF` Compute impulse response functions from LP model.
            %
            %   `irfObj = IRF(obj, maxHorizon, varargin)` computes IRFs up to
            %   `maxHorizon` based on the LP model.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %   - `maxHorizon` (integer): Maximum forecast horizon.
            %   - `varargin`: Name-value pairs for options:
            %     - `identificationMethod` (optional): Identification method.
            %
            %   ## Returns
            %   - `irfObj` (IRFContainer): Container with computed IRFs.
            %
            %   ## Notes
            %   - If `identificationMethod` is provided, LP is refitted first.
            %   - The IRFs have dimensions `(k x k x (maxHorizon+1))`:
            %       - First dimension: Responding variables.
            %       - Second dimension: Shocks.
            %       - Third dimension: Horizon.
            %
            %   See also `fit`, `IRFContainer`

            opts.identificationMethod = missing; 
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if ~isequal(obj.horizons, 0:maxHorizon)
                error("IRF horizons differ from LP horizons.");
            end

            if ~ismissing(opts.identificationMethod)
                obj.fit(opts.identificationMethod);
                irfObj = obj.IRF(maxHorizon);
                return;
            end

            requireFitted(obj);
            data = obj.getInputData();
            k = size(obj.Y, 2);
            idxTreatment = findVariableIndex(data, obj.treatment);

            irfs = nan(k, k, maxHorizon + 1);
            B = obj.coeffs(true);
            irfs(:, idxTreatment, :) = B(:, idxTreatment, :);
            varnames = data.Properties.VariableNames;
            irfObj = IRFContainer(irfs, varnames, obj, opts.identificationMethod);
        end

        function effects = transmission(obj, shock, condition, order, maxHorizon, varargin)
            % `transmission` Compute transmission effects in an LP model.
            %
            %   `effects = transmission(obj, shock, condition, order, maxHorizon, varargin)`
            %   computes transmission effects for a `shock` satisfying a
            %   `condition`, based on the ordering `order`, up to `maxHorizon`.
            %
            %   ## Arguments
            %   - `obj` (LP): LP model object.
            %   - `shock` (integer): Index of the shock variable.
            %   - `condition` (Q): Transmission condition object.
            %   - `order` (cell array of char): Variable transmission ordering.
            %   - `maxHorizon` (integer): Maximum horizon.
            %   - `varargin`: Name-value pairs for options:
            %     - `identificationMethod` (optional): Identification method.
            %
            %   ## Returns
            %   - `effects` (3D array): Transmission effects over horizons:
            %       - First dimension: Endogenous variables (responses).
            %       - Second dimension: Selected shock.
            %       - Third dimension: Horizon.
            %
            %   ## Notes
            %   - If `identificationMethod` is provided, the LP model is refitted.
            %
            %   See also `LP.through`, `LP.notThrough`, `LP.IRF`

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
                error("Shock must be provided as integer for LP models.");
            end
            if ~isa(condition, 'Q')
                error("The provided transmission condition is not valid.");
            end

            shockIdx = shock;
            orderIdx = obj.vars2idx_(order);
            irfsStructural = obj.IRF(maxHorizon).irfs;
            irfsStructural = irfsStructural(orderIdx, :, :);

            data = obj.getInputData();
            k = size(irfsStructural, 1);
            irfsOrthogonal = nan(k, k, maxHorizon + 1);
            for treatment = 1:k
                modelTmp = LP(data(:, orderIdx), treatment, obj.p, 0:maxHorizon, 'includeConstant', obj.includeConstant);
                modelTmp.fit(Recursive());
                irfsTmp = modelTmp.IRF(maxHorizon).irfs;
                irfsOrthogonal(:, treatment, :) = irfsTmp(:, treatment, :);
            end

            irfsStructural = toTransmissionIrfs(irfsStructural);
            irfsOrthogonal = toTransmissionIrfs(irfsOrthogonal);
            effects = transmission(shockIdx, irfsStructural, irfsOrthogonal, condition, "irf", orderIdx); 
        end
    end
end
