classdef ExternalInstrument < IdentificationMethod
    % `ExternalInstrument` Identify IRFs using external instruments.
    %
    %   This class identifies structural IRFs from (S)VARs using the
    %   external instrument approach proposed by Stock and Watson (2018).
    %
    %   ## Notes
    %   - IRFs are computed relative to the reaction of the
    %     `treatment` variable at `normalisingHorizon`.
    %   - `instruments` must be a matrix or table spanning the same
    %     period as the model data.
    %   - Defaults:
    %     - `normalisingHorizon = 0`
    %
    %   ## Reference
    %   - Stock, J. H., & Watson, M. W. (2018).
    %     "Identification and Estimation of Dynamic Causal Effects in
    %     Macroeconomics Using External Instruments."
    %     *The Economic Journal*, 128(610), 917-948.
    %     [https://doi.org/10.1111/ecoj.12593]
    %
    %   See also `IdentificationMethod`, `Recursive`, `InternalInstrument`

    properties
        treatment           % number of variale name
        instruments         % matrix or table
        normalisingHorizon  % number
    end

    methods (Static)
        function B = fit2SLS_(X, Y, Z)
            Xhat = Z * ((Z' * Z) \ (Z' * X));
            B = (Xhat' * Xhat) \ (Xhat' * Y);
        end
    end

    methods
        function obj = ExternalInstrument(treatment, instruments, varargin)
            % `ExternalInstrument` Construct an external instrument identification method.
            %
            %   `obj = ExternalInstrument(treatment, instruments)` instruments the
            %   `treatment` variable with `instruments`.
            %   `obj = ExternalInstrument(treatment, instruments, 'normalisingHorizon', h)`
            %   also sets the horizon at which the treatment's response is normalised.
            %
            %   ## Arguments
            %   - `treatment` (char or integer): Treatment variable, by name or index.
            %     Must match the treatment of the `LP`.
            %   - `instruments` (matrix or table): Instruments, spanning the same
            %     period as the model data.
            %   - `varargin`: Name-value pairs:
            %     - `'normalisingHorizon'` (integer): Defaults to `0`.
            %
            %   ## Returns
            %   - `obj` (ExternalInstrument): Identification method.
            opts.normalisingHorizon = 0;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if ~ismatrix(instruments)
                error("Instruments must be provided as matrix or table.");
            end

            obj.treatment = treatment; 
            obj.instruments = instruments; 
            obj.normalisingHorizon = opts.normalisingHorizon;
        end

        function varargout = identify(obj, model)
            % `identify` Estimate an LP using external instruments.
            %
            %   `coeffs = identify(obj, model)` estimates the coefficients of the `LP`
            %   `model` per horizon by two-stage least squares, instrumenting the
            %   treatment with `obj.instruments`.
            %
            %   ## Arguments
            %   - `obj` (ExternalInstrument): Identification method.
            %   - `model` (LP): Local projection model. Other models are not supported.
            %
            %   ## Returns
            %   - `coeffs` (3D array): Coefficient estimates per horizon.
            %
            %   See also `LP.fit`
            switch class(model)
                case 'LP'
                    varargout{1} = obj.identifyLP_(model);
                otherwise
                    error(class(model) + " not supported.");
            end
        end

        function coeffs = identifyLP_(obj, model)
            data = model.getInputData();
            idxTreatment = findVariableIndex(data, model.treatment);
            idxInstrumentTreatment = findVariableIndex(data, obj.treatment);
            if idxTreatment ~= idxInstrumentTreatment
                error("Treatment of LP and ExternalInstrument differ.");
            end
            if size(obj.instruments, 1) ~= size(data, 1)
                error("Instruments must be observed over the same period as the data.");
            end

            Z = obj.instruments;
            if istable(Z)
                Z = table2array(Z);
            end
            Z = Z((model.p+1):end, :);
            
            m = model.includeConstant;
            if obj.normalisingHorizon > 0
                % lead the treatment column in X to adjust for which horizon 
                % the unit effect normalisation applies to
                nlead = obj.normalisingHorizon;
                model.X(:, idxTreatment+m) = makeLeadMatrix(model.X(:, idxTreatment+m), nlead);
                % remove NaNs at end of data
                model.X = model.X(1:(end-nlead), :);
                model.Y = model.Y(1:(end-nlead), :, :);
                Z = Z(1:(end-nlead), :);
            end

            % Adding all other exogenous variables to Z
            % these are all variables that are not the treatment variable
            Z = [Z model.X(:, setdiff(1:size(model.X, 2), idxTreatment + m))];

            k = size(model.Y, 2);
            numCoeffs = size(model.X, 2);
            coeffs = nan(k, numCoeffs, length(model.horizons));
            for i = 1:length(model.horizons)
                h = model.horizons(i);
                X = model.X(1:(end-h), :);
                Y = model.Y(1:(end-h), :, i);
                Zi = Z(1:(end-h), :);
                coeffs(:, :, i) = ExternalInstrument.fit2SLS_(X, Y, Zi)';
            end
        end
    end
end
