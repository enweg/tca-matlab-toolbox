classdef ExternalInstrument < IdentificationMethod
    properties
        treatment           % number of variale name
        instruments         % numbers or variable names (if names must be cell)
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
            opts.normalisingHorizon = 0;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if ischar(instruments)
                instruments = {instruments};
            end

            obj.treatment = treatment; 
            obj.instruments = instruments; 
            obj.normalisingHorizon = opts.normalisingHorizon;
        end

        function varargout = identify(obj, model)
            data = model.getInputData();
            idxTreatment = findVariableIndex(data, model.treatment);
            idxInstrumentTreatment = findVariableIndex(data, obj.treatment);
            if idxTreatment ~= idxInstrumentTreatment
                error("Treatment of LP and ExternalInstrument differ.");
            end
            if iscell(obj.instruments)
                % instruments were provided by name
                idxInstruments = cellfun(@(x) findVariableIndex(data, x), obj.instruments);
            else
                % instruments are already provided by index
                idxInstruments = obj.instruments;
            end
            if ~all(idxInstruments < idxTreatment)
                error("Instruments must come before treatment in data.");
            end

            m = model.includeConstant;
            % model.X = [constant contemporaneous treatment lag]
            Z = model.X(:, idxInstruments + m);
            % we assume instruments are completely excluded from model
            X = model.X(:, setdiff(1:size(model.X, 2), idxInstruments + m));
            % we also need to adjus the treatment index
            model.treatment = idxTreatment - length(idxInstruments);


            if obj.normalisingHorizon > 0
                % lead the treatment column in X to adjust for which horizon 
                % the unit effect normalisation applies to
                nlead = obj.normalisingHorizon;
                X(:, model.treatment+m) = makeLeadMatrix(X(:, model.treatment+m), nlead);
                % remove NaNs at end of data
                X = X(1:(end-nlead), :);
                model.Y = model.Y(1:(end-nlead), :, :);
                Z = Z(1:(end-nlead), :);
            end
            model.X = X;

            % Adding all other exogenous variables to Z
            % these are all variables that are not the treatment variable
            Z = [Z X(:, setdiff(1:size(X, 2), model.treatment + m))];

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

            varargout{1} = coeffs;
        end

    end
end
