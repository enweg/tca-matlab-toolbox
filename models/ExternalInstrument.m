classdef ExternalInstrument < IdentificationMethod
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
