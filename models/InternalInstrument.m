classdef InternalInstrument < IdentificationMethod
    properties
        instrument
        normalisingVariable
        normalisingHorizon
    end

    methods
        function obj = InternalInstrument(normalisingVariable, varargin)
            opts.instrument = 1;
            opts.normalisingHorizon = 0;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error(varargin{i} + " is not a valid option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            obj.instrument = opts.instrument;
            obj.normalisingVariable = normalisingVariable; 
            obj.normalisingHorizon = opts.normalisingHorizon;
        end

        function irfs = identifyVARIrfs_(obj, B, SigmaU, p, maxHorizon)
            if isa(obj.instrument, 'string') || isa(obj.instrument, 'char')
                error("Instrument must be defined as an integer.")
            end
            if isa(obj.normalisingVariable, 'string') || isa(obj.normalisingVariable, 'char')
                error("Normalising variable must be defined as an integer.")
            end
            k = size(B, 1);

            irfsCholesky = Recursive.identifyVARIrfs_(B, SigmaU, p, maxHorizon);
            normaliseBy = irfsCholesky(obj.normalisingVariable, obj.instrument, obj.normalisingHorizon + 1);
            irfsCholesky = irfsCholesky ./ normaliseBy;
            irfsCholesky(:, setdiff(1:k, obj.instrument), :) = NaN;
            irfs = irfsCholesky;
        end

        function irfs = identifyIrfs(obj, model, maxHorizon)

            if ~isnumeric(obj.instrument)
                varnames = model.getVariableNames();
                idxInstrument = find(cellfun(@(c) isequal(c, obj.instrument), varnames), 1, 'first');  
                obj = InternalInstrument(obj.normalisingVariable, 'instrument', idxInstrument, 'normalisingHorizon', obj.normalisingHorizon);
            end
            if ~isnumeric(obj.normalisingVariable)
                varnames = model.getVariableNames();
                idxNormalisingVariable = find(cellfun(@(c) isequal(c, obj.normalisingVariable), varnames), 1, 'first');
                obj = InternalInstrument(idxNormalisingVariable, 'instrument', obj.instrument, 'normalisingHorizon', obj.normalisingHorizon);
            end

            switch class(model)
                case 'VAR'
                    B = model.coeffs(true);
                    SigmaU = model.SigmaU;
                    p = model.p;
                    irfs = obj.identifyVARIrfs_(B, SigmaU, p, maxHorizon);
                otherwise
                    error("InternalInstrument identification of IRFs has not been implemented for model " + class(model));
            end
        end

        function varargout = identify(obj, model)
            switch class(model)
                case 'VAR'
                    error("InternalInstruments can only be used to identify IRFs of SVARs but not the entire SVAR.");
                otherwise
                    error("InternalInstrument identification of " + class(model) + " is not implemented.");
            end
        end
    end
end
