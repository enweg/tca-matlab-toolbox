classdef InternalInstrument < IdentificationMethod
    % `InternalInstrument` Identify IRFs using internal instruments.
    %
    %   This class identifies impulse responses in (S)VARs using the
    %   internal instruments method proposed by Plagborg-Moller and
    %   Wolf (2021).
    %
    %   ## Notes
    %   - IRFs are computed from Cholesky-orthogonalized shocks with
    %     respect to the selected `instrument`.
    %   - Shocks are normalized by the response of the
    %     `normalisingVariable` at `normalisingHorizon`.
    %   - The `instrument` can be specified as an integer (index) or
    %     character (variable name).
    %   - Defaults:
    %     - `instrument = 1` (first variable)
    %     - `normalisingHorizon = 0`
    %
    %   ## Reference
    %   - Plagborg-Moller, M., & Wolf, C. K. (2021).
    %     "Local Projections and VARs Estimate the Same Impulse Responses."
    %     *Econometrica*, 89(2), 955-980.
    %     [https://doi.org/10.3982/ecta17813]
    %
    %   See also `IdentificationMethod`, `Recursive`

    properties
        instrument          % Instrument variable (name or index)
        normalisingVariable % Variable whose response normalises the shock (name or index)
        normalisingHorizon  % Horizon at which the normalising response is one
    end

    methods
        function obj = InternalInstrument(normalisingVariable, varargin)
            % `InternalInstrument` Construct an internal instrument identification method.
            %
            %   `obj = InternalInstrument(normalisingVariable)` identifies the shock to
            %   the first variable, normalised such that the response of
            %   `normalisingVariable` at horizon `0` is one.
            %   `obj = InternalInstrument(normalisingVariable, 'instrument', instrument, 'normalisingHorizon', h)`
            %   sets the instrument and the normalising horizon.
            %
            %   ## Arguments
            %   - `normalisingVariable` (char or integer): Variable whose response is
            %     normalised to one.
            %   - `varargin`: Name-value pairs:
            %     - `'instrument'` (char or integer): Instrument variable. Defaults to `1`.
            %     - `'normalisingHorizon'` (integer): Horizon of the normalisation.
            %       Defaults to `0`.
            %
            %   ## Returns
            %   - `obj` (InternalInstrument): Identification method.
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

        function effects = identifyVARTransmission_(obj, B, SigmaU, p, data, trendExponents, shock, condition, order, maxHorizon)

            irfsStructural = obj.identifyVARIrfs_(B, SigmaU, p, maxHorizon);
            irfsStructural = irfsStructural(order, :, :);
            % Getting orthogonal IRFs using a temporary model. This is 
            % not the most efficient, but the most robust. 
            modelTmp = SVAR(data(:, order), p, 'trendExponents', trendExponents);
            modelTmp.fit(Recursive());
            irfsOrthogonal = modelTmp.IRF(maxHorizon).irfs();

            irfsStructural = toTransmissionIrfs(irfsStructural);
            irfsOrthogonal = toTransmissionIrfs(irfsOrthogonal);

            effects = transmission(shock, irfsStructural, irfsOrthogonal, condition, "irf", order);
        end

        function irfs = identifyIrfs(obj, model, maxHorizon)
            % `identifyIrfs` Identify structural IRFs from a reduced-form model.
            %
            %   `irfs = identifyIrfs(obj, model, maxHorizon)` returns IRFs up to
            %   horizon `maxHorizon` to the shock identified by the internal instrument. IRFs to all other shocks are `NaN`.
            %
            %   ## Arguments
            %   - `obj` (InternalInstrument): Identification method.
            %   - `model` (Model): Estimated reduced-form model. Currently only `VAR`
            %     models are supported.
            %   - `maxHorizon` (integer): Maximum IRF horizon.
            %
            %   ## Returns
            %   - `irfs` (3D array): IRFs with dimensions (response variable, shock,
            %     horizon).
            %
            %   See also `IdentificationMethod`, `VAR.IRF`

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
            % `identify` Not supported for internal instruments.
            %
            %   Always throws an error, since internal instruments only identify IRFs
            %   and not the entire structural model. Use `identifyIrfs` instead.
            switch class(model)
                case 'VAR'
                    error("InternalInstruments can only be used to identify IRFs of SVARs but not the entire SVAR.");
                otherwise
                    error("InternalInstrument identification of " + class(model) + " is not implemented.");
            end
        end

        function effects = identifyTransmission(obj, model, shock, condition, order, maxHorizon)
            % `identifyTransmission` Compute transmission effects using this identification method.
            %
            %   `effects = identifyTransmission(obj, model, shock, condition, order, maxHorizon)`
            %   computes the effect of `shock` through the channel defined by
            %   `condition` under the ordering `order`, up to horizon `maxHorizon`.
            %
            %   ## Arguments
            %   - `obj` (InternalInstrument): Identification method.
            %   - `model` (Model): Estimated reduced-form model. Currently only `VAR`
            %     models are supported.
            %   - `shock` (integer): Index of the shock.
            %   - `condition` (Q): Transmission condition.
            %   - `order` (cell array or vector): Ordering of the variables given by
            %     the transmission matrix, as names or indices.
            %   - `maxHorizon` (integer): Maximum horizon.
            %
            %   ## Returns
            %   - `effects` (3D array): Transmission effects with dimensions (response
            %     variable, shock, horizon).
            %
            %   See also `VAR.transmission`, `transmission`
            if ~isnumeric(shock)
                error("Shock must be provided as integer.");
            end
            if ~isa(condition, 'Q')
                error("Invalid transmission condition.");
            end
            shockIdx = shock; 
            orderIdx = model.vars2idx_(order);

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
                    data = model.getInputData();
                    trendExponents = model.trendExponents;
                    effects = obj.identifyVARTransmission_(B, SigmaU, p, data, trendExponents, shockIdx, condition, orderIdx, maxHorizon);
                otherwise
                    error(class(model) + " not supported.")
            end
        end
    end
end
