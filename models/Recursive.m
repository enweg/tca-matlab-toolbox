classdef Recursive < IdentificationMethod 
    % `Recursive` Identify reduced-form models using recursivity.
    %
    %   This class implements identification under a recursive
    %   (Cholesky) structure, assuming conditional ignorability.
    %
    %   ## Notes
    %   - For SVARs, this corresponds to Cholesky identification.
    %   - For LPs, it assumes conditioning on prior ordered variables.
    %   - Commonly used in macroeconomics when contemporaneous shocks
    %     are assumed to have a lower triangular impact.
    %
    %   See also `IdentificationMethod`, `SVAR`, `LP`
    methods (Static)

        function [A0, APlus] = identifyVAR_(B, SigmaU)
            L = chol(SigmaU)';
            A0 = inv(L);
            APlus = A0 * B;
        end

        function irfs = identifyVARIrfs_(B, SigmaU, p, maxHorizon)
            [A0, Aplus] = Recursive.identifyVAR_(B, SigmaU);
            Phi0 = inv(A0);
            irfs = VAR.IRF_(B, p, maxHorizon);
            for h = 0:maxHorizon
                irfs(:, :, h+1) = irfs(:, :, h+1) * Phi0;
            end
        end

        function effects = identifyVARTransmission_(B, SigmaU, shock, condition, order, maxHorizon)
            [A0, Aplus] = Recursive.identifyVAR_(B, SigmaU);
            Phi0 = inv(A0);
            Bs = VAR.coeffsToCellArray_(B);
            Psis = cell(0, 0);
            [B, Omega] = makeSystemsForm(Phi0, Bs, Psis, order, maxHorizon);
            effects = transmission(shock, B, Omega, condition, "BOmega", order);
        end

        function coeffs = identifyLP_(X, Y, horizons)
            k = size(Y, 2);
            numCoeffs = size(X, 2);
            numHorizons = length(horizons);

            coeffs = nan(k, numCoeffs, numHorizons);
            for i = 1:numHorizons
                h = horizons(i);
                XTmp = X(1:(end-h), :);
                YTmp = Y(1:(end-h), :, i);

                coeffs(:, :, i) = (YTmp' * XTmp) / (XTmp' * XTmp);
            end
        end

    end

    methods
        
        function irfs = identifyIrfs(obj, model, maxHorizon)
            % `identifyIrfs` Identify structural IRFs from a reduced-form model.
            %
            %   `irfs = identifyIrfs(obj, model, maxHorizon)` returns IRFs up to
            %   horizon `maxHorizon` using a Cholesky decomposition of the reduced-form error covariance.
            %
            %   ## Arguments
            %   - `obj` (Recursive): Identification method.
            %   - `model` (Model): Estimated reduced-form model. Currently only `VAR`
            %     models are supported.
            %   - `maxHorizon` (integer): Maximum IRF horizon.
            %
            %   ## Returns
            %   - `irfs` (3D array): IRFs with dimensions (response variable, shock,
            %     horizon).
            %
            %   See also `IdentificationMethod`, `VAR.IRF`
            switch class(model)
                case 'VAR'
                    B = model.coeffs(true);
                    SigmaU = model.SigmaU;
                    p = model.p;
                    irfs = Recursive.identifyVARIrfs_(B, SigmaU, p, maxHorizon);
                otherwise
                    error("Recursive identification of IRFs has not been implemented for model " + class(model));
            end
        end

        function varargout = identify(obj, model)
            % `identify` Identify the structural model using a recursive scheme.
            %
            %   `[A0, APlus] = identify(obj, model)` identifies the structural matrices
            %   of an SVAR from the estimated `VAR` `model` using a Cholesky
            %   decomposition. `coeffs = identify(obj, model)` estimates the
            %   coefficients of an `LP` `model` per horizon, treating the treatment
            %   variable as recursively identified.
            %
            %   ## Arguments
            %   - `obj` (Recursive): Identification method.
            %   - `model` (VAR or LP): Model to identify.
            %
            %   ## Returns
            %   - For a `VAR`: `A0` (contemporaneous matrix) and `APlus` (stacked
            %     coefficient matrix `[C A1 ... Ap]`).
            %   - For an `LP`: the coefficient estimates per horizon.
            %
            %   See also `IdentificationMethod`, `SVAR.fit`, `LP.fit`
            switch class(model)
                case 'VAR'
                    B = model.coeffs();
                    SigmaU = model.SigmaU;
                    [A0, APlus] = Recursive.identifyVAR_(B, SigmaU);
                    varargout{1} = A0;
                    varargout{2} = APlus;
                case 'LP'
                    X = model.X;
                    Y = model.Y; 
                    horizons = model.horizons; 
                    varargout{1} = Recursive.identifyLP_(X, Y, horizons);
                otherwise
                    error("Recursive identification of " + class(model) + " is not implemented.");
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
            %   - `obj` (Recursive): Identification method.
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

            switch class(model)
                case 'VAR'
                    B = model.coeffs(true);
                    SigmaU = model.SigmaU;
                    p = model.p;
                    effects = Recursive.identifyVARTransmission_(B, SigmaU, shockIdx, condition, orderIdx, maxHorizon);
                otherwise
                    error(class(model) + " not supported.")
            end
        end

    end

end
