classdef Recursive < IdentificationMethod 
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
