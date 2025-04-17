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
                    % Something
                    B = model.coeffs();
                    SigmaU = model.SigmaU;
                    [A0, APlus] = Recursive.identifyVAR_(B, SigmaU);
                    varargout{1} = A0;
                    varargout{2} = APlus;
                otherwise
                    error("Recursive identification of " + class(model) + " is not implemented.");
            end
        end

    end

end
