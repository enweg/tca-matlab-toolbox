classdef SVAR < handle & Model 
    properties
        A0
        APlus
        p
        trendExponents
        VARModel
    end

    methods (Static)
        function Y = simulate(shocksOrT, A0, APlus, varargin)
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
            % APlus exlucludes coefficients for deterministic components.
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
            requireFitted(obj);
            Yhat = obj.VARModel.fitted();
        end

        function U = residuals(obj)
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
            requireFitted(obj);
            C = obj.VARModel.makeCompanionMatrix();
        end

        function rho = spectralRadius(obj)
            requireFitted(obj);
            rho = obj.VARModel.spectralRadius();
        end

        function flag = isStable(obj)
            requireFitted(obj);
            flag = obj.VARModel.isStable();
        end

        function val = aic(obj)
            requireFitted(obj);
            val = obj.VARModel.aic();
        end
        function val = bic(obj)
            requireFitted(obj);
            val = obj.VARModel.bic();
        end
        function val = sic(obj)
            requireFitted(obj);
            val = obj.VARModel.sic();
        end
        function val = hqc(obj)
            requireFitted(obj);
            val = obj.VARModel.hqc();
        end

        function fit(obj, identificationMethod)
            obj.VARModel.fit();
            [A0, APlus] = identificationMethod.identify(obj.VARModel);
            obj.A0 = A0;
            obj.APlus = APlus;
        end

        function [modelBest, icTable] = fitAndSelect(obj, identificationMethod, icFunction)
            if nargin < 3
                icFunction = @VAR.aic_;
            end

            [modelVARBest, icTable] = fitAndSelect(obj.VARModel, icFunction);
            [A0, APlus] = identificationMethod.identify(modelVARBest);

            data = obj.getInputData();
            modelBest = SVAR(data, modelVARBest.p, 'trendExponents', obj.trendExponents, 'A0', A0, 'APlus', APlus, 'VARModel', modelVARBest);
        end

        function irfObj = IRF(obj, maxHorizon)
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

        % TODO: test
        function effects = transmission(obj, shock, condition, order, maxHorizon, varargin)
            opts.identificationMethod = missing
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
