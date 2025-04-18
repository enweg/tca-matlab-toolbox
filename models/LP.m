classdef LP < handle & Model
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
            requireFitted(obj);
            if nargin < 2
                excludeDeterministic = false;
            end

            if ~excludeDeterministic
                B = obj.B;
                return;
            end

            B = obj.B(:, (obj.includeConstant+1):end, :)
        end

        function Yhat = fitted(obj)
            requireFitted(obj);
            Yhat = obj.Yhat;
        end

        function U = residuals(obj)
            requireFitted(obj);
            U = obj.U;
        end

        function n = nobs(obj)
            n = size(obj.Y, 1) - obj.p - obj.horizons;
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
            irfs = obj.coeffs(true)
            data = obj.getInputData();
            idxTreatment = findVariableIndex(data, obj.treatment);
            irfs = irfs(:, idxTreatment, :);
            varnames = data.Properties.VariableNames;
            irfObj = IRFContainer(irfs, varnames, obj, opts.identificationMethod);
        end
    end
end
