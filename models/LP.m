classdef LP < handle & Model
    properties
        data 
        treatment
        p
        horizons
        includeConstant

        coeffs
        Y
        X
        U
        Yhat
    end

    methods (Static)
        function [X, Y] = createXY_(data, treatment, p, horizons, varargin)
            opts.includeConstant = true
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
                X = [ones(T-p, 1), X]
            end

            maxHorizon = max(horizons);
            Y = makeLeadMatrix(dataMatrix, maxHorizon);
            Y = [dataMatrix Y];
            Y = Y((p+1):end, :);
            Y = reshape(Y, T-p, k, []);
            Y = Y(:, :, horizons + 1)
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

            [Y, X] = LP.createXY_(data, treatment, p, horizons, opts.includeContant);

            if ~istable(data)
                data = array2table(data, 'VariableNames', arrayfun(@(i) "Y" + i, 1:size(data, 2)));
            end
            obj.data = data; 
            obj.treatment = treatment; 
            obj.p = p;
            obj.horizons = horizons; 
            obj.includeConstant = opts.includeConstant;
            obj.coeffs = [];
            obj.Y = Y;
            obj.X = X;
            obj.U = [];
            obj.Yhat = [];
        end

        function flag = isFitted(obj)
            flag = size(obj.Yhat) > 0
        end

        function B = coeffs(obj, excludeDeterministic)
            requireFitted(obj);
            if nargin < 2
                excludeDeterministic = false;
            end

            if ~excludeDeterministic
                B = obj.coeffs;
                return;
            end

            B = obj.coeffs(:, (obj.includeConstant+1):end, :)
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

        % TODO: implement fit, fitAndSelect, IRF

    end
end
