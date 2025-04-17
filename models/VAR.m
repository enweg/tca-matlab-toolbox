classdef VAR < handle & Model
    properties
        B                   % Coefficient matrices [C, B_1, ..., B_p]
        SigmaU              % Covariance of error term
        p                   % Lag-length
        trendExponents      % Deterministic trend exponents for time trend

        inputData           % Data used to fit the model
        Y                   % LHS matrix
        X                   % RHS matrix
        U                   % Residual Matrix
        Yhat                % Fitted values
    end

    methods(Static)
        function val = ic_(SigmaU, nCoeffs, ct)
            val = logdet(SigmaU) + ct * nCoeffs; 
        end
        function val = aic_(SigmaU, nCoeffs, T)
            ct = 2 / T;
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = hqc_(SigmaU, nCoeffs, T)
            ct = 2 * log(log(T)) / T;
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = sic_(SigmaU, nCoeffs, T)
            ct = log(T) / T;
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = bic_(SigmaU, nCoeffs, T)
            val = VAR.sic_(SigmaU, nCoeffs, T);
        end

        function Y = simulate(errorsOrT, B, varargin)

            % Default values
            opts.trendExponents = [0];
            opts.initial = [];  % will become zeros as default
            opts.SigmaU = [];   % will become identity matrix by default
            % Parsing user options (name-value pairs)
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error("VAR.simulate: " + varargin{i} + " is not a valid simulation option.");
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            if isscalar(errorsOrT)
                errors = randn(size(B, 1), errorsOrT);
                if ~isempty(opts.SigmaU)
                    errors = chol(opts.SigmaU)' * errors;
                end
            elseif ismatrix(errorsOrT)
                errors = errorsOrT;
            else
                error("VAR.simulate: First argument must either be a k*T matrix of errors terms or an integer denoting the number of periods to simulate.")
            end

            [k, T] = size(errors);
            m = length(opts.trendExponents);
            kp = size(B, 2) - m;

            if mod(kp, k) ~= 0
                error("VAR.simulate: Dimensions of B are wrong.")
            end

            if isempty(opts.initial)
                opts.initial = zeros(kp, 1);
            end

            Y = zeros(size(errors));
            Zt = ones(m + kp, 1);
            Zt((m+1):end) = opts.initial;
            for t = 1:T
                % trend terms
                Zt(1:m) = arrayfun(@(x) t^x, opts.trendExponents);
                Y(:, t) = B * Zt + errors(:,  t);
                if kp > 0
                    % Rotating in the new data point into Zt
                    Zt((m+1+k):end) = Zt((m+1):(end-k));
                    Zt((m+1):(m+k)) = Y(:, t);
                end
            end

            Y = Y';
        end

        function irfs = IRF_(B, p, maxHorizon)
            % B must exclude coefficients for deterministic components
            k = size(B, 1);
            irfs = zeros(k, k, maxHorizon + 1);
            irfs(:, :, 1) = eye(k);

            for h = 1:maxHorizon
                for j = 1:min(h, p)
                    Bj = B(:, ((j-1)*k+1):(j*k));
                    irfs(:, :, h+1) = irfs(:, :, h+1) + Bj * irfs(:, :, h+1-j);
                end
            end
        end

        function C = makeCompanionMatrix_(B, p, m)
            % m is the number of exogenous components
            B = B(:, (m+1):end);
            k = size(B, 1);
            C = diag(ones(k * p - k, 1), -k);
            C(1:k, :) = B;
        end

        function rho = spectralRadius_(C)
            rho = max(abs(eig(C)));
        end
    end

    methods

        function obj = VAR(data, p, varargin)
            % Default values
            opts.trendExponents = [0];
            opts.B = [];
            opts.SigmaU = [];

            % Parsing options (name-value pairs)
            for i=1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error("VAR: " + varargin{i} + " is not a valid option.")
                end

                opts.(varargin{i}) = varargin{i+1};
            end

            if istable(data)
                dataMatrix = table2array(data);
            else
                dataMatrix = data;
                varnames = arrayfun(@(i) "Y" + i, 1:size(data, 2));
                data = array2table(data, 'VariableNames', varnames);
            end

            U = [];
            Yhat = [];
            Y = dataMatrix((p+1):end, :);
            X = makeLagMatrix(dataMatrix, p);
            X = X((p+1):end, :);
            for i = 1:length(opts.trendExponents)
                te = opts.trendExponents(i);
                trend = arrayfun(@(x) x^te, (p+1):size(dataMatrix, 1))';
                X = [trend X];
            end

            obj.B = opts.B;
            obj.SigmaU = opts.SigmaU;
            obj.p = p;
            obj.trendExponents = opts.trendExponents;
            obj.inputData = data;
            obj.Y = Y;
            obj.X = X;
            obj.U = U;
            obj.Yhat = Yhat;
        end

        function flag = isFitted(obj)
            if size(obj.Yhat) > 0
                flag = true;
                return;
            end
            flag = false;
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

            m = length(obj.trendExponents);
            B = obj.B(:, (m+1):end);
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
            n = size(obj.Y, 1);
        end

        function n = ncoeffs(obj)
            requireFitted(obj);
            n = prod(size(coeffs(obj)));
        end

        function Y = getDependent(obj)
            Y = obj.Y;
        end

        function X = getIndependent(obj)
            X = obj.X;
        end

        function data = getInputData(obj)
            data = obj.inputData;
        end

        function varnames = getVariableNames(obj)
            data = getInputData(obj);
            varnames = data.Properties.VariableNames;
        end

        function flag = isStructural(obj)
            flag = false;
        end

        function C = makeCompanionMatrix(obj)
            requireFitted(obj);
            m = length(obj.trendExponents);
            C = VAR.makeCompanionMatrix_(coeffs(obj), obj.p, m); 
        end

        function rho = spectralRadius(obj)
            requireFitted(obj);
            rho = VAR.spectralRadius_(makeCompanionMatrix(obj));
        end

        function flag = isStable(obj)
            requireFitted(obj);
            flag = (spectralRadius(obj) < 1);
        end

        function val = aic(obj)
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.aic_(SigmaU, nCoeffs, T);
        end
        function val = hqc(obj)
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.hqc_(SigmaU, nCoeffs, T);
        end
        function val = sic(obj)
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.sic_(SigmaU, nCoeffs, T);
        end
        function val = bic(obj)
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.bic_(SigmaU, nCoeffs, T);
        end

        function fit(obj)
            X = obj.X;
            Y = obj.Y;

            obj.B = (Y' * X) / (X' * X);
            obj.Yhat = X * obj.B';
            obj.U = Y - obj.Yhat;
            obj.SigmaU = obj.U' * obj.U / nobs(obj);
        end

        function [modelBest, icTable] = fitAndSelect(obj, icFunction)
            if nargin < 2
                icFunction = @VAR.aic_;
            end

            pMax = obj.p;
            ps = 0:pMax;
            ics = zeros(length(ps), 1);

            obj.fit();
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            T = nobs(obj);
            icBest = icFunction(SigmaU, nCoeffs, T);
            ics(pMax+1) = icBest;
            modelBest = obj;

            for p = (pMax-1):-1:0
                modelTmp = VAR(getInputData(obj), p, 'trendExponents', obj.trendExponents);
                modelTmp.fit();
                nCoeffsTmp = ncoeffs(modelTmp);
                UTmp = residuals(modelTmp);
                UTmp = UTmp((pMax-p+1):end, :);
                SigmaUTmp = UTmp' * UTmp / T;

                icTmp = icFunction(SigmaUTmp, nCoeffsTmp, T);
                ics(p+1) = icTmp;
                if icTmp < icBest
                    icBest = icTmp;
                    modelBest = modelTmp;
                end
            end

            icTable = table(ps', ics, 'VariableNames', {'p', 'IC'});
        end

        function irfObj = IRF(obj, maxHorizon, varargin)
            opts.identificationMethod = missing;
            for i = 1:2:length(varargin)
                if ~isfield(opts, varargin{i})
                    error("VAR.IRF: " + varargin{i} + " is not a valid option.")
                end
                opts.(varargin{i}) = varargin{i+1};
            end

            requireFitted(obj);
            if ismissing(opts.identificationMethod)
                B = obj.coeffs(true);
                irfs = VAR.IRF_(B, obj.p, maxHorizon);
            else
                % This way users can easily implement new identification methods.
                irfs = opts.identificationMethod.identifyIrfs(obj, maxHorizon);
            end
            varnames = getVariableNames(obj);
            irfObj = IRFContainer(irfs, varnames, obj, opts.identificationMethod);
        end
    end
end
