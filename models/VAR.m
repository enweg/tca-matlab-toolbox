classdef VAR < Model
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
            ct = 2 * log(log(T)) / T
            val = VAR.ic_(SigmaU, nCoeffs, ct);
        end
        function val = sic_(SigmaU, nCoeffs, T)
            ct = log(T) / T
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

            [k, T] = size(errors);
            m = length(opts.trendExponents);
            kp = size(B, 2) - m;

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

            if mod(kp, k) ~= 0
                error("VAR.simulate: Dimensions of B are wrong.")
            end

            if isempty(opts.initial)
                opts.initial = zeros(kp);
            end

            Y = zeros(size(errors));
            Zt = ones(m + kp);
            Zt((m+1):end) = opts.initial
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
        end

        function irfs = IRF_(B, p, max_horizon)
            % B must exclude coefficients for deterministic components
            k = size(B, 1);
            irfs = zeros(k, k, max_horizon + 1);
            irfs(:, :, 1) = eye(k);

            for h = 1:max_horizon
                for j = 1:min(h, p)
                    Bj = B(:, ((j-1)*k+1):(j*k));
                    irfs(:, :, h+1) = irfs(:, :, h+1) + Bj * irfs(:, :, h+1-j);
                end
            end
        end
    end

    methods
        function flag = isFitted(obj)
            if size(obj.Yhat) > 0
                flag = true;
                return;
            end
            flag = false;
        end

        function B = coeffs(obj, excludeDeterministic)
            if nargin < 2
                excludeDeterministic = false;
            end
            requireFitted(obj);
            if ~excludeDeterministic
                B = obj.B;
                return;
            end

            m = length(obj.trendExponents);
            B = B(:, (m+1):end);
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

        function flag = isStructural(obj)
            flag = false;
        end

        function C = makeCompanionMatrix(obj)
            requireFitted(obj);
            m = length(obj.trendExponents);
            C = makeCompanionMatrix(coeffs(obj), obj.p, m); 
        end

        function rho = spectralRadius(obj)
            requireFitted(obj);
            rho = spectralRadius(makeCompanionMatrix(obj));
        end

        function flag = isStable(obj)
            requireFitted(obj);
            flag = (spectralRadius(obj) < 1);
        end

        function val = aic(obj)
            T = nobs(obj);
            SigmaU = obj.SigmaU;
            nCoeffs = ncoeffs(obj);
            val = VAR.aic_(SigmaU, nCoeffs, T)
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
                icFunction = VAR.aic_;
            end

            pMax = obj.p;
            ps = 0:pMax;
            ics = zeros(length(ps));

            obj.fit();
            SigmaU = obj.SigmaU;
            nCoeffs = ncoefs(obj);
            T = nobs(obj);
            icBest = icFunction(SigmaU, nCoeffs, T);
            ics(pMax+1) = icBest;
            modelBest = obj;

            for p = (pMax-1):-1:0
                modelTmp = VAR(getInputData(obj), p, obj.trendExponents);
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

            icTable = [ps ics];
        end

        % TODO: should return an IRF object !!!!
        function irfs = IRF(obj, max_horizon)
            B = obj.coeffs(true);
            irfs = VAR.IRF_(B, obj.p, max_horizon);
        end

    end
end
