classdef (Abstract) Model < handle
    methods (Abstract)
        isFitted(obj);
        coeffs(obj);
        fitted(obj);
        residuals(obj);
        nobs(obj);
        getDependent(obj);
        getIndependent(obj);
        getInputData(obj);
        isStructural(obj);
        fit(obj);
        fitAndSelect(obj);
    end

    methods
        function requireFitted(obj)
            if ~isFitted(obj)
                errId = class(obj) + ":NotFitted";
                errMsg = class(obj) + " has not been estimated.";
                error(errId, errMsg);
            end
        end

        function varsIdx = vars2idx_(obj, vars)
            if ischar(vars)
                vars = {vars};
            end
            if ~iscell(vars) && ~isnumeric(vars)
                error("Variables must either be given using indices or variable names.");
            end
            if iscell(vars) && ~ischar(vars{1})
                error("If variables are given by names, they must be `char`.");
            end

            if isnumeric(vars)
                varsIdx = vars;
                return;
            end

            varnames = obj.getVariableNames();
            varsIdx = zeros(length(vars), 1);
            for i = 1:length(vars)
                idx = find(cellfun(@(c) isequal(c, vars{i}), varnames), 1, 'first');
                if isempty(idx)
                    error(vars{i} + " is not a valid variable name.");
                end
                varsIdx(i) = idx;
                % varsIdx(i) = find(varnames == vars, 1, 'first');
            end
        end

        % TODO: test all three functions below
        function q = notThrough(obj, vars, horizons, order)
            varsIdx = obj.vars2idx_(vars);
            orderIdx = obj.vars2idx_(order);
            q = notThrough(varsIdx, horizons, orderIdx);
        end
        function q = through(obj, vars, horizons, order)
            varsIdx = obj.vars2idx_(vars);
            orderIdx = obj.vars2idx_(order);
            q = through(varsIdx, horizons, orderIdx);
        end
        function orderIdx = defineOrder(obj, order)
            % This is just a more user friendly name.
            orderIdx = obj.vars2idx_(order);
        end
        % TODO: how can we handle `transmission`? Do we need to code it for all
        % models separately? 
        
    end
end
