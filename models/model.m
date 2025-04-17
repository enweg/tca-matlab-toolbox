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
    end
end
