classdef (Abstract) Model
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
            if !isFitted(obj)
                error(class(obj) + "must first be estimated.")
            end
        end
    end
end
