classdef IRFContainer
    properties
       irfs
       varnames
       model
       identificationMethod
    end

    methods
        function obj = IRFContainer(irfs, varnames, model, identificationMethod)
            if ~isa(model, 'Model')
                error("The provided model is not a Model.");
            end
            if nargin == 3
                identificationMethod = missing;
            end
            if ~isa(identificationMethod, 'IdentificationMethod') && ~ismissing(identificationMethod)
                error("The provided identification method is not a valid.")
            end
        
            obj.irfs = irfs; 
            obj.varnames = varnames; 
            obj.model = model; 
            obj.identificationMethod = identificationMethod;
        end

        % Customising indexing
        function out = subsref(obj, S)
            switch S(1).type
                case '()'
                    % Allow obj(i) to return irfs(i)
                    out = obj.irfs(S(1).subs{:});
                case '.'
                    % Allow normal property/method access
                    out = builtin('subsref', obj, S);
                otherwise
                    error('Unsupported subscripted reference type');
            end
        end

        function obj = plus(obj, scalar)
            if isnumeric(scalar)
                obj.irfs = obj.irfs + scalar;
            else
                error('Addition only supported with scalars.');
            end
        end

        function obj = minus(obj, scalar)
            if isnumeric(scalar)
                obj.irfs = obj.irfs - scalar;
            else
                error('Subtraction only supported with scalars.');
            end
        end

        function obj = times(obj, scalar)
            if isnumeric(scalar)
                obj.irfs = obj.irfs .* scalar;
            else
                error('Multiplication only supported with scalars.');
            end
        end

        function obj = rdivide(obj, scalar)
            if isnumeric(scalar)
                obj.irfs = obj.irfs ./ scalar;
            else
                error('Division only supported with scalars.');
            end
        end

        function disp(obj)
            disp(obj.irfs);
        end
    end
end
