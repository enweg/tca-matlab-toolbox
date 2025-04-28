classdef IRFContainer
    % `IRFContainer` Container for storing impulse response functions.
    %
    %   This class holds IRFs, variable names, the originating model,
    %   and the identification method used (if any).
    %
    %   ## Properties
    %   - `irfs` (3D array): IRFs with dimensions:
    %       - First: Response variable.
    %       - Second: Shock.
    %       - Third: Horizon.
    %   - `varnames` (cell array): Names of response variables.
    %   - `model` (Model): Model used to compute the IRFs.
    %   - `identificationMethod` (IdentificationMethod):
    %       - Identification method if the model is reduced form.

    properties
       irfs                 % 3d array of IRFs with diemensions (response variable, shock, horizon)
       varnames             % variable names for the response variables
       model                % model used to compute the IRFs
       identificationMethod % identification method if `model` is reduced form and IRFs are structural
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

        function irfs = getIrfArray(obj)
            irfs = obj.irfs; 
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
