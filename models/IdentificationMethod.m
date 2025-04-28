classdef (Abstract) IdentificationMethod 
    % `IdentificationMethod` Abstract class for structural identification.
    %
    %   This abstract class specifies the interface for identification
    %   methods that recover structural models and IRFs from reduced-form
    %   models.
    %
    %   ## Required Methods
    %   1. `irfs = identifyIrfs(obj, model, maxHorizon)`
    %      - Identifies IRFs from the reduced form `model`.
    %      - Returns a 3D array with dimensions:
    %          (response variable, shock, horizon).
    %      - IRFs should be computed up to `maxHorizon`.
    %
    %   2. `[varargout] = identify(obj, model)`
    %      - Identifies the structural form of a reduced form `model`.
    %      - For SVARs from VARs:
    %          - `varargout{1}` = A0 (contemporaneous matrix)
    %          - `varargout{2}` = APlus (lag polynomial matrix)
    %      - For LPs:
    %          - `varargout{1}` = coefficient estimates per horizon.
    %
    %   ## Notes
    %   - See the `Recursive` class for an example implementation.
    %
    %   See also `Recursive`, `SVAR`, `LP`
end
