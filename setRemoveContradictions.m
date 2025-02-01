function setRemoveContradictions(bool)
    % `setRemoveContradictions` Set the global flag for removing contradictions in transmission conditions.
    %
    %   `setRemoveContradictions(bool)` controls whether contradictions in transmission 
    %   conditions are removed when calling `setRemoveContradictions`.
    %
    %   ## Arguments
    %   - `bool` (logical): If `true`, contradictions of the form "x_i & !x_i" will 
    %     be removed from transmission conditions. If `false`, contradictions will 
    %     not be removed.
    %
    %   ## Example
    %   ```
    %   % Enable contradiction removal
    %   setRemoveContradictions(true);
    %
    %   % Disable contradiction removal
    %   setRemoveContradictions(false);
    %   ```
    %   
    %   See also `removeContradictions`
    global REMOVECONTRADICTIONS
    REMOVECONTRADICTIONS = bool;
end
