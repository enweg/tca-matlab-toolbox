function s = stringAnd(s1, s2)
    % `stringAnd` Combine two Boolean conditions with AND.
    %
    %   `s = stringAnd(s1, s2)` joins `s1` and `s2` using `' & '`. Each
    %   variable occurs only once in the result, and variables are sorted in
    %   descending order, so negated variables come after all non-negated
    %   variables. If one of the inputs is `'T'` (true), the other input is
    %   returned.
    %
    %   ## Arguments
    %   - `s1` (string): Boolean condition using systems form variables `x`.
    %   - `s2` (string): Boolean condition using systems form variables `x`.
    %
    %   ## Returns
    %   - `s` (string): The combined condition.
    %
    %   ## Example
    %   ```
    %   stringAnd('x1 & x2', 'x1 & !x2')  % returns 'x2 & x1 & !x2'
    %   ```
    %
    %   ## Notes
    %   - This function is intended for internal use only.

    if isequal(s2, 'T')
        s = s1;
        return;
    end
    if isequal(s1, 'T')
        s = s2; 
        return;
    end
    
    % Combine the strings with ' & '
    combined = strjoin({s1, s2}, ' & ');
    
    % Extract unique matches of the pattern `(!{0,1}x\d+)`
    pattern = '!?x\d+';
    matches = regexp(combined, pattern, 'match');
    
    % Sort the matches in reverse order
    uniqueMatches = unique(matches, 'stable');
    % same as sort descend but matlab does not support sort descend for 
    % character cell arrays
    sortedMatches = flip(sort(uniqueMatches));
    
    % Combine the sorted matches with ' & '
    s = strjoin(sortedMatches, ' & ');
end

