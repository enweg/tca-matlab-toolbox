function s = stringAnd(s1, s2)
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

