function s = stringAnd(s1, s2)
    if s2 == 'T'
        s = s1;
        return;
    end
    if s1 == 'T'
        s = s2; 
        return;
    end
    
    % Combine the strings with ' & '
    combined = strjoin({s1, s2}, ' & ');
    
    % Extract unique matches of the pattern `(!{0,1}x\d+)`
    pattern = '!?x\d+';
    matches = regexp(combined, pattern, 'match');
    
    % Sort the matches in reverse order
    unique_matches = unique(matches, 'stable');
    % same as sort descend but matlab does not support sort descend for 
    % character cell arrays
    sorted_matches = flip(sort(unique_matches));
    
    % Combine the sorted matches with ' & '
    s = strjoin(sorted_matches, ' & ');
end

