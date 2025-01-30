function q = makeCondition(s)
    matches = regexp(s, '(x\d+)', 'match');
    vars = unique(matches);
    
    for i = 1:length(vars)
        v = vars{i};
        eval(sprintf('%s = Q(''%s'');', v, v)); 
    end
    
    % matlab uses ~ for not, while we use the more standard !. 
    % Users could just use ~ and it would work, but we will make sure that !
    % also works
    q = eval(strrep(s, '!', '~'));
end
