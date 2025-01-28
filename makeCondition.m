function q = makeCondition(s)
    matches = regexp(s, '(x\d+)', 'match');
    vars = unique(matches);
    
    for i = 1:length(vars)
        v = vars{i};
        eval(sprintf('%s = Q(''%s'');', v, v)); 
    end
    
    q = eval(s);
end
