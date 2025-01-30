function q = makeCondition(s)
    % makeCondition Create a transmission condition from a Boolean string.
    %
    %   q = makeCondition(s) constructs a transmission condition Q(b) from a 
    %   Boolean statement given as a string. The Boolean expression should use 
    %   variables of the systems form `x<num>`, where `<num>` represents a 
    %   variable index.
    %
    %   Arguments:
    %   - s (string): A Boolean condition string where variables must be 
    %     represented using `x<num>`.
    %
    %   Returns:
    %   - q (Q): A transmission condition.
    %
    %   Example:
    %   s = "x2 & !x3";
    %   cond = makeCondition(s);
    %
    %   Notes:
    %   - Boolean conditions can include AND (&), NOT (! or ~), OR (|), and 
    %     parentheses.
    %   - The resulting transmission condition can be used in transmission to 
    %     calculate the transmission effect.
    %
    %   See also transmission, makeConditionY

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
