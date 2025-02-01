function q = makeConditionY(strY, order)
    % `makeConditionY` Create a transmission condition from a Boolean string using dyanmic form variables.
    %
    %   `q = makeConditionY(strY, order)` constructs a transmission condition $Q(b)$ 
    %   from a Boolean statement specified in terms of dynamic system variables 
    %   (i.e., `y_{i,t}` notation). 
    %
    %   ## Arguments
    %   - `strY` (string): A Boolean condition string where variables are represented 
    %     as `y_{i,t}`, with `i` as the variable index and `t >= 0` as the time period.
    %     `t=0` corresponds to the contemporaneous horizon.
    %   - `order` (vector of integers): The variable ordering defined by the transmission matrix.
    %
    %   ## Returns
    %   - `q` (Q): A transmission condition object.
    %
    %   ## Example
    %   ```
    %   s_y = "y_{1,0} & !y_{1,1}";
    %   order = [3,1,2];
    %   cond = makeConditionY(s_y, order);
    %   ```
    %
    %   ## Notes
    %   - Boolean conditions can include AND (&), NOT (! or ~), OR (|), and 
    %     parentheses.
    %   - The resulting transmission condition can be used in `transmission` to 
    %     calculate the transmission effect.
    %
    %   See also `transmission`, `makeCondition`

    strX = mapY2X(strY, order);
    q = makeCondition(strX);
end
