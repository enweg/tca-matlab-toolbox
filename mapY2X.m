function strX = mapY2X(strY, order)
    % `mapY2X` Map a transmission condition from the variables `y` of the dynamic system to systems form variables `x`.
    %
    %   `strX = mapY2X(strY, order)` replaces every variable `y_{i, t}` in
    %   `strY` by the corresponding systems form variable `xj`.
    %
    %   ## Arguments
    %   - `strY` (string): Transmission condition using the variables `y` of
    %     the dynamic system, e.g. `'y_{1, 2} & y_{3, 1}'`.
    %   - `order` (vector): Ordering of variables given by the transmission matrix.
    %
    %   ## Returns
    %   - `strX` (string): The same transmission condition using systems form
    %     variables `x`.
    %
    %   ## Example
    %   ```
    %   order = [3, 1, 2];
    %   mapY2X('y_{1, 2} & y_{3, 1}', order)  % returns 'x8 & x4'
    %   ```
    %
    %   See also `mapX2Y`, `mapY2XInt`, `makeConditionY`

    strX = strY;
    pattern = 'y_{(\d+),\s*(\d+)}';
    matches = regexp(strY, pattern, 'tokens');
    
    for k = 1:length(matches)
        i = str2double(matches{k}{1});
        t = str2double(matches{k}{2});
        xi = mapY2XInt(i, t, length(order), order);
        strX = regexprep(strX, sprintf('y_{%i,\\s*%i}', i, t), ['x', num2str(xi)]);
    end
end
