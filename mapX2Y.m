function strY = mapX2Y(strX, order)
    % `mapX2Y` Map a transmission condition from systems form variables `x` to the variables `y` of the dynamic system.
    %
    %   `strY = mapX2Y(strX, order)` replaces every systems form variable `xi`
    %   in `strX` by the corresponding variable `y_{i, t}` of the dynamic
    %   system. This is the reverse of `mapY2X`.
    %
    %   ## Arguments
    %   - `strX` (string): Transmission condition using systems form variables
    %     `x`, e.g. `'x8 & x4'`.
    %   - `order` (vector): Ordering of variables given by the transmission matrix.
    %
    %   ## Returns
    %   - `strY` (string): The same transmission condition using the variables
    %     `y` of the dynamic system.
    %
    %   ## Example
    %   ```
    %   order = [3, 1, 2];
    %   mapX2Y('x8 & x4', order)  % returns 'y_{1, 2} & y_{3, 1}'
    %   ```
    %
    %   See also `mapY2X`, `mapX2YInt`, `makeConditionY`

    strY = strX;
    pattern = 'x(\d+)';
    matches = regexp(strX, pattern, 'tokens');
    
    for k = 1:length(matches)
        xi = str2double(matches{k}(1));
        [i, t] = mapX2YInt(xi, order);
        strY = regexprep(strY, sprintf('x%i', xi), ['y_{', num2str(i), ', ', num2str(t), '}']);
    end
end
