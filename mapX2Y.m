function strY = mapX2Y(strX, order)
    strY = strX;
    pattern = 'x(\d+)';
    matches = regexp(strX, pattern, 'tokens');
    
    for k = 1:length(matches)
        xi = str2double(matches{k}(1));
        [i, t] = mapX2YInt(xi, order);
        strY = regexprep(strY, sprintf('x%i', xi), ['y_{', num2str(i), ', ', num2str(t), '}']);
    end
end
