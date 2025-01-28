function strX = mapY2X(strY, order)
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
