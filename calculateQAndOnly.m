function effect = calculateQAndOnly(from, irfs, irfsOrtho, vars, multiplier)
    
    if isempty(vars)
        % indicating TRUE
        effect = irfs(:, from) * multiplier;
        return;
    end
    
    vars = sort(vars);
    
    effect = zeros(size(irfs, 1), 1);
    effect(:) = multiplier * irfs(vars(1), from);
    
    for i = 1:(length(vars)-1)
        effect = effect .* (irfsOrtho(vars(i+1), vars(i)) / irfsOrtho(vars(i), vars(i)));
    end
    
    effect = effect .* (irfsOrtho(:, vars(end)) / irfsOrtho(vars(end), vars(end)));
end
