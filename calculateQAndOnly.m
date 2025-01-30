function effect = calculateQAndOnly(from, irfs, irfsOrtho, vars, multiplier)
    % calculateQAndOnly Compute the transmission effect of a condition with only ANDs.
    %
    %   result = calculateQAndOnly(from, irfs, irfsOrtho, vars, multiplier) 
    %   calculates the transmission effect of a transmission condition/query 
    %   that involves only AND operations.
    %
    %   Arguments:
    %   - from (integer): Index of the shock.
    %   - irfs (matrix): IRFs in transmission form. See also toTransmissionIrfs.
    %   - irfsOrtho (matrix): Orthogonalized / Cholesky IRFs in transmission form. 
    %     See also toTransmissionIrfs.
    %   - vars (vector of integers): Indices of variables involved in the condition.
    %   - multiplier (number): Multiplier.
    %
    %   Returns:
    %   - result (vector): A vector where entry i corresponds to the transmission 
    %     effect on variable x_i.
    %
    %   Notes:
    %   - This function is intended for internal use only.
    %
    %   See also toTransmissionIrfs

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
