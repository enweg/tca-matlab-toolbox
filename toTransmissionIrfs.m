function irfs = toTransmissionIrfs(irfs)
    [k, ~, maxHorizon] = size(irfs);
    maxHorizon = maxHorizon - 1;
    
    % Reshape the 3D array into a 2D matrix by stacking along the third dimension
    irfs = reshape(permute(irfs, [1, 3, 2]), k * (maxHorizon + 1), []);
    
    % Construct the output matrix
    irfsNew = [];
    for h = 0:maxHorizon
        padding = zeros(k * h, k);
        truncated = irfs(1:end - k * h, :);
        irfsNew = [irfsNew, [padding; truncated]];
    end
    
    irfs = irfsNew;
end
