function leadMatrix = makeLeadMatrix(Y, p)
    [T, k] = size(Y);
    leadMatrix = nan(T, k * p);
    for l = 1:p
        leadMatrix(1:(end-l), ((l-1)*k+1):(l*k)) = Y((l+1):end, :);
    end
end
