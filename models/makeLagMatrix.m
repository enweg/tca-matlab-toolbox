function lagMatrix = makeLagMatrix(Y, p)
    [T, k] = size(Y);
    lagMatrix = nan(T, k * p);
    for l = 1:p
        lagMatrix((l+1):end, ((l-1)*k+1):(l*k)) = Y(1:(end-l), :);
    end
end
