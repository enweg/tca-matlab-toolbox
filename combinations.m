function combs = combinations(arr)
    % provides all possible combinations of length 1 through K=length(arr)
    % of the elements in arr

    combs = {};
    for kk = 1:length(arr)
        combs = [combs; num2cell(nchoosek(arr, kk), 2)];
    end
end
