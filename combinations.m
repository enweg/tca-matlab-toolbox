function combs = combinations(arr)
    % combinations Generate all possible subsets of an array.
    %
    %   combs = combinations(arr) returns all possible combinations of elements 
    %   in arr, for subset lengths ranging from 1 to K, where K is the length of arr.
    %
    %   Arguments:
    %   - arr (vector): Input array of elements.
    %
    %   Returns:
    %   - combs (cell array): A cell array where each element is a subset of arr 
    %     with varying lengths (from 1 to length(arr)).
    %
    %   Example:
    %   arr = [1, 2, 3];
    %   combs = combinations(arr);
    %   % Output:
    %   %   { [1], [2], [3], [1 2], [1 3], [2 3], [1 2 3] }

    combs = {};
    for kk = 1:length(arr)
        combs = [combs; num2cell(nchoosek(arr, kk), 2)];
    end
end
