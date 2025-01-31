function [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q)
    % `getVarNumsAndMultiplier` Extract AND and NOT expressions from a transmission condition.
    %
    %   `[andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q)` extracts 
    %   the AND and NOT expressions from a transmission condition and returns 
    %   the corresponding variable numbers. It also returns the multiplier for 
    %   each term.
    %
    %   A valid transmission condition consists of terms that involve only AND 
    %   and NOT expressions. For each term, the AND and NOT expressions are 
    %   collected, and vectors of the respective variable numbers are returned.
    %
    %   ## Arguments
    %   - `q` (Q): A transmission condition. See also `Q`.
    %
    %   ## Returns
    %   - `andNums` (cell array of vectors): Each cell contains a vector of variable 
    %     numbers included via AND in the term.
    %   - `andNotNums` (cell array of vectors): Each cell contains a vector of 
    %     variable numbers included via NOT in the term.
    %   - `multiplier` (vector of numbers): Contains the multiplier for each term.
    %
    %   ##Example
    %   ```
    %   q = Q({"x1", "!x2", "x1 & !x2"}, [1, 2, 3]);
    %   [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q);
    %   % Output:
    %   % andNums = { [1], [], [1] }
    %   % andNotNums = { [], [2], [2] }
    %   % multiplier = [1, 2, 3]
    %   ```
    %
    %   See also `Q`

    % Initialize cell arrays for AND and NOT numbers
    andNums = cell(1, numel(q.vars));
    andNotNums = cell(1, numel(q.vars));

    for i = 1:numel(q.vars)
        % Extract numbers matching "x<num>" (positive matches)
        andMatches = regexp(q.vars{i}, '(?<!\!)x(\d+)', 'tokens');
        andNums{i} = cellfun(@(x) str2double(x{1}), andMatches);

        % Extract numbers matching "!x<num>" (negated matches)
        andNotMatches = regexp(q.vars{i}, '!x(\d+)', 'tokens');
        andNotNums{i} = cellfun(@(x) str2double(x{1}), andNotMatches);
    end

    % Extract the multiplier
    multiplier = q.multiplier;
end

