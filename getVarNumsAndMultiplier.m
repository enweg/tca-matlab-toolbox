function [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q)

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

