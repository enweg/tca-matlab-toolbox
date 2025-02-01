classdef Q
    % `Q` Represents a transmission condition.
    %
    %   The `Q` class is used to define transmission conditions based on Boolean 
    %   statements. A transmission condition is denoted as $Q(b)$, where $b$ is a 
    %   Boolean statement involving variables `x<num>`, such as `x1`, `x2`, etc.
    %   Boolean statements should only contain AND (&) and NOT (!) operations.
    %
    %   ## Properties
    %   - `vars` (cell array of strings): Contains the Boolean variable expressions.
    %   - `multiplier` (vector of numbers): Multipliers associated with each term.
    %
    %   ## Methods
    %   - `Q`: Constructor to create a transmission condition.
    %   - `and`: Overloads `&` (logical AND) for `Q` objects.
    %   - `or`: Overloads `|` (logical OR) for `Q` objects.
    %   - `not`: Overloads `~` (logical NOT) for `Q` objects.
    %   - `disp`: Custom display function.
    %   - `display`: Calls `disp` for better formatting.
    %
    %   ## Usage
    %   ```
    %   % Define variables as transmission conditions
    %   x = arrayfun(@(i) Q(sprintf('x%d', i)), 1:10);
    %   q = (x(1) | x(2)) & ~x(3);
    %
    %   % Alternatively, define variables separately
    %   x1 = Q('x1');
    %   x2 = Q('x2');
    %   x3 = Q('x3');
    %   q = (x1 | x2) & ~x3;
    %
    %   % Creating Q objects with multipliers
    %   q = Q('x1 & !x3', 1);
    %   q = Q({'x1', 'x2', 'x1 & x2'}, [1, 1, -1]);
    %   ```
    %
    %   ## Notes
    %   - The recommended constructor is `Q(i)`, where `i` is an integer representing a variable index.
    %   - Other constructors are for internal use and may lead to incorrect results if misused.
    %   - DO NOT use OR (`|`) inside the string input for `Q`, as it is not supported.

    properties
        vars
        multiplier
    end
    
    methods
        % Constructor
        function obj = Q(vars, multiplier)
            % `Q` Construct a transmission condition.
            %
            %   `obj = Q(vars)` constructs a transmission condition with the given variable.
            %   `obj = Q(vars, multiplier)` constructs a transmission condition with 
            %   a specified multiplier.
            %
            %   ## Arguments
            %   - `vars` (string, cell array of strings, or integer): The variable(s) 
            %     in the Boolean condition. Must be formatted as `x<num>`.
            %   - `multiplier` (number or vector): Multiplier(s) associated with 
            %     each term. 
            %
            %   ## Returns
            %   - `obj` (Q): A transmission condition.
            %
            %   ## Example
            %   ```
            %   q = Q('x1');            % Single variable
            %   q = Q({'x1', 'x2'}, [1, -1]);  % Multiple variables with multipliers
            %   ```
            %
            %   ## Notes
            %   - The recommended way to define a variable is using `Q(i)`, where `i` 
            %     is an integer representing a variable index.
            %   - Users should not directly specify OR (`|`) inside the variable strings.

            if nargin == 0
                obj.vars = {};
                obj.multiplier = [];
            elseif nargin == 1
                if ischar(vars)
                    obj.vars = {vars};
                    obj.multiplier = 1.0;
                elseif isnumeric(vars) & length(vars) == 1
                    obj.vars = {sprintf('x%d', vars)};
                    obj.multiplier = 1.0;
                elseif iscell(vars) & (ischar(vars{1}) | isstring(vars{1}))
                    obj.vars = vars;
                    obj.multiplier = ones(size(vars));
                elseif isvector(vars) & isnumeric(vars(1))
                    obj.vars = arrayfun(@(x) sprintf("x%d", x), vars, 'UniformOutput', false);
                    obj.multiplier = ones(size(vars));
                else
                    error('Invalid input type');
                end
            elseif nargin == 2
                if ischar(vars)
                    obj.vars = {vars};
                    obj.multiplier = multiplier;
                elseif isnumeric(vars) & length(vars) == 1
                    obj.vars = {sprintf('x%d', vars)};
                    obj.multiplier = multiplier;
                elseif iscell(vars) & (ischar(vars{1}) | isstring(vars{1}))
                    obj.vars = vars;
                    obj.multiplier = multiplier;
                elseif isvector(vars) & isnumeric(vars(1))
                    obj.vars = arrayfun(@(x) sprintf("x%d", x), vars, 'UniformOutput', false);
                    obj.multiplier = multiplier;
                else
                    error('Invalid input type');
                end
            end
        end

        % Overload the & operator (logical AND)
        function result = and(obj1, obj2)
            % `&` combine two transmission conditions using logical AND.
            %
            %   `result = q1 & q2` performs a logical AND operation between two 
            %   transmission conditions, returning a new `Q` object. This is the 
            %   same as $Q(b \land b')$ where $b$ and $b'$ are the Boolean conditions for
            %   `q1` and `q2` respectively.
            %
            %   ## Arguments
            %   - `obj1` (Q): First transmission condition.
            %   - `obj2`(Q): Second transmission condition.
            %
            %   ## Returns
            %   - `result`(Q): The combined transmission condition.
            %
            %   Example:
            %   ```
            %   q1 = Q(1);
            %   q2 = Q(2);
            %   q = q1 & q2;  % Returns Q("x2 & x1")
            %   ```
            %
            %   See also `or` (`|`), `not` (`~`)
            if isa(obj2, 'Q')
              if length(obj1.vars) == 1
                % Perform operation when obj1 has only one variable
                q = collectTerms(Q( ...
                    arrayfun(@(i) stringAnd(obj2.vars{i}, obj1.vars{1}), 1:length(obj2.vars), 'UniformOutput', false), ...
                    obj1.multiplier(1) * obj2.multiplier));
                result = removeContradictions(q);
              elseif length(obj2.vars) == 1
                % Perform operation when obj2 has only one variable
                q = collectTerms(Q( ...
                    arrayfun(@(i) stringAnd(obj1.vars{i}, obj2.vars{1}), 1:length(obj1.vars), 'UniformOutput', false), ...
                    obj2.multiplier(1) * obj1.multiplier));
                result = removeContradictions(q);
              else
                % General case when both have multiple variables
                qs = arrayfun(@(i) obj1 & Q(obj2.vars{i}, obj2.multiplier(i)), 1:length(obj2.vars), 'UniformOutput', false);
                vars = {};
                mults = [];
                for ii=1:numel(qs)
                  vars = [vars, qs{ii}.vars];
                  mults = [mults, qs{ii}.multiplier];
                end
                result = collectTerms(Q(vars, mults));
              end
            else
                error('Operand must be an instance of class Q');
            end
        end

        % Overload the | operator (logical OR)
        function result = or(obj1, obj2)
            % `|` Combine two transmission conditions using a logical OR.
            %
            %   `result = q1 | q2` performs a logical OR operation between two 
            %   transmission conditions, returning a new `Q` object.
            %
            %   ## Arguments
            %   - `obj1` (Q): First transmission condition.
            %   - `obj2` (Q): Second transmission condition.
            %
            %   ## Returns
            %   - `result` (Q): The combined transmission condition.
            %
            %   ## Example
            %   ```
            %   q1 = Q(1);
            %   q2 = Q(2);
            %   q = q1 | q2;
            %   ```
            %
            %   See also `and` (`&`), `not` (`~`)
            vars = [obj1.vars, obj2.vars];
            q = obj1 & obj2;
            vars = [vars, q.vars];
            mults = [obj1.multiplier, obj2.multiplier];
            mults = [mults, -1 * q.multiplier];
            result = collectTerms(Q(vars, mults));
        end

         % Overload the ~ operator (logical NOT)
        function result = not(obj)
            % `~` Negate a boolean condition using the logican NOT.
            %
            %   `result = ~q` negates a transmission condition, creating a condition 
            %   where the given Boolean statement does not hold.
            %
            %   ## Arguments
            %   - `obj` (Q): A transmission condition to negate.
            %
            %   ## Returns
            %   - `result` (Q): The negated transmission condition.
            %
            %   ## Example
            %   ```
            %   q1 = Q(1);
            %   q = ~q1; 
            %   ```
            %
            %   ## Notes
            %   - If the condition consists of a single variable, it is simply negated.
            %   - If the condition is more complex, an auxiliary `"T"` (true) condition 
            %     is used and the returned condition is equivalent to $Q(T) - Q(b)$ 
            %     where $b$ is the origional Boolean condition.
            %
            %   See also `and` (`&`), `or` (`|`)

            if length(obj.vars) == 1 && length(regexp(obj.vars{1}, 'x\d+', 'match')) == 1
                % If there's a single variable and it matches the pattern "x<number>"
                vars = {['!', obj.vars{1}]}; % Prepend "!" to the variable
                result = Q(vars, obj.multiplier(1));
            else
                % General case
                vars = [{'T'}, obj.vars]; % T denotes true
                mults = [1.0, -1 * obj.multiplier]; 
                result = collectTerms(Q(vars, mults));
            end
        end

        % Overloading the display option
        function disp(obj, order)
            s = "";
            for i = 1:length(obj.multiplier)
                m = obj.multiplier(i);
                v = obj.vars{i};

                mSign = sign(m);
                if mod(m, 1) == 0
                    ms = num2str(floor(abs(m)));
                else
                    ms = num2str(abs(m));
                end
                if abs(m) == 1
                    ms = "";
                end

                vs = string(v);
                sTmp = sprintf("%sQ(%s)", ms, vs);

                if isequal(s, "") && mSign < 0
                    s = join("-", sTmp, "");
                elseif isequal(s, "") && mSign > 0
                    s = sTmp;
                elseif mSign < 0
                    s = join([s, sTmp], " - ");
                else 
                    s = join([s, sTmp], " + ");
                end
            end

            if nargin == 2
                s = mapX2Y(s, order);
                disp(s);
            else
                disp(s);
            end
        end
        function display(obj)
            disp(obj);
        end
    end
end

