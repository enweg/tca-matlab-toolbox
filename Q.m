classdef Q
    properties
        vars
        multiplier
    end
    
    methods
        % Constructor
        function obj = Q(vars, multiplier)
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
            vars = [obj1.vars, obj2.vars];
            q = obj1 & obj2;
            vars = [vars, q.vars];
            mults = [obj1.multiplier, obj2.multiplier];
            mults = [mults, -1 * q.multiplier];
            result = collectTerms(Q(vars, mults));
        end

         % Overload the ~ operator (logical NOT)
        function result = not(obj)
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
        function disp(obj)
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
            disp(s);
        end
        function display(obj)
            disp(obj);
        end
    end
end

