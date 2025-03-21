function checkDynare()
    % First check if Dynare is loaded
    if exist("dynare", "file") == 2
        pathDynare = fileparts(which("dynare.m"));
        if isempty(pathDynare)
            error("Could not determine Dynare installation path.");
        end

        if exist("kalman_transition_matrix", "file") == 2
            % Older versions of Dynare have the function in the root directory.
            % The function is therefore loaded as soon as base Dynare is loaded.
            return; 
        end

        % Newer versions of Dynare include the file as part of the 
        % stochastic_solver subdirectory. The function is thus no-longer loaded 
        % as soon as Dynare itself is added to the path. We must add 
        % stochastic_solver manually to the path. 

        pathStochasticSolver = fullfile(pathDynare, 'stochastic_solver');

        % Check if the folder exists
        if exist(pathStochasticSolver, 'dir')
            addpath(pathStochasticSolver);
            disp('[INFO]: Dynare exists and is ready.');
        else
            error("Could not find 'kalman_transition_matrix' in Dynare functions.");
        end
    else
        error("Dynare does not exist in path. Please add Dynare to your path first.");
    end
end
