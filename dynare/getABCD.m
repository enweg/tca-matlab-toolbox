function [A, B, C, D]=getABCD(M_, oo_, options_)
    % `getABCD` Obtain the ABCD state-space representation of a DSGE model.
    %
    %   `[A, B, C, D] = getABCD(M_, oo_, options_)` computes the state-space
    %   representation 
    %   $$
    %   \begin{split}
    %   x_t &= Ax_{t-1} + B\varepsilon_t \\ 
    %   y_t &= Cx_{t-1} + D\varepsilon_t
    %   \end{split}
    %   of a DSGE model estimated using Dynare. Only the minimal state 
    %   representation is returned. 
    %
    %   ## Arguments
    %   - `M_` (struct): Returned by Dynare. 
    %   - `oo_` (struct): Returned by Dynare. 
    %   - `options_` (struct): Returned by Dynare. 
    %
    %   ## Returns
    %   - `A` (matrix): State transition matrix. See above equation.
    %   - `B` (matrix): Control input matrix capturing exogenous shocks. See above equation.
    %   - `C` (matrix): Observation matrix mapping state variables to observed variables. See above equation.
    %   - `D` (matrix): Observation noise matrix. See above equation.
    %
    %   ## Notes
    %   - Requires MATLAB's Control Toolbox. 
    %

    if ~isfield(options_,'varobs_id')
        warning('getABCD: No observables have been defined using the varobs-command.')
        return;    
    end

    % Dynare re-orders variables into the order static, backward, mixed, forward. 
    % The state variables are the backward and mixed variables. 
    % Thus, in the DR (internal order) ordering, the state variables are given by 
    % the following indices, where nspred is the number of state variables.
    ipred = M_.nstatic+(1:M_.nspred)';
    % options_.varobs_id is in declaration order. Need to change this to internal DR
    % order for ABCD matrices. 
    obs_var=oo_.dr.inv_order_var(options_.varobs_id);

    % get state transition matrices
    [A,B] = kalman_transition_matrix(oo_.dr,ipred,1:M_.nspred,M_.exo_nbr);
    % get observation equation matrices
    [C,D] = kalman_transition_matrix(oo_.dr,obs_var,1:M_.nspred,M_.exo_nbr);

    % We need the minimum state representation
    if user_has_matlab_license('control_toolbox')
        [A,B,C,D]=minreal(A,B,C,D); % Matlab control toolbox; TODO: find alternative
    else
        error('Control Toolbox is missing')
    end
end
