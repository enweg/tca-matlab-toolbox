% References: 
%   - Morris, S. D. (2016). VARMA representation of DSGE models. 
%     Economics Letters, 138, 30–33. https://doi.org/10.1016/j.econlet.2015.11.027
function [A, B, C, D]=getABCD(M_, oo_, options_)

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
