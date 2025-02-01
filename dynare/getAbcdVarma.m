function [Ax, Ay, B, Cx, Cy, D, vars]=getAbcdVarma(M_, options_, oo_, obsVar)
    if ~isfield(options_,'varobs_id')
        warning('getAbcdVarma: No observables have been defined using the varobs-command.')
        return;    
    end

    %get state indices
    ipred = M_.nstatic+(1:M_.nspred)';
    %get observable position in decision rules
    if nargin == 3
        warning("getAbcdVarma: Using observed variables from .mod file");
        obsVar=oo_.dr.inv_order_var(options_.varobs_id);
    end

    vars = [ipred; obsVar];

    %get state transition matrices
    [A,B] = kalman_transition_matrix(oo_.dr,ipred,1:M_.nspred,M_.exo_nbr);
    %get observation equation matrices
    [C,D] = kalman_transition_matrix(oo_.dr,obsVar,1:M_.nspred,M_.exo_nbr);

    % Stack them into a VAR
    ns = size(A, 1);
    nc = size(C, 1);
    A1 = [A zeros(ns, nc); C zeros(nc, nc)];
    Psi = [B; D];
    % Now sort them into observed and unobserved
    % that is, whenever a variable is observed and a state variable 
    % is will be twice in the current system. We need to remove the second 
    % appearance of it, and move the state equation into the observed block.
    for vo=obsVar'
        if any(ismember(vo, ipred))
            % observed variable is also state variable 
            ixLocation = find(vars==vo);
            ixS = ixLocation(1);
            ixO = ixLocation(2);

            % copy over row and column of A1
            A1(ixO, :) = A1(ixS, :);
            A1(:, ixO) = A1(:, ixS);;;;
            % then remove ixS column and row 
            A1 = A1([1:(ixS-1), (ixS+1):size(A1, 1)], [1:(ixS-1), (ixS+1):size(A1, 1)]);

            % copy over row in Psi
            Psi(ixO, :) = Psi(ixS, :);
            % remove ixS row
            Psi = Psi([1:(ixS-1), (ixS+1):size(Psi, 1)], :);

            % remove entry in vars
            vars = vars([1:(ixS-1), (ixS+1):size(vars, 1)]);
        end
    end


    n = size(A1, 1);
    no = size(obsVar, 1);
    ns = n - no;
    % removing all unnecessary state variables
    keep = repelem(true, size(A1, 1));
    for i=1:ns
        keep(i) = any(abs(A1([1:(i-1), (i+1):size(A1, 1)], i)) > 1e-10);
        if ~keep(i)
            ns = ns - 1;
            warning("getAbcdVarma: Removing state variable " + get_variable_ordering(M_, options_, oo_, vars(i)));
        end
    end
    A1 = A1(keep, keep);
    Psi = Psi(keep, :);
    vars = vars(keep);

    D = Psi((ns+1):end, :);
    if size(D, 1) ~= size(D, 2)
        error("D is not square.");
    end
    Ax = A1(1:ns, 1:ns);
    Ay = A1(1:ns, (ns+1):end);
    Cx = A1((ns+1):end, 1:ns);
    Cy = A1((ns+1):end, (ns+1):end);
    B = Psi(1:ns, :);
end
