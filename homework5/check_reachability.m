function [N, XN, K] = check_reachability(model, x0, xf_constr, umin, umax, flag, x_range)
% This function computes the reachable set of a control system and checks
% if a given initial state is contained in it.
%
% Inputs:
% - model: a control system model
% - x0: the initial state
% - xf_constr: constraints on the final state
% - umin: the lower bound on the control input
% - umax: the upper bound on the control input
% - flag: a flag indicating whether there is a limitation for the set of x
%   (1 for yes, 0 for no)
% - x_range: the set of x values that satisfy the limitation (empty for no
%   limitation)
%
% Outputs:
% - N: the number of iterations used to compute the reachable set
% - XN: the final reachable set
% - K: the sequence of reachable sets computed during the iteration


    % Initialize variables
    N = 0;
    XN = xf_constr;
    K(1) = xf_constr;

    % Compute the reachable set and check if the initial state is contained
    % in it
    while N < 100
        % Update the iteration count
        N = N + 1;

        % Compute the reachable set for the current input constraints
        XN = model.reachableSet('X', XN, 'U', Polyhedron('lb', umin, 'ub', umax), 'direction', 'backward');

        % If there is limitation for the set of x
        if flag ~= 0
            XN = XN.intersect(x_range);
        end

        % Add the current reachable set to the invariant set
        K(N+1) = XN;

        % Check if the initial state is contained in the reachable set
        if XN.contains(x0)
            break;
        end
    end

    % If the initial state is not contained in the reachable set after
    % the maximum number of iterations, display a warning
    if N == 100
        warning('The initial state may not be reachable within the maximum number of iterations.');
    end
% Remove redundant inequalities from the computed sets
for i = 1:length(K)
    K(i) = K(i).minHRep();
end
XN = XN.minHRep();

end