function reachable_set = pre_operation(sys_model, initial_set, horizon_length)
    % Define state and input constraints
    state_constr = Polyhedron('lb', sys_model.x.min, 'ub', sys_model.x.max);
    input_constr = Polyhedron('lb', sys_model.u.min, 'ub', sys_model.u.max);

    % Initialize reachable set to the initial set
    reachable_set = initial_set;

    % Compute reachable set for each time step
    for i = 1:horizon_length
        % Compute reachable set for one time step
        one_step_reachable_set = sys_model.reachableSet('X', reachable_set, ...
            'U', input_constr, 'N', 1, 'direction', 'backward');

        % Intersect reachable set with state constraint set
        reachable_set = state_constr.intersect(one_step_reachable_set);
    end
end
