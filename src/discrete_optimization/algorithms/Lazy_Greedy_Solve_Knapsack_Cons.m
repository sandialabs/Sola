function [sensors, val] = Lazy_Greedy_Solve_Knapsack_Cons(f, costs, num_sensors, budget)
    %% Initialize lists of marginal gains
    marginal_gains = RedBlackTree();
    weighted_gains = RedBlackTree();
    for i = 1:num_sensors
        gain = f(i);
        marginal_gains.Insert(gain, struct('sensor', i, 'time', 1));
        weighted_gains.Insert(gain / costs(i), struct('sensor', i, 'time', 1));
    end
    %% Find greedy solution
    greedy_sensors = [];
    greedy_cost = 0;
    val = 0;
    t = 1;
    while marginal_gains.size > 0
        [gain, data] = marginal_gains.PopMax();
        v = data.sensor;
        if greedy_cost + costs(v) > budget
            continue
        elseif data.time == t
            greedy_sensors(t) = v;
            val = val + gain;
            t = t + 1;
            greedy_cost = greedy_cost + costs(v);
        else
            gain = f([greedy_sensors, v]) - val;
            marginal_gains.Insert(gain, struct('sensor', v, 'time', t));
        end
    end
    %% Find cost-benefit greedy solution
    cb_sensors = [];
    cb_cost = 0;
    cb_val = 0;
    t = 1;
    while weighted_gains.size > 0
        [cb_gain, data] = weighted_gains.PopMax();
        v = data.sensor;
        if cb_cost + costs(v) > budget
            continue
        elseif data.time == t
            cb_sensors(t) = v;
            cb_val = cb_val + cb_gain * costs(v);
            t = t + 1;
            cb_cost = cb_cost + costs(v);
        else
            cb_gain = (f([cb_sensors, v]) - cb_val) / costs(v);
            weighted_gains.Insert(cb_gain, struct('sensor', v, 'time', t));
        end
    end
    %% Pick the better of the two solutions
    if val > cb_val
        sensors = greedy_sensors;
    else
        sensors = cb_sensors;
        val = cb_val;
    end
end