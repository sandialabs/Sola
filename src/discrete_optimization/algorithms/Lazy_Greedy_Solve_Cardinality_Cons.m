function [sensors, val] = Lazy_Greedy_Solve_Cardinality_Cons(f, num_sensors, budget)
    sensors = [];
    marginal_gains = RedBlackTree();
    val = 0;
    t = 1;
    for v = 1:num_sensors
        marginal_gains.Insert(f(v), struct('sensor', v, 'time', t));
    end
    while t <= budget
        [gain, data] = marginal_gains.PopMax();
        v = data.sensor;
        if data.time == t
            sensors(t) = v;
            val = val + gain;
            t = t + 1;
        else
            gain = f([sensors, v]) - val;
            marginal_gains.Insert(gain, struct('sensor', v, 'time', t));
        end
    end
end