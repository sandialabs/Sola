function [sensors, val] = Lazy_Stochastic_CB_Greedy_Solve(f, costs, num_sensors, budget, eps)
    %% Initialize list of weighted gains
    cb_gains = dictionary();
    for v = 1:num_sensors
        if costs(v) <= budget
            cb_gains(v) = struct('cb_gain', inf, 'time', 0);
        end
    end
    sensors = [];
    val = 0;
    cost = 0;
    num_samples = ceil(sum(costs(1:num_sensors)) / budget * log(1 / eps));
    t = 1;
    %% Use cost benefit greedy to add sensors until we can't
    while cb_gains.numEntries > 0
        % sample a subset of sensors
        T = min(num_samples, cb_gains.numEntries);
        indices = datasample(cb_gains.keys, T, 'Replace', false);

        % Create sorted list of weighted gains
        sorted_cb_gains = RedBlackTree();
        for idx = 1:numel(indices)
            v = indices(idx);
            data = cb_gains(v);
            sorted_cb_gains.Insert(data.cb_gain, struct('sensor', v, 'time', data.time));
        end

        % Choose the best sensor
        while true
            [cb_gain, data] = sorted_cb_gains.PopMax();
            v = data.sensor;
            if data.time == t
                sensors(t) = v;
                val = val + cb_gain * costs(v);
                t = t + 1;
                cost = cost + costs(v);
                cb_gains = cb_gains.remove(v);
                break
            else
                cb_gain = (f([sensors, v]) - val) / costs(v);
                sorted_cb_gains.Insert(cb_gain, struct('sensor', v, 'time', t));
                cb_gains(v) = struct('cb_gain', cb_gain, 'time', t);
            end
        end

        % Remove expensive sensors
        remaining_sensors = cb_gains.keys;
        for i = 1:numel(remaining_sensors)
            v = remaining_sensors(i);
            if cost + costs(v) > budget
                cb_gains = cb_gains.remove(v);
            end
        end
    end
end
