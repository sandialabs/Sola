function [D] = Evaluate_Discrepancy(con_hifi, con_lofi, Z)
    % D is not preallocated so perhaps less efficient (generally, N is small though)
    for i = 1:size(Z, 2)
        D(:, i) = con_hifi.State_Solve(Z(:, i)) - con_lofi.State_Solve(Z(:, i));
    end
end
