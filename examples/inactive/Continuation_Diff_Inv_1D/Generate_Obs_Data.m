%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Generate_Obs_Data(con)

    z_true = 1 + sin(2 * pi * con.x).^2;
    theta_true = 1 + 0.3 * (1 - con.x).^2;
    u_true = con.Parameterized_State_Solve(z_true, theta_true);
    u_data = u_true .* (1 + 0.02 * randn(length(u_true), 1));

    save('Obs_Data.mat', 'u_data', 'u_true', 'z_true');
end
