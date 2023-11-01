function [] = Generate_Obs_Data(con, z_true, forcing)
    con.forcing = forcing;
    u_true = con.State_Solve(z_true);
    u_data = u_true .* (1 + 0.02 * randn(size(u_true)));
    save('Obs_Data.mat', 'u_data', 'u_true', 'z_true', 'forcing');
end
