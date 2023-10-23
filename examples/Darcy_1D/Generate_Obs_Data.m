function [] = Generate_Obs_Data(con,forcing_fun)
    z_true = 2 + cos(2*pi*con.x);
    forcing = forcing_fun(con.x);
    con.forcing = forcing;
    u_true = con.State_Solve(z_true);
    u_data = u_true.*(1+0.02*randn(size(u_true)));
    save('Obs_Data.mat','u_data','u_true','z_true','forcing')
end

