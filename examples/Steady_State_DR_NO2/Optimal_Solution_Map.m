function [u, z, D] = Optimal_Solution_Map(T, z0)
    m = length(T);
    diff_coeff = 1;
    react_coeff = -1;
    reg_coeff = 1.e-4;
    obj = Diff_React_Objective(m, reg_coeff);
    obj.T = T;
    con = Diff_React_Constraint(m, diff_coeff, react_coeff);

    opt = Reduced_Space_Optimization(obj, con);

    [u, z] = opt.Optimize(z0);

    tmp = con.c_u_Transpose_Inverse_Apply(con.M, u, z);
    B = con.c_z_Transpose_Apply(tmp, u, z);

    [val, grad, hessian_data] = opt.Jhat(z);
    H = opt.Jhat_hessVec(hessian_data, eye(m));
    D = -linsolve(H, B);

end
