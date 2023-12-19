function [] = Plot_States(u, con, fig_nums)

    n_y = con.n_y / 4;
    n_t = con.n_t;

    u_rs = reshape(u, 4 * n_y, n_t);
    T = con.I_T * u_rs;
    u1 = con.I_u1 * u_rs;
    v1 = con.I_v1 * u_rs;
    v2 = con.I_v2 * u_rs;

    figure(fig_nums(1));
    surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', T);
    xlabel('Space');
    ylabel('Time');
    title('T');
    colorbar();
    set(gca, 'fontsize', 18);

    figure(fig_nums(2));
    surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', u1);
    xlabel('Space');
    ylabel('Time');
    title('u1');
    colorbar();
    set(gca, 'fontsize', 18);

    figure(fig_nums(3));
    surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v1);
    xlabel('Space');
    ylabel('Time');
    title('v1');
    colorbar();
    set(gca, 'fontsize', 18);

    figure(fig_nums(4));
    surf(con.fe.x * ones(1, n_t), ones(n_y, 1) * con.t_mesh', v2);
    xlabel('Space');
    ylabel('Time');
    title('v2');
    colorbar();
    set(gca, 'fontsize', 18);

end
