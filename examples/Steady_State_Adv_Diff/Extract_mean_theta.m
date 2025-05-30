function [theta] = Extract_mean_theta(z_prior_interface, post_data)
    m = size(post_data.D, 1);
    n = size(post_data.Z, 1);
    p = m * (n + 1);
    N = size(post_data.D, 2);
    theta = zeros(p, N);
    I1 = 1:m;
    I2 = (m + 1):p;
    for ell = 1:N

        theta(I1, ell) = post_data.a_ell(ell) * post_data.u_ell(:, ell);
        tmp = z_prior_interface.Apply_M_z_Inverse(z_prior_interface.Apply_W_z_Inverse(post_data.Z(:, ell) - post_data.z_opt));
        theta(I2, ell) = kron(post_data.u_ell(:, ell), tmp);

        for i = 1:N
            yi = post_data.Z * post_data.g_vecs(:, i) - sum(post_data.g_vecs(:, i)) * post_data.z_opt;
            si = sum(post_data.g_vecs(:, i)) - yi' * post_data.W_z_inv_z_opt;
            theta(I1, ell) = theta(I1, ell) - post_data.b_i_ell(i, ell) * si * post_data.u_i_ell{i}(:, ell);
            tmp = z_prior_interface.Apply_M_z_Inverse(z_prior_interface.Apply_W_z_Inverse(yi));
            theta(I2, ell) = theta(I2, ell) - post_data.b_i_ell(i, ell) * kron(post_data.u_i_ell{i}(:, ell), tmp);
        end

    end

    theta = (1 / post_data.alpha_d) * sum(theta, 2);
end
