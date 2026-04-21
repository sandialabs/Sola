%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Plot_Control(x, z, con)
    zt = con.B * reshape(z, con.num_space_control_nodes, length(con.z_time_mesh));
    ymin = min(zt(:)) - .05 * abs(min(zt(:)));
    ymax = max(zt(:)) + .05 * abs(max(zt(:)));
    figure;
    hold on;
    for j = 1:length(con.z_time_mesh)
        t = con.z_time_mesh(j);
        plot(x, zt(:, j), 'LineWidth', 3, 'Color', (1 - .8 * t) * ones(3, 1));
        ylim([ymin, ymax]);
        title(['Time = ', num2str(t)]);
        pause(.05);
    end
end
