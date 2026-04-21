%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val] = State_Norm(u, con)
    u_reshape = reshape(u, con.n_y, con.n_t);
    val = diag(u_reshape' * con.M * u_reshape);
    w = ones(con.n_t, 1);
    w(1) = 1 / 2;
    w(end) = 1 / 2;
    w = con.T * w / sum(w);
    val = sqrt(w' * val);
end
