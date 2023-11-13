function [f] = Hess_Time_Instance_RHS_AD(this, yz, t, lambda)
    y = yz(1:this.n_y);
    z = yz((this.n_y + 1):end);
    f = lambda' * this.Time_Instance_RHS_AD(y, z, t);
end
