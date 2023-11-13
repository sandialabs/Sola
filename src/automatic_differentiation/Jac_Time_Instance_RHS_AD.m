function [f] = Jac_Time_Instance_RHS_AD(this, yz, t)
    y = yz(1:this.n_y);
    z = yz((this.n_y + 1):end);
    f = this.Time_Instance_RHS_AD(y, z, t);
end
