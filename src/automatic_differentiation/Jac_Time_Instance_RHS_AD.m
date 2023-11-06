function [f] = Jac_Time_Instance_RHS_AD(this, yz, t)
    y = yz(1:this.m);
    z = yz((this.m + 1):end);
    f = this.Time_Instance_RHS_AD(y, z, t);
end
