function [f] = Hess_Time_Instance_RHS_AD(this, yz, t, lambda)
    y = yz(1:this.m);
    z = yz((this.m + 1):end);
    f = lambda' * this.Time_Instance_RHS_AD(y, z, t);
end
