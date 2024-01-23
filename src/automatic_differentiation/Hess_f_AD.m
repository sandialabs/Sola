function [f] = Hess_f_AD(this, yz, t, lambda)
    y = yz(1:this.n_y);
    z = yz((this.n_y + 1):end);
    f = lambda' * this.f_AD(y, z, t);
end
