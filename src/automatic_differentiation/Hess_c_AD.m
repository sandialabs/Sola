function [c] = Hess_c_AD(this, uz, lambda)
    u = uz(1:this.n_u);
    z = uz((this.n_u + 1):end);
    c = lambda' * this.c_AD(u, z);
end
