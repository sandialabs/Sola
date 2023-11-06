function [c] = Jac_c_AD(this, uz)
    u = uz(1:this.n_u);
    z = uz((this.n_u + 1):end);
    c = this.c_AD(u, z);
end
