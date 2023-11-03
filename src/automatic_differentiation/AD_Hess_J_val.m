function [val] = AD_Hess_J_val(this, uz)
    u = uz(1:this.n_u);
    z = uz((this.n_u + 1):end);
    val = this.J_val(u, z);
end
