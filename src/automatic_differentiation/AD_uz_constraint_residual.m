function [c] = AD_uz_constraint_residual(this, uz, lambda)
    u = uz(1:this.n_u);
    z = uz((this.n_u + 1):end);
    c = lambda' * this.constraint_residual(u, z);
end
