function [c] = AD_zz_constraint_residual(this, u, z, lambda)
    c = lambda' * this.constraint_residual(u, z);
end
