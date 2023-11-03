function [c] = AD_u_constraint_residual(this, u, z)
    c = this.constraint_residual(u, z);
end
