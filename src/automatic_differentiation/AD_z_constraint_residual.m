function [c] = AD_z_constraint_residual(this, u, z)
    c = this.constraint_residual(u, z);
end
