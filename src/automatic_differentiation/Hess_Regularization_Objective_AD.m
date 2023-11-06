function [val] = Hess_Regularization_Objective_AD(this, z)
    val = this.Regularization_Objective_AD(z);
end
