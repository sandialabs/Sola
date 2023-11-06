function [h] = Hess_Initial_Condition_AD(this, z, lambda)
    h = lambda' * this.Initial_Condition_AD(z);
end
