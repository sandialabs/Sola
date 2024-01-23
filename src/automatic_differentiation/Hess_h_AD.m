function [h] = Hess_h_AD(this, z, lambda)
    h = lambda' * this.h_AD(z);
end
