function [val] = Hess_Time_Instance_Objective_AD(this, y, t)
    val = this.Time_Instance_Objective_AD(y, t);
end
