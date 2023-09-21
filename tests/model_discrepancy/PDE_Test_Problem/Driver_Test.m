clear
close all
clc
addpath(genpath('../../src'))

m = 200;
diff_coeff = 1;
vel_coeff = 1/2;
robin_coeff = 2; 
reg_coeff = 10;
obj_hifi = Adv_Diff(m,diff_coeff,vel_coeff,robin_coeff,reg_coeff);
obj_lofi = Diff(obj_hifi);
x = obj_hifi.x;
z0 = rand(m,1);

mms_check = true;
finite_diff_check = true;

if mms_check
    c = 4*robin_coeff/(4+robin_coeff);
    y = 1-c*(x-.5).^2;
    z = (10^-2)*(2*c*diff_coeff - 2*c*vel_coeff*(x-.5));
    u = obj_hifi.State_Solve(z);
    
    figure,
    hold on
    plot(x,u,'LineWidth',3)
    plot(x,y,'--','LineWidth',3)
end

if finite_diff_check
    diffs = obj_hifi.Finite_Difference_Gradient_Check(z0);
    diffs = obj_hifi.Finite_Difference_Hessian_Check(z0);
end

if mms_check
    c = 4*robin_coeff/(4+robin_coeff);
    y = 1-c*(x-.5).^2;
    z = (10^-2)*(2*c*diff_coeff + 0*x);
    u = obj_lofi.State_Solve(z);
    
    figure,
    hold on
    plot(x,u,'LineWidth',3)
    plot(x,y,'--','LineWidth',3)
end

if finite_diff_check
    diffs = obj_lofi.Finite_Difference_Gradient_Check(z0);
    diffs = obj_lofi.Finite_Difference_Hessian_Check(z0);
end
