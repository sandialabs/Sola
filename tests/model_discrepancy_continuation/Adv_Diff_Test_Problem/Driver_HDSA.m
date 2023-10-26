%%
clear
close all
clc
addpath(genpath('../../../src'))
rng(1234423)

m = 200;
diff_coeff = 1;
vel_coeff = 1;
robin_coeff = 2; 
reg_coeff = 10;
obj = Adv_Diff_Objective(m,reg_coeff);
con_hifi = Adv_Diff_Constraint(m,diff_coeff,vel_coeff,robin_coeff);
con_lofi = Diff_Constraint(obj,con_hifi);
opt_hifi = Reduced_Space_Optimization(obj,con_hifi);
opt_lofi = Reduced_Space_Optimization(obj,con_lofi);
x = con_hifi.x;

%%
alpha_u = 1/(2^2);
alpha_z = 1/(600^2);
md_interface = HDSA_Sabl_MD_Interface_Elliptic_Prior_PDE_Test_Prob(opt_lofi,alpha_u,alpha_z);
md_continuation_interface = HDSA_Sabl_MD_Continuation_Interface(md_interface,con_lofi);

%%
continuation_step_sweep = 1:5;
p = length(continuation_step_sweep);
u = cell(p,1);
z = cell(p,1);

for k = 1:p
    num_continuation_steps = continuation_step_sweep(k);
    md_update = HDSA_MD_Continuation_Update(md_continuation_interface,num_continuation_steps);
    alpha_d = 1.e-5;
    num_post_samples = 100;
    md_update.Compute_Posterior_Data(alpha_d,num_post_samples);
    [u{k},z{k}] = md_update.Posterior_Update_Mean();
end

%%
z_hifi = load('z_hifi.mat').z_hifi;
figure,
hold on
plot(x,md_update.z_opt,'color','black','LineWidth',3)
plot(x,z_hifi,'color','cyan','LineWidth',3)
lg = cell(p+2,1);
lg{1} = 'Low-fidelity solution';
lg{2} = 'High-fidelity solution';
for k = 1:p
    plot(x,z{k}(:,end),'--','LineWidth',3)
    lg{k+2} = [num2str(k),' step posterior mean update'];
end
legend(lg)

