clear
close all
clc

load Truth_Results.mat
load Std_OED_Results.mat
load Seq_OED_Results.mat

% Set Hi-Fi and Lo-Fi Objectives and Constraints
m = 200;
diff_coeff = 1;
vel_coeff = 2;
robin_coeff = 2;
reg_coeff = 10;

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt = Reduced_Space_Optimization(obj, con);
x = con.x;
M = con.M;

A = con.diff_coeff * con.S + con.robin_coeff * con.robin_bc;
AinvM = (1.e2) * linsolve(A,M);
delta_int = @(theta) theta(1:m);
delta_jac = @(theta)  reshape(theta((m+1):end),m,m)' * M;
T = opt_prob_interface.sabl_opt.obj.T;
reg_mat = opt_prob_interface.sabl_opt.obj.reg_mat;

opt_sol_map = @(theta) linsolve( (AinvM + delta_jac(theta))'*M*(AinvM + delta_jac(theta)) + reg_coeff * reg_mat , (AinvM + delta_jac(theta))'*M*(T-delta_int(theta)));

u_lofi = opt_prob_interface.State_Solve(z_lofi);
Im = eye(m);
B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * M)] * x, u_lofi, z_lofi), z_lofi);
B2 = @(x) [zeros(m, m) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', M)] * x;
B = @(x) B1(x) + B2(x);
PHinvB = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B(x));
opt_sol_lin_map = @(theta) z_lofi - PHinvB(theta);

% Compare [Im kron(Im, z_lofi' * M)] * seq_oed_mean_theta(:,k) with 
% (1 / this.md_post_sampling.post_data.alpha_d) * u in Line 66 of
% MD_Update.m

figure,
plot(x,seq_oed_mean_z(:,end),x,opt_sol_lin_map(seq_oed_mean_theta(:,end)))

seq_oed_hifi_obj = zeros(6,1);
for k = 1:6
    seq_oed_hifi_obj(k) = opt_hifi.Jhat(opt_sol_lin_map(seq_oed_mean_theta(:,k)));
end

% L = 100;
% t = linspace(0,1,L)';
% opt_sol = zeros(m,L);
% opt_sol_mag = zeros(L,1);
% for k = 1:L
%     thetak = t(k) * best_theta;
%     opt_sol(:,k) = opt_sol_map(thetak);
%     opt_sol_mag(k) = sqrt(opt_sol(:,k)'*M*opt_sol(:,k));
% end

% L = 100;
% p = length(best_theta);
% theta_pert = .1 * norm(best_theta) * randn(p,L)/m;
% opt_sol = zeros(m,L);
% opt_sol_lin = zeros(m,L);
% opt_sol_mag = zeros(L,1);
% opt_sol_lin_mag = zeros(L,1);
% opt_sol_obj = zeros(L,1);
% opt_sol_lin_obj = zeros(L,1);
% for k = 1:L
%     thetak = best_theta + theta_pert(:,k);
%     opt_sol(:,k) = opt_sol_map(thetak);
%     opt_sol_mag(k) = sqrt(opt_sol(:,k)'*M*opt_sol(:,k));
%     opt_sol_obj(k) = opt_hifi.Jhat(opt_sol(:,k));
%     opt_sol_lin(:,k) = opt_sol_lin_map(thetak);
%     opt_sol_lin_mag(k) = sqrt(opt_sol_lin(:,k)'*M*opt_sol_lin(:,k));
%     opt_sol_lin_obj(k) = opt_hifi.Jhat(opt_sol_lin(:,k));
% end
% 
% seq_oed_hifi_obj = opt_hifi.Jhat(opt_sol_lin_map(seq_oed_mean_theta(:,6)));
% std_oed_hifi_obj = opt_hifi.Jhat(opt_sol_lin_map(std_oed_mean_theta(:,6)));
% best_theta_hifi_obj = opt_hifi.Jhat(opt_sol_lin_map(best_theta));