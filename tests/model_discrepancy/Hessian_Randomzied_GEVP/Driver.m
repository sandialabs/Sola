clear
close all
addpath(genpath('../../../src'))
load Optimization_Results.mat
rng(1432)

suppress_figures = true;

obj = Diff_React_Objective(m,reg_coeff);
con = Diff_React_Constraint(m,diff_coeff,react_coeff);
opt_lofi = Reduced_Space_Optimization(obj,con);
x = con.x;

alpha_u = (1.75/1)^2;
alpha_z = (1/200000)^2;
md_interface = Diff_React_HDSA(opt_lofi,alpha_u,alpha_z);

md_update = HDSA_MD_Update(md_interface);

num_evals = 20;
oversampling = 20;
md_update.Compute_Hessian_GEVP(num_evals,oversampling);

H = md_interface.Apply_RS_Hessian(eye(m),z_lofi);
W = md_update.gevp.Apply_Weighting_Operator(eye(m));
[V,D] = eig(H,W,'vector');
for k = 1:m
    V(:,k) = V(:,k)/sqrt(V(:,k)'*W*V(:,k));
end

if ~suppress_figures
    k = num_evals;
    vec = sign(V(1,k))*sign(md_update.gevp.evecs(1,k))*V(:,k);
    figure,
    plot(x,vec,x,md_update.gevp.evecs(:,k),'--','LineWidth',3)
end

if ~suppress_figures
    vec = md_interface.Apply_W_z_Inverse_Factor(randn(m,1));
    tmp1 = md_update.gevp.Apply_Projected_RS_Hessian_Inverse(vec);
    tmp2 = V(:,1:num_evals)*V(:,1:num_evals)'*W*linsolve(H,vec);
    figure,
    plot(x,tmp1,x,tmp2,'--','LineWidth',3)
end

evecs = load('reference_solution.mat','evecs').evecs;
evals = load('reference_solution.mat','evals').evals;
ref_diff = max(norm(evals-md_update.gevp.evals),norm(evecs-md_update.gevp.evecs));

if ref_diff>1.e-14
    disp('Hessian_Randomized_GEVP difference:')
    disp(ref_diff)
end
