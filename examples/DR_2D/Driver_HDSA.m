%% Set up
clear
close all
clc
addpath(genpath('../../src'))
load Optimization_Results.mat

surpress_figures = false;

h = .1;
mesh = PDE_Meshing(h);
x = mesh.x;
y = mesh.y;
M = mesh.M;
diff_react_lofi = Diff_React_Lofi(mesh);
diff_react_hifi = Diff_React_Hifi(mesh);
m = length(x);

reg_coeff = 1.e-4;
obj = Diff_React_Objective(diff_react_lofi,reg_coeff);
con = Diff_React_Constraint(diff_react_lofi);
opt = Reduced_Space_Optimization(obj,con);

alpha_u = (2/3)^2;
alpha_z = (1/50000)^2; %(1/200000)^2;
md_interface = Diff_React_HDSA(opt,alpha_u,alpha_z);

%%
num_prior_samples = 100;
md_prior_sampling = HDSA_MD_Prior_Sampling(md_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~surpress_figures
    for k = 1:10
        name = ['Sample ',num2str(k)];
        mesh.Plot_Field(delta_samples(:,k),name)
    end
end
%%
z = zeros(m,3);
z(:,1) = .5*Z(:,1)+.5*250*(1+x).*(1-x).*(1+y).*(1-y);
z(:,2) = Z(:,1).*(1+.02*cos(20*pi*x).*cos(20*pi*y));
z(:,3) = 300*exp(-50*(x.^2+y.^2));
if ~surpress_figures
    for k = 1:3
        name = ['z_',num2str(k)];
        mesh.Plot_Field(z(:,k),name)
    end
end

delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z,num_prior_samples);
if ~surpress_figures
    for j = 1:2
        for k = 1:3
            name = ['Prior discrepancy sample evaluated at z_',num2str(k)];
            mesh.Plot_Field(delta_prior_samples{j}(:,k),name)
        end
    end
end

%%
md_update = HDSA_MD_Update(md_interface);

num_evals = 15;
oversampling = 20;
md_update.Compute_Hessian_GEVP(num_evals,oversampling);

alpha_d = 1.e-2;
num_post_samples = 100;
md_update.Compute_Posterior_Data(alpha_d,num_post_samples);
Z_test = zeros(m,3);
Z_test(:,1:2) = Z;
Z_test(:,3) = z(:,3);
[delta_mean,delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

if ~surpress_figures
    name = 'Y_1';
    mesh.Plot_Field(D(:,1),name)
    
    name = 'Y_1 discrepancy mean';
    mesh.Plot_Field(delta_mean{1},name)
    
    for k = 1:5
        name = 'Y_1 discrepancy sample';
        mesh.Plot_Field(delta_samples{1}(:,k),name)
    end
    
    diff = D(:,1)-delta_mean{1};
    normalize = sqrt(D(:,1)'*md_interface.Apply_M_u(D(:,1)));
    mean_diff = sqrt(diff'*md_interface.Apply_M_u(diff))/normalize;
    sample_diff = zeros(num_post_samples,1);
    for k = 1:num_post_samples
        diff = delta_mean{1} - delta_samples{1}(:,k);
        sample_diff(k) = sqrt(diff'*md_interface.Apply_M_u(diff))/normalize;
    end
    figure,
    hold on
    histogram(sample_diff)
    title(['Mean discrepancy error = ',num2str(mean_diff)])
    
    for k = 1:5
        name = 'Discrepancy at Z_3 sample';
        mesh.Plot_Field(delta_samples{3}(:,k),name)
    end
    
    normalize = sqrt(delta_mean{3}'*md_interface.Apply_M_u(delta_mean{3}));
    sample_diff = zeros(num_post_samples,1);
    for k = 1:num_post_samples
        diff = delta_mean{3} - delta_samples{3}(:,k);
        sample_diff(k) = sqrt(diff'*md_interface.Apply_M_u(diff))/normalize;
    end
    figure,
    hold on
    histogram(sample_diff)
    
end

%%
num_evals = 20;
oversampling = 10;
md_update.Compute_Hessian_GEVP(num_evals,oversampling);

[z_update_mean,z_update_samples] = md_update.Posterior_Update_Samples();

if ~surpress_figures
    name = 'Low-fidelity control';
    mesh.Plot_Field(z_lofi,name)
    
    name = 'Updated control mean';
    mesh.Plot_Field(z_update_mean,name)
    
    for k = 1:10
        name = 'Updated control sample';
        mesh.Plot_Field(z_update_mean,name)
    end
end

%%
u_lofi = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(z_lofi));
u_update_mean = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(z_update_mean));

val_lofi = obj.J(u_lofi,z_lofi);
val_update = obj.J(u_update_mean,z_update_mean);

disp(['Objective at low-fidelity solution = ',num2str(val_lofi)])
disp(['Objective at mean update solution = ',num2str(val_update)])