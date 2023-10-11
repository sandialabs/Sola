%% Set up
clear
close all
clc
addpath(genpath('../../src'))
load Optimization_Results.mat

suppress_figures = false;

obj = Adv_Diff_Objective(m,reg_coeff);
con_hifi = Adv_Diff_Constraint(m,diff_coeff,vel_coeff,robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj,con_hifi);
opt_lofi = Reduced_Space_Optimization(obj,con_lofi);
x = con_lofi.x;

alpha_u = (1/2)^2;
alpha_z = (1/100)^2;
md_interface = Diff_HDSA(opt_lofi,alpha_u,alpha_z);

%%
num_prior_samples = 100;
md_prior_sampling = HDSA_MD_Prior_Sampling(md_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
if ~suppress_figures
    figure,
    plot(x,delta_samples(:,1:10),'LineWidth',3)
    
    figure,
    plot(x,delta_samples,'LineWidth',3,'color',[.9,.9,.9])
end

%%
z = zeros(m,3);
z(:,1) = Z(:,1) + .1*x.*(1-x);
z(:,2) = Z(:,1).*(1+.02*cos(20*pi*x));
z(:,3) = exp(-1*(x-0.5).^2);
if ~suppress_figures
    figure,
    plot(x,z,'LineWidth',3)
end

delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z,num_prior_samples);
if ~suppress_figures
    for k = 1:10
        figure,
        hold on
        plot(x,delta_prior_samples{k},'LineWidth',3)
    end
end

%%
md_update = HDSA_MD_Update(md_interface);
alpha_d = 1.e-2;
num_post_samples = 100;
md_update.Compute_Posterior_Data(alpha_d,num_post_samples);
Z_test = zeros(m,3);
Z_test(:,1:2) = Z;
Z_test(:,3) = z(:,3);
[delta_mean,delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures
    figure,
    hold on
    plot(x,md_update.post_data.D(:,1),'color','black','LineWidth',3)
    plot(x,delta_mean{1},'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,delta_samples{1}(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,md_update.post_data.D(:,1),'color','black','LineWidth',3)
    plot(x,delta_mean{1},'--','color','red','LineWidth',3)
    
    figure,
    hold on
    plot(x,md_update.post_data.D(:,2),'color','black','LineWidth',3)
    plot(x,delta_mean{2},'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,delta_samples{2}(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,md_update.post_data.D(:,2),'color','black','LineWidth',3)
    plot(x,delta_mean{2},'--','color','red','LineWidth',3)
    
    figure,
    hold on
    plot(x,delta_mean{3},'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,delta_samples{3}(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,delta_mean{3},'--','color','red','LineWidth',3)
    
end

%%
[z_update_mean,z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    figure,
    hold on
    plot(x,z_lofi,'color','black','LineWidth',3)
    plot(x,z_hifi,'color','cyan','LineWidth',3)
    plot(x,z_update_mean,'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,z_update_samples(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,md_update.z_opt,'color','black','LineWidth',3)
    plot(x,z_hifi,'color','cyan','LineWidth',3)
    plot(x,z_update_mean,'--','color','red','LineWidth',3)
    legend({'Low-fidelity control','High-fidelity control','Update'})
end