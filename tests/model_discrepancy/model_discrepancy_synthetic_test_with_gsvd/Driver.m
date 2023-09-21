%%
clear
close all
clc
addpath(genpath('../../src'))
rng(121234)

suppress_figures = true;

m = 51;
x = linspace(0,1,m)';

alpha_u = 1/(2^2);
alpha_z = 1/(100^2);
md_interface = HDSA_MD_Interface_synthetic_test_with_gsvd(m,alpha_u,alpha_z);
num_prior_samples = 100;
md_prior_sampling = HDSA_MD_Prior_Sampling(md_interface);

%%
delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~suppress_figures
    figure,
    plot(x,delta_samples,'LineWidth',3,'color',[.9,.9,.9])
end

z = zeros(m,3);
z(:,1) = x;
z(:,2) = x.^2+1;
z(:,3) = sin(2*pi*x);
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
alpha_d = 1.e-5;
num_post_samples = 100;
md_update.Compute_Posterior_Data(alpha_d,num_post_samples);
Z_test = randn(m,3);
Z_test(:,1:2) = md_update.post_data.Z;
Z_test(:,3) = 1.5*ones(m,1);
[delta_mean,delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures
    figure,
    hold on
    plot(x,md_update.post_data.Y(:,1),'color','black','LineWidth',3)
    plot(x,delta_mean{1},'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,delta_samples{1}(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,md_update.post_data.Y(:,1),'color','black','LineWidth',3)
    plot(x,delta_mean{1},'--','color','red','LineWidth',3)
    
    figure,
    hold on
    plot(x,md_update.post_data.Y(:,2),'color','black','LineWidth',3)
    plot(x,delta_mean{2},'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,delta_samples{2}(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,md_update.post_data.Y(:,2),'color','black','LineWidth',3)
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
    plot(x,(1+x)/(1.2^(1/3)),'color','black','LineWidth',3)
    plot(x,1+x,'color','cyan','LineWidth',3)
    plot(x,z_update_mean,'--','color','red','LineWidth',3)
    for k = 1:num_post_samples
        plot(x,z_update_samples(:,k),'color',[.9,.9,.9],'LineWidth',3)
    end
    plot(x,(1+x)/(1.2^(1/3)),'color','black','LineWidth',3)
    plot(x,1+x,'color','cyan','LineWidth',3)
    plot(x,z_update_mean,'--','color','red','LineWidth',3)
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref-z_update_mean)/norm(z_update_mean),norm(z_update_samples-z_samples_ref)/norm(z_update_samples));
if ref_diff>1.e-14
    disp('model_discrepancy_sythetic_test difference:')
    disp(ref_diff)
end