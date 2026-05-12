%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

m = 100;
[M, S, nodes] = Assemble_Mass_and_Stiffness(m);

%%
alpha_u = 3.4;
beta_u = 0.026;

laplacian_op_properties = MD_Laplacian_Like_Operator_Properties();
e = laplacian_op_properties.Get_Rectangular_Domain_Laplacian_Like_Operator_Evals(beta_u, nodes);

h = 1/(m-1);
C = 1/6;
gamma_bound = h^2 * C * (1/beta_u) * sum((e-1)./(e.^2));
gamma_var = alpha_u * laplacian_op_properties.Get_Rectangular_Domain_Squared_Inv_Operator_Trace(beta_u, nodes);
scaling = gamma_var / (gamma_var + gamma_bound);

E = beta_u * S + M;
Einv = E \ eye(m);
Winv = alpha_u * Einv * M * Einv;

M_lumped = diag(M*ones(m,1));
Winv_bar = alpha_u * Einv * M_lumped * Einv;

var_diff = scaling * trace(Winv_bar*M) - trace(Winv*M);

%%
data_interface = MD_Data_Interface_Test(m);
u_hyperparam_interface = MD_u_Hyperparameter_Interface(false);
u_hyperparam_interface.Set_alpha_u(alpha_u)
u_hyperparam_interface.Set_beta_u(beta_u);

u_prior_lumped_mass = MD_Lumped_Mass_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);
u_prior_bilaplacian = MD_Bilaplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

num_samples = 100;
lumped_mass_samples = u_prior_lumped_mass.Sample_with_Covariance_W_u_Inverse(num_samples);
bilaplacian_samples = u_prior_bilaplacian.Sample_with_Covariance_W_u_Inverse(num_samples);