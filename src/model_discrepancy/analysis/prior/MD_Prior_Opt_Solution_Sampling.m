%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Prior_Opt_Solution_Sampling < handle

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        opt_prob_interface
        hessian_analysis
        u_opt
        z_opt
        scaling
    end

    methods

        function this = MD_Prior_Opt_Solution_Sampling(data_interface, u_prior_interface, z_prior_interface, opt_prob_interface, hessian_analysis)
            arguments
                data_interface MD_Data_Interface
                u_prior_interface MD_u_Prior_Interface
                z_prior_interface MD_z_Prior_Interface
                opt_prob_interface MD_Opt_Prob_Interface
                hessian_analysis MD_Hessian_Analysis
            end
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.opt_prob_interface = opt_prob_interface;
            this.hessian_analysis = hessian_analysis;
            this.u_opt = data_interface.Load_Optimal_u();
            this.z_opt = data_interface.Load_Optimal_z();
            J_grad = this.opt_prob_interface.Misfit_Gradient(this.u_opt, this.z_opt);
            this.scaling = sqrt(J_grad' * this.u_prior_interface.Apply_W_u_Inverse(J_grad));
        end

        function [opt_sol_samps] = Generate_Prior_Opt_Solution_Samples(this, num_samps)
            u_samps = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
            z_samps = this.z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samps);

            u_tmp = this.opt_prob_interface.Apply_Misfit_Hessian(u_samps, this.u_opt, this.z_opt);
            z_1 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp, this.z_opt);

            z_2 = this.scaling * this.z_prior_interface.Apply_M_z(z_samps);

            opt_sol_samps = this.z_opt - this.hessian_analysis.Apply_Projected_RS_Hessian_Inverse(z_1 + z_2);
        end

    end
end
