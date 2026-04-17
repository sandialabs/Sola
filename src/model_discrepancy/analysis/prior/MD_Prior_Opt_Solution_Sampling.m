classdef MD_Prior_Opt_Solution_Sampling < handle

    properties
        md_data_interface
        md_u_prior_interface
        md_z_prior_interface
        md_opt_prob_interface
        md_hessian_analysis
        u_opt
        z_opt
        scaling
    end

    methods

        function this = MD_Prior_Opt_Solution_Sampling(md_data_interface, md_u_prior_interface, md_z_prior_interface, md_opt_prob_interface, md_hessian_analysis)
            arguments
                md_data_interface MD_Data_Interface
                md_u_prior_interface MD_u_Prior_Interface
                md_z_prior_interface MD_z_Prior_Interface
                md_opt_prob_interface MD_Opt_Prob_Interface
                md_hessian_analysis MD_Hessian_Analysis
            end
            this.md_data_interface = md_data_interface;
            this.md_u_prior_interface = md_u_prior_interface;
            this.md_z_prior_interface = md_z_prior_interface;
            this.md_opt_prob_interface = md_opt_prob_interface;
            this.md_hessian_analysis = md_hessian_analysis;
            this.u_opt = md_data_interface.Load_Optimal_u();
            this.z_opt = md_data_interface.Load_Optimal_z();
            J_grad = this.md_opt_prob_interface.Misfit_Gradient(this.u_opt, this.z_opt);
            this.scaling = sqrt(J_grad' * this.md_u_prior_interface.Apply_W_u_Inverse(J_grad));
        end

        function [opt_sol_samps] = Generate_Prior_Opt_Solution_Samples(this, num_samps)
            u_samps = this.md_u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
            z_samps = this.md_z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samps);

            u_tmp = this.md_opt_prob_interface.Apply_Misfit_Hessian(u_samps, this.u_opt, this.z_opt);
            z_1 = this.md_opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp, this.z_opt);

            z_2 = this.scaling * this.md_z_prior_interface.Apply_M_z(z_samps);

            opt_sol_samps = this.z_opt - this.md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(z_1 + z_2);
        end

    end
end
