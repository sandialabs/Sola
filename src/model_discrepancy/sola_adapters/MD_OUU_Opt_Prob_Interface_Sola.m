%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_OUU_Opt_Prob_Interface_Sola < MD_OUU_Opt_Prob_Interface

    properties
        sola_opt
        rs_opt_per_sample
        n_r
    end

    %% Implementation of base class virtual functions
    methods

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(this, u_in, z, s)

            if size(u_in, 2) > 1
                disp('This funcion has not been implemented to support block matrix multiplies');
            end

            u = this.State_Solve_Per_Sample(z, s);
            tmp = this.sola_opt.cons{s}.c_u_Transpose_Inverse_Apply(u_in, u, z);
            z_out = -this.sola_opt.cons{s}.c_z_Transpose_Apply(tmp, u, z);
        end

        function [z_out] = Apply_RS_Hessian_Per_Sample(this, z_in, z, s)
            [~, ~, hessian_data] = this.sola_opt.Jhat(z);
            z_out = this.rs_opt_per_sample{s}.Jhat_hessVec(hessian_data(:, s), z_in);
        end

        function [grad_u] = Misfit_Gradient_Per_Sample(this, u, z, s)
            [~, grad_u] = this.sola_opt.obj.J(u, z);
        end

        function [u_out] = Apply_Misfit_Hessian_Per_Sample(this, u_in, u, z, s)
            u_out = this.sola_opt.obj.J_uu_Apply(u_in, u, z);
        end

        function [u] = State_Solve_Per_Sample(this, z, s)
            u = this.sola_opt.cons{s}.State_Solve(z);
        end

    end

    %% Constructor and helper function
    methods

        function this = MD_OUU_Opt_Prob_Interface_Sola(data_interface, sola_opt)
            arguments
                data_interface MD_Data_Interface
                sola_opt Reduced_Space_Optimization_Under_Uncertainty
            end
            this@MD_OUU_Opt_Prob_Interface(data_interface);
            this.sola_opt = sola_opt;
            this.n_r = length(sola_opt.cons);
            this.rs_opt_per_sample = cell(this.n_r, 1);
            for s = 1:this.n_r
                this.rs_opt_per_sample{s} = Reduced_Space_Optimization(sola_opt.obj, sola_opt.cons{s});
            end
        end

    end

end
