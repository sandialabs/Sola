classdef MD_OUU_Opt_Prob_Interface_Sabl < MD_OUU_Opt_Prob_Interface

    properties
        sabl_opt
        rs_opt_per_sample
        z_current
        u_current
        hessian_data
        m
        n_r
    end

    %% Implementation of base class virtual functions
    methods

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(this, u_in, z, s)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sabl_opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m, :);
            end
            if size(u_in, 2) > 1
                disp('This funcion has not been implemented to support block matrix multiplies');
            end
            tmp = this.sabl_opt.cons{s}.c_u_Transpose_Inverse_Apply(u_in, this.u_current(:, s), z);
            z_out = -this.sabl_opt.cons{s}.c_z_Transpose_Apply(tmp, this.u_current(:, s), z);
        end

        function [z_out] = Apply_RS_Hessian_Per_Sample(this, z_in, z, s)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sabl_opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m, :);
            end
            z_out = this.rs_opt_per_sample{s}.Jhat_hessVec(this.hessian_data(:, s), z_in);
        end

        function [grad_u] = Misfit_Gradient_Per_Sample(this, u, z, s)
            grad_u = this.sabl_opt.obj.J(u, z);
        end

        function [u_out] = Apply_Misfit_Hessian_Per_Sample(this, u_in, u, z, s)
            u_out = this.sabl_opt.obj.J_uu_Apply(u_in, u, z);
        end

    end

    %% Constructor and helper function
    methods

        function this = MD_OUU_Opt_Prob_Interface_Sabl(md_data_interface, sabl_opt)
            arguments
                md_data_interface MD_Data_Interface
                sabl_opt Reduced_Space_Optimization_Under_Uncertainty
            end
            this@MD_OUU_Opt_Prob_Interface(md_data_interface);
            this.sabl_opt = sabl_opt;
            this.z_current = md_data_interface.z_opt;
            [~, ~, this.hessian_data] = this.sabl_opt.Jhat(this.z_current);
            this.m = (size(this.hessian_data, 1) - length(this.z_current)) / 2;
            this.u_current = this.hessian_data(1:this.m, :);
            this.n_r = size(md_data_interface.Xi, 2);
            this.rs_opt_per_sample = cell(this.n_r, 1);
            for s = 1:this.n_r
                this.rs_opt_per_sample{s} = Reduced_Space_Optimization(sabl_opt.obj, sabl_opt.cons{s});
            end
        end

    end

end
