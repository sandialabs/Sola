classdef HDSA_Sabl_MD_Interface_Elliptic_Prior < HDSA_MD_Interface_Elliptic_Prior

    properties
        opt
        z_current
        u_current
        hessian_data
        m
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_out] = Apply_E_u_Inverse(this, u_in)

        [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)

        [u_out] = Apply_M_u(this, u_in)

        [u_out] = Apply_M_u_Inverse(this, u_in)

        [z_out] = Apply_E_z_Inverse(this, z_in)

        [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)

        [z_out] = Apply_M_z(this, z_in)

        [u_out] = Apply_E_d(this, u_in)

        [u_out] = Apply_E_d_Transpose(this, u_in)

        [u_opt] = Load_Optimal_u(this)

        [z_opt] = Load_Optimal_z(this)

        [Z] = Load_Z_Data(this)

        [D] = Load_d_Data(this)

    end

    methods

        function this = HDSA_Sabl_MD_Interface_Elliptic_Prior(opt, alpha_u, alpha_z)
            this@HDSA_MD_Interface_Elliptic_Prior(alpha_u, alpha_z);
            this.opt = opt;
            this.z_current = this.Load_Optimal_z();
            [~, ~, this.hessian_data] = this.opt.Jhat(this.z_current);
            this.m = (length(this.hessian_data) - length(this.z_current)) / 2;
            this.u_current = this.hessian_data(1:this.m);
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            tmp = this.opt.con.c_u_Transpose_Inverse_Apply(u_in, this.u_current, z);
            z_out = -this.opt.con.c_z_Transpose_Apply(tmp, this.u_current, z);
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            z_out = this.opt.Jhat_hessVec(this.hessian_data, z_in);
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            [~, grad_u, ~] = this.opt.obj.J(u, z);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.opt.obj.J_uu_Apply(u_in, u, z);
        end

    end

end
