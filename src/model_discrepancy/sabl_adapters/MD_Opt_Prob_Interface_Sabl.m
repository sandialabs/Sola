classdef MD_Opt_Prob_Interface_Sabl < MD_Opt_Prob_Interface

    properties
        sabl_opt
        z_current
        u_current
        hessian_data
        m
    end

    methods

        function this = MD_Opt_Prob_Interface_Sabl(sabl_opt, md_interface_data)
            this@MD_Opt_Prob_Interface();
            this.sabl_opt = sabl_opt;
            this.z_current = md_interface_data.z_init;
            [~, ~, this.hessian_data] = this.sabl_opt.Jhat(this.z_current);
            this.m = (length(this.hessian_data) - length(this.z_current)) / 2;
            this.u_current = this.hessian_data(1:this.m);
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sabl_opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            tmp = this.sabl_opt.con.c_u_Transpose_Inverse_Apply(u_in, this.u_current, z);
            z_out = -this.sabl_opt.con.c_z_Transpose_Apply(tmp, this.u_current, z);
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sabl_opt.Jhat(this.z_current);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            z_out = this.sabl_opt.Jhat_hessVec(this.hessian_data, z_in);
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            [~, grad_u, ~] = this.sabl_opt.obj.J(u, z);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.sabl_opt.obj.J_uu_Apply(u_in, u, z);
        end

        function [u] = State_Solve(this, z)
            if norm(z - this.z_current) == 0
                u = this.u_current;
            else
                u = this.sabl_opt.con.State_Solve(z);
                this.u_current = u;
                this.z_current = z;
            end
        end

        function [u_out] = Apply_Solution_Operator_z_Jacobian(this, z_in, z)
            u = this.State_Solve(z);
            tmp = this.sabl_opt.con.c_z_Apply(z_in, u, z);
            u_out = -this.sabl_opt.con.c_u_Inverse_Apply(tmp, u, z);
        end

    end

end
