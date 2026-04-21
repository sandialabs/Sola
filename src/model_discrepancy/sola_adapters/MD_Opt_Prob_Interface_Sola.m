%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Opt_Prob_Interface_Sola < MD_Opt_Prob_Interface

    properties
        sola_opt
        z_current
        u_current
        hessian_data
        m
    end

    %% Implementation of base class virtual functions
    methods

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sola_opt.Jhat(z);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            tmp = this.sola_opt.con.c_u_Transpose_Inverse_Apply(u_in, this.u_current, z);
            z_out = -this.sola_opt.con.c_z_Transpose_Apply(tmp, this.u_current, z);
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sola_opt.Jhat(z);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            z_out = this.sola_opt.Jhat_hessVec(this.hessian_data, z_in);
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            [~, grad_u, ~] = this.sola_opt.obj.J(u, z);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.sola_opt.obj.J_uu_Apply(u_in, u, z);
        end

        function [u_out] = Apply_Solution_Operator_z_Jacobian(this, z_in, z)
            if norm(z - this.z_current) ~= 0
                [~, ~, this.hessian_data] = this.sola_opt.Jhat(z);
                this.z_current = z;
                this.u_current = this.hessian_data(1:this.m);
            end
            tmp = this.sola_opt.con.c_z_Apply(z_in, this.u_current, z);
            u_out = -this.sola_opt.con.c_u_Inverse_Apply(tmp, this.u_current, z);
        end

        function [val, grad_u, grad_z] = Objective_Function(this, u, z)
            [val, grad_u, grad_z] = this.sola_opt.obj.J(u, z);
        end

    end

    %% Constructor and helper function
    methods

        function this = MD_Opt_Prob_Interface_Sola(sola_opt, data_interface)
            arguments
                sola_opt Reduced_Space_Optimization
                data_interface MD_Data_Interface
            end
            this@MD_Opt_Prob_Interface();
            this.sola_opt = sola_opt;
            this.z_current = data_interface.Load_Optimal_z();
            [~, ~, this.hessian_data] = this.sola_opt.Jhat(this.z_current);
            this.m = (length(this.hessian_data) - length(this.z_current)) / 2;
            this.u_current = this.hessian_data(1:this.m);
        end

        function [u] = State_Solve(this, z)
            if norm(z - this.z_current) == 0
                u = this.u_current;
            else
                u = this.sola_opt.con.State_Solve(z);
            end
        end

    end

end
