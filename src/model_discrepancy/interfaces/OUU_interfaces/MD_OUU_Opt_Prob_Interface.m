classdef MD_OUU_Opt_Prob_Interface < MD_Opt_Prob_Interface

    properties
        data_interface
    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(this, u_in, z, s)

        [z_out] = Apply_RS_Hessian_Per_Sample(this, z_in, z, s)

        [grad_u] = Misfit_Gradient_Per_Sample(this, u, z, s)

        [u_out] = Apply_Misfit_Hessian_Per_Sample(this, u_in, u, z, s)

    end

    %% Virtual functions for user implementation
    methods

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            n = size(u_in, 2);
            z_out = zeros(length(z), n);
            for k = 1:n
                z_out_k = zeros(length(z), this.data_interface.n_r);
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                for s = 1:this.data_interface.n_r
                    z_out_k(:, s) = this.Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(u_in_k(:, s), z, s);
                end
                z_out(:, k) = sum(z_out_k, 2);
            end
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            n = size(z_in,2);
            z_out = zeros(length(z),n);
            for k = 1:n
                z_out_k = zeros(length(z), this.data_interface.n_r);
                for s = 1:this.data_interface.n_r
                    z_out_k(:, s) = this.Apply_RS_Hessian_Per_Sample(z_in(:,k), z, s);
                end
                z_out(:,k) = mean(z_out_k, 2);
            end
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            u = this.data_interface.Reshape_State_to_Mat(u);
            grad_u = 0 * u;
            for s = 1:this.data_interface.n_r
                grad_u(:, s) = this.Misfit_Gradient_Per_Sample(u(:, s), z, s);
            end
            grad_u = (1/this.data_interface.n_r) * this.data_interface.Reshape_State_to_Vec(grad_u);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = 0 * u_in;
            n = size(u_in, 2);
            u = this.data_interface.Reshape_State_to_Mat(u);
            for k = 1:n
                u_in_k = this.data_interface.Reshape_State_to_Mat(u_in(:, k));
                u_out_tmp = 0 * u_in_k;
                for s = 1:this.data_interface.n_r
                    u_out_tmp(:, s) = this.Apply_Misfit_Hessian_Per_Sample(u_in_k(:, s), u(:, s), z, s);
                end
                u_out(:, k) = (1/this.data_interface.n_r) * this.data_interface.Reshape_State_to_Vec(u_out_tmp);
            end
        end

    end

    %% Constructor
    methods

        function this = MD_OUU_Opt_Prob_Interface(data_interface)
            arguments
                data_interface MD_OUU_Data_Interface
            end
            this.data_interface = data_interface;
        end

    end

end
