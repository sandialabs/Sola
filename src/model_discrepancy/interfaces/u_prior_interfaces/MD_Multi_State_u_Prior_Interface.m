classdef MD_Multi_State_u_Prior_Interface < MD_u_Prior_Interface

    properties
        u_prior_interface_cell
        u_hyperparam_interface
        n_c
        I
        total_dofs
    end

    %% Constructor
    methods

        function this = MD_Multi_State_u_Prior_Interface(data_interface, u_prior_interface_cell, u_hyperparam_interface_cell)
            this.u_prior_interface_cell = u_prior_interface_cell;
            this.u_hyperparam_interface = MD_Multi_State_u_Hyperparameter_Interface(u_hyperparam_interface_cell);
            this.n_c = length(u_prior_interface_cell);
            this.I = cell(this.n_c, 1);
            this.total_dofs = 0;
            for i = 1:this.n_c
                this.I{i} = data_interface.Separate_State_Components(i);
                this.total_dofs = this.total_dofs + length(this.I{i});
            end
        end

    end

    %% Implementation of base class functions
    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = 0 * u_in;
            for i = 1:this.n_c
                if size(u_in,1) == this.total_dofs
                    J = this.I{i};
                else
                    J = this.I{i}(1:this.u_prior_interface_cell{i}.transient_prior_cov.n_y);
                end
                u_out(J, :) = this.u_prior_interface_cell{i}.Apply_M_u(u_in(J, :));
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0 * u_in;
            for i = 1:this.n_c
                u_out(this.I{i}, :) = this.u_prior_interface_cell{i}.Apply_W_u_Plus_scalar_M_u_Inverse(u_in(this.I{i}, :), scalar);
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = 0 * u_in;
            for i = 1:this.n_c
                u_out(this.I{i}, :) = this.u_prior_interface_cell{i}.Apply_W_u_Inverse(u_in(this.I{i}, :));
            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            u_out = zeros(this.total_dofs, num_samples);
            for i = 1:this.n_c
                u_out(this.I{i}, :) = this.u_prior_interface_cell{i}.Sample_with_Covariance_W_u_Inverse(num_samples);
            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            u_out = zeros(this.total_dofs, num_samples);
            for i = 1:this.n_c
                u_out(this.I{i}, :) = this.u_prior_interface_cell{i}.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samples, scalar);
            end
        end

    end

end
