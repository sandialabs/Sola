classdef MD_Multi_State_u_Prior_Interface < MD_u_Prior_Interface

    properties
        u_prior_interface_cell
        n_c
        I
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = 0 * u_in;
            for i = 1:this.n_c
                u_out(this.I{i},:) = this.u_prior_interface_cell{i}.Apply_M_u(u_in(this.I{i},:));
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0 * u_in;
            for i = 1:this.n_c
                u_out(this.I{i},:) = this.u_prior_interface_cell{i}.Apply_W_u_Plus_scalar_M_u_Inverse(u_in(this.I{i},:),scalar);
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = 0 * u_in;
            for i = 1:this.n_c
                u_out(this.I{i},:) = this.u_prior_interface_cell{i}.Apply_W_u_Inverse(u_in(this.I{i},:));
            end
        end

        % Compute samples from a mean zero Gaussian with covariance W_u^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            dim = 0;
            for i = 1:this.n_c
                dim = dim + length(this.I{i});
            end
            u_out = zeros(dim,num_samples);
            for i = 1:this.n_c  
                u_out(this.I{i},:) = this.u_prior_interface_cell{i}.Sample_with_Covariance_W_u_Inverse(num_samples);
            end
        end

        % Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            dim = 0;
            for i = 1:this.n_c
                dim = dim + length(this.I{i});
            end
            u_out = zeros(dim,num_samples);
            for i = 1:this.n_c  
                u_out(this.I{i},:) = this.u_prior_interface_cell{i}.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samples,scalar);
            end
        end

        function this = MD_Multi_State_u_Prior_Interface(u_prior_interface_cell)
            this.u_prior_interface_cell = u_prior_interface_cell;
            this.n_c = length(u_prior_interface_cell);
            this.I = cell(this.n_c,1);
            for i = 1:this.n_c
                this.I{i} = this.u_prior_interface_cell{i}.hyperparams.data_interface.Separate_State_Components(i);
            end
        end

    end

end