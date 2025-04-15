classdef MD_Vector_z_Prior_Interface < MD_Scaled_z_Prior_Interface

    properties
        M_z
        R_z
        n_z
        z_hyperparam_interface
        determine_z_hyperparams
    end

    %% Implementation of base class virtual functions
    methods (Access = public)

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M_z * z_in;
        end

        function [z_out] = Apply_W_z_Acute_Inverse(this, z_in)
            tmp = linsolve(this.R_z', z_in);
            z_out = linsolve(this.R_z, tmp);
        end

        function [z_out] = Sample_with_Covariance_W_z_Acute_Inverse(this, num_samples)
            z_out = linsolve(this.R_z, randn(this.n_z,num_samples));
        end

        function [z_out] = Apply_W_z_Acute(this, z_in)
            z_out = this.M_z * z_in;
        end

    end

    %% Constructor
    methods

        function this = MD_Vector_z_Prior_Interface(M_z, data_interface, z_hyperparam_interface, u_prior_interface)
            arguments
                M_z (:,:) {mustBeNumeric}
                data_interface MD_Data_Interface
                z_hyperparam_interface MD_z_Hyperparameter_Interface
                u_prior_interface MD_u_Prior_Interface
            end
            this@MD_Scaled_z_Prior_Interface(z_hyperparam_interface.alpha_z);
            this.M_z = M_z;
            this.R_z = chol(M_z);
            this.n_z = size(M_z,1);

            this.z_hyperparam_interface = z_hyperparam_interface;
            this.determine_z_hyperparams = MD_Determine_z_Hyperparameters(data_interface, z_hyperparam_interface, u_prior_interface);

            if this.z_hyperparam_interface.alpha_z == 0.0
                this.determine_z_hyperparams.Determine_alpha_z(this);
            end
            this.Set_alpha_z(this.z_hyperparam_interface.alpha_z);

        end

    end

end
