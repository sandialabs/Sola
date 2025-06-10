classdef E_z_Inv_GSVD < Randomized_GSVD

    properties
        z_prior_interface
    end

    %% Implementation of base class virtual functions
    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            tmp = this.z_prior_interface.Apply_M_z(vec_in);
            vec_out = this.z_prior_interface.Apply_E_z_Inverse(tmp);
        end

        function [vec_out] = Apply_Operator_Transpose(this, vec_in)
            tmp = this.z_prior_interface.Apply_E_z_Inverse_Transpose(vec_in);
            vec_out = this.z_prior_interface.Apply_M_z(tmp);
        end

        function [vec_out] = Apply_Input_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.z_prior_interface.Apply_M_z_Inverse(vec_in);
        end

        function [vec_out] = Apply_Output_Weighting_Operator(this, vec_in)
            vec_out = this.z_prior_interface.Apply_M_z(vec_in);
        end

    end

    %% Constructor
    methods

        function this = E_z_Inv_GSVD(z_prior_interface, z)
            arguments
                z_prior_interface MD_z_Prior_Interface
                z (:, 1) double
            end
            this@Randomized_GSVD(z, z);
            this.z_prior_interface = z_prior_interface;
        end

    end

end
