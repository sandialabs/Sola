classdef M_z_Sqrt < Matrix_Sqrt

    properties
        z_prior_interface
    end

    methods (Access = public)

        function [vec_out] = Matrix_Apply(this, vec_in)
            vec_out = this.z_prior_interface.Apply_M_z(vec_in);
        end

    end

    methods

        function this = M_z_Sqrt(z_prior_interface)
            this.z_prior_interface = z_prior_interface;
        end

    end

end
