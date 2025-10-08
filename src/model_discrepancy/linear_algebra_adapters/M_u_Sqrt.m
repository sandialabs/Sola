classdef M_u_Sqrt < Matrix_Sqrt

    properties
        u_prior_interface
    end

    %% Implementation of base class functions
    methods (Access = public)

        function [vec_out] = Matrix_Apply(this, vec_in)
            vec_out = this.u_prior_interface.Apply_M_u(vec_in);
        end

    end

    %% Constructor
    methods

        function this = M_u_Sqrt(u_prior_interface)
            arguments
                u_prior_interface MD_u_Prior_Interface
            end
            this.u_prior_interface = u_prior_interface;
        end

    end

end
