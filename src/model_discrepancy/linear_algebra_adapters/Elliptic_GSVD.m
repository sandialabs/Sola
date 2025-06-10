classdef Elliptic_GSVD < Randomized_GSVD

    properties
        u_prior_interface
    end

    %% Implementation of base class virtual functions
    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            vec_out = this.u_prior_interface.Apply_E_u_Inverse(vec_in);
        end

        function [vec_out] = Apply_Operator_Transpose(this, vec_in)
            vec_out = this.u_prior_interface.Apply_E_u_Inverse_Transpose(vec_in);
        end

        function [vec_out] = Apply_Input_Weighting_Operator_Inverse(this, vec_in)
            vec_out = this.u_prior_interface.Apply_M_u(vec_in);
        end

        function [vec_out] = Apply_Output_Weighting_Operator(this, vec_in)
            vec_out = this.u_prior_interface.Apply_M_u(vec_in);
        end

    end

    %% Constructor
    methods

        function this = Elliptic_GSVD(u_prior_interface, u_in, u_out)
            arguments
                u_prior_interface MD_u_Prior_Interface
                u_in (:, 1) double
                u_out (:, 1) double
            end
            this@Randomized_GSVD(u_in, u_out);
            this.u_prior_interface = u_prior_interface;
        end

    end

end
