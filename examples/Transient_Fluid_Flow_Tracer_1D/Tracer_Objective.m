classdef Tracer_Objective < Objective

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            % Call the Python function J
            if nargout == 1
                val = py.fluid_flow_1d_lofi.J(z, u);
                grad_u = NaN;
                grad_z = NaN;
            elseif nargout == 3
                val = NaN;
                grad_u = NaN;
                grad_z = double(py.fluid_flow_1d_lofi.Jz(z, u));
                if isvector(grad_z)
                    grad_z = grad_z';
                end
            end
            % val = double(result{1}); % Objective value
            % grad_z = double(result{2})'; % Gradient with respect to z
            % grad_u = double(result{3})'; % Gradient with respect to u
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            % Call the Python function J_uu_apply
            Mv = py.fluid_flow_1d_lofi.J_uu_apply(z, u, v);
            Mv = double(Mv); % Convert Python array to MATLAB array
            if isvector(Mv)
                Mv = Mv';
            end
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(length(u), 1);
            disp(size(Mv));
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(length(z), 1);
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            % Call the Python function J_zz_apply
            Mv = py.fluid_flow_1d_lofi.J_zz_apply(z, u, v);
            Mv = double(Mv);
            if isvector(Mv)
                Mv = Mv';
            end
        end

    end
end
