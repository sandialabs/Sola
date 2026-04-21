%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Tutorial_2B_Constraint_AD < Dynamic_Constraint_AD
    % Constraint for an orbital dynamics problem.
    %
    % dy/dt = (y2, y1 y4^2 - k/y1^2 + q2, y4, (q1 -2 y2 y4) / y1)
    % y(0) = z
    %
    % Here y = (r, r', w, w') represents radius, change in radius, angular
    % momentum, and change in angular momentum of an object. The constants
    % r and w are the target radius and angular momentum. The controls
    % (q1(t), q2(t)) are thrust from the satellite in the radial and
    % tangential directions, respectively.

    properties
        k            % Constant of proportionality.
        n_q          % Number of controls at each fixed time.
        y0           % Initial condition.
    end

    methods (Access = public)

        function this = Tutorial_2B_Constraint_AD(T, n_t, r0, w0, k)
            % Constructor.
            %
            % Parameters
            % ----------
            % T
            %   Final simulation time.
            % n_t
            %   Number of time steps.
            % r0
            %   Initial radius of the satellite.
            % w0
            %   Initial angular velocity of the satellite.
            n_q = 2;
            n_z = n_q * (n_t - 1);
            this = this@Dynamic_Constraint_AD(4, n_z, T, n_t);
            this.y0 = [r0; 0; 0; w0];
            this.k = k;
            this.n_q = n_q;
        end

        function [f] = f_AD(this, y, z, t)
            I = this.Input_Indices(t);
            q = z(I);
            f = [y(2)
                 y(1) * y(4)^2 - (this.k / y(1)^2) + q(2)
                 y(4)
                 (q(1) - 2 * y(2) * y(4)) / y(1)];
        end

        function [h] = h_AD(this, ~)
            h = this.y0;
        end

    end

    methods (Access = private)

        function [mask] = Input_Indices(this, t)
            % Get the indices of the control at time t, i.e.,
            % I = Input_Indices(t) --> q(t) = z(I).

            [~, t_index] = min(abs(t - this.t_mesh));
            if t_index == 1
                error('no control at initial time!');
            end
            idx = t_index - 1;
            mask = (this.n_q * (idx - 1) + 1):(this.n_q * idx);
        end

    end
end
