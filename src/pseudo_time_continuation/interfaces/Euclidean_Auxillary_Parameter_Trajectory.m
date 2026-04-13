classdef Euclidean_Auxillary_Parameter_Trajectory < Auxillary_Parameter_Trajectory

    properties
        theta_bar
        theta_star
        dtheta
    end

    methods

        function this = Euclidean_Auxillary_Parameter_Trajectory(N, theta_bar, theta_star)
            this@Auxillary_Parameter_Trajectory(N);
            this.theta_bar = theta_bar;
            this.theta_star = theta_star;
            this.dtheta = theta_star - theta_bar;
        end

        function [dtheta] = Get_dtheta(this)
            dtheta = this.dtheta;
        end

        function [theta_n] = Get_theta_n(this, n)
            if (n < 0) || (n > this.N)
                disp('The time index n is out of range');
            end

            theta_n = this.theta_bar + (n / this.N) * this.dtheta;
        end

    end
end
