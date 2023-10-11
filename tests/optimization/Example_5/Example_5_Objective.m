classdef Example_5_Objective < Dynamic_Objective

    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the ordinary differential equation
    % dy/dt = [2*t*y_1 ; 3*t^2*y_2 + t^2 - z(t)]
    % y(0) = [1, 1]
    % g(y) = (y_1-exp(t^2))^2 + (y_2-exp(t^3))^2
    % R(z) = int_0^T (z(t) - t^2)^2dt

    properties
        z_time_mesh;
        weights;
        beta_reg;
    end

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this,y,t)
            val = (y(1)-exp(t^2))^2 + (y(2)-exp(t^3))^2;
            grad_y = zeros(2,1);
            grad_y(1) = 2*(y(1)-exp(t^2));
            grad_y(2) = 2*(y(2)-exp(t^3));
        end

        function [val,grad_z] = Regularization_Objective(this,z)
            val = this.beta_reg*(z-this.z_time_mesh.^2)'*diag(this.weights)*(z-this.z_time_mesh.^2);
            grad_z = this.beta_reg*2*diag(this.weights)*(z-this.z_time_mesh.^2);
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this,v,y,t)
            A = 2*eye(2);
            Mv = A*v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(this,v,z)
            Mv = this.beta_reg*2*diag(this.weights)*v;
        end

    end

    methods (Access = public)
        function this = Example_5_Objective(m,n,T,N)
            this = this@Dynamic_Objective(m,n,T,N);
            this.z_time_mesh = linspace(0,T,n+1)';
            this.z_time_mesh = this.z_time_mesh(2:end);

            weights = ones(this.n+1,1);
            weights(1) = .5; weights(end) = .5;
            weights = T*weights/sum(weights);
            weights = weights(2:end);
            this.weights = weights;
            this.beta_reg = 10^-4;
        end

        function [w] = Temporal_Weights(this,t)
           w = (this.z_time_mesh-t)/(this.z_time_mesh(2)-this.z_time_mesh(1));
           Im = intersect(find(w<=0),find(abs(w)<=1));
           Ip = intersect(find(w>0),find(abs(w)<=1));
           I = find(abs(w)>1);
           w(I) = 0;
           w(Im) = 1+w(Im);
           w(Ip) = 1-w(Ip);
        end

    end

end
