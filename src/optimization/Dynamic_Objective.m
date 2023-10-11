% Define the objective function int_0^T g(y(z)(t),t)dt + R(z) 
% where
% T > 0
% y(t) in R^m
% z in R^n
% g:R^m \times [0,T] \to R
% R:R^n \to R
%
% We derive off Objective so that J(u,z) corresponds to applying the trapazoid rule 
% to approximate the integral of g in the objective
classdef Dynamic_Objective < Objective

    properties
        m; % dimension of the ODE state y(t)
        n; % dimension of the control z
        N; % number of nodes in the time mesh
        t_mesh; % time mesh (a vector of length N)
        w; % vector of quadrature weights for time integration
    end 

    methods (Abstract, Access = public)

        % Input:
        % y: the ODE state y(t) in R^m
        % t: the time in the interval [0,T]
        % Output:
        % val: g(y(t),t) in R
        % grad_u: \nabla_y g(y(t),t) in R^m
        [val, grad_y] = Time_Instance_Objective(this,y,t);

        % Input:
        % z: the control z in R^n
        % Output:
        % val: R(z) in R
        % grad_z: \nabla_z R(z) in R^n
        [val,grad_z] = Regularization_Objective(this,z);

        % Input:
        % v: a direction v in R^m
        % y: the ODE state y(t) in R^m
        % t: the time in the interval [0,T]
        % Output:
        % Mv: \nabla_{y,y} g(y(t),t)v in R^m
        [Mv] = Time_Instance_Objective_yy_Apply(this,v,y,t);

        % Input:
        % v: a direction v in R^n
        % z: the control z in R^n
        % Output:
        % Mv: \nabla_{z,z} R(z)v in R^n
        [Mv] = Regularization_Objective_zz_Apply(this,v,z);
        
    end

    methods (Access = public)

        %% Instantiation of base class pure virtual functions
        function [val, grad_u, grad_z] = J(this,u,z)
            val = 0;
            grad_u = 0*u;
            for k = 1:this.N
                I = ((k-1)*this.m+1):(k*this.m);          % y_{k} = u(I)
               [valk,gradk] = this.Time_Instance_Objective(u(I),this.t_mesh(k));
               val = val + this.w(k)*valk;
               grad_u(I) = this.w(k)*gradk;
            end
            [valk,grad_z] = this.Regularization_Objective(z);
            val = val + valk;
        end

        function [Mv] = J_uu_Apply(this,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(this.m*this.N,num_vecs);
            for k = 1:this.N
               I = ((k-1)*this.m+1):(k*this.m);          % y_{k} = u(I)
               Mv(I,:) = this.w(k)*this.Time_Instance_Objective_yy_Apply(v(I,:),u(I),this.t_mesh(k));
            end
        end

        function [Mv] = J_uz_Apply(this,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(this.m*this.N,num_vecs);
        end

        function [Mv] = J_zu_Apply(this,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(this.n,num_vecs);
        end

        function [Mv] = J_zz_Apply(this,v,u,z)
            Mv = this.Regularization_Objective_zz_Apply(v,z);
        end

    end

    methods (Access = public)

        % Input:
        % m: the dimension of the ODE state y(t)
        % n: the dimension of the control z
        % T: the final time
        % N: the number of nodes in the time mesh
        function this = Dynamic_Objective(m,n,T,N)
            this.m = m;                      % ODE state dimension
            this.n = n;                      % control dimension
            this.N = N;                      % Number of time nodes
            this.t_mesh = linspace(0,T,N)';  % Discrete time domain
            w = ones(N,1);
            w(2:end-1) = 2;
            this.w = T*w/sum(w);             % Spatial weights
        end

    end

end
