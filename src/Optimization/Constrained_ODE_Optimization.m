%%
% Solve the optimiation problem
%
% $$\min_z \int_0^T g(y(z)(t),t)dt + R(z)$$
%
% where
% $y(z):[0,T] \to \mathbf{R}^m$ solves the odinary differential equation
%
% $\frac{du}{dt}=f(u,z)$
%
% $u(0)=h(z)$
%
% where
%
% * $T > 0$
% * $u(t) \in \mathbf{R}^m$
% * $z \in \mathbf{R}^n$
% * $g:\mathbf{R}^m \times [0,T] \to \mathbf{R}$
% * $R:\mathbf{R}^n \to \mathbf{R}$
%
% We derive off <Constrained_Optimization.html Constrained_Optimization> so that
%
% $c(u,z)=0$ corresponds to the discretization of the ODE system and
% $J(u,z)$ corresponds to applying the trapazoid rule to approximate the
% integral of $g$ in the objective
classdef Constrained_ODE_Optimization < Constrained_Optimization

    %%
    % Member properties:
    %
    % * m: the dimension of the state $u(t)$
    % * n: the dimension of the control $z$
    % * $T$: the final time
    % * N: the number of nodes in the time mesh
    % * t_mesh: the time mesh (a vector of length N)
    % * w: a vector of quadrature weights for time integration
    % * time_step_solver_options: options set for fsolve in time step
    properties
        m;
        n;
        T;
        N;
        t_mesh;
        w;
        time_step_solver_options;
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions for gradient computation

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * t: the time in the interval $[0,T]$
        %
        % Desciption:
        %
        % * Evaluate the time instance objective and its gradient.
        %
        % Output:
        %
        % * val: $g(u(t),t) \in \mathbf{R}$
        % * grad_u: $\nabla_u g(u(t),t) \in \mathbf{R}^m$
        [val, grad_u] = Time_Instance_Objective(obj,u,t);

        %%
        % Input:
        %
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Evaluate the regularization objective and gradient.
        %
        % Output:
        %
        % * val: $R(z) \in \mathbf{R}$
        % * grad_z: $\nabla_z R(z) \in \mathbf{R}^n$
        [val,grad_z] = Regularization_Objective(obj,z);

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Evaluate the ODE right hand side and its jacobians.
        %
        % Output:
        %
        % * f: $f(u(t),z) \in \mathbf{R}^m$
        % * f_u: $f_u(u(t),z,t) \in \mathbf{R}^{m \times m}$
        % * f_z: $f_z(u(t),z,t) \in \mathbf{R}^{m \times n}$
        [f, f_u, f_z] = Time_Instance_RHS(obj,u,z,t);

        %%
        % Input
        %
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Evaluate the initial condition and its jacobian.
        %
        % Output:
        %
        % * h: $h(z)\in \mathbf{R}^m$
        % * h_z: $h_z(z)\in \mathbf{R}^{m \times n}$
        [h, h_z] = Initial_Condition(obj,z);

        %% Pure virtual functions for hessian-vector product computation

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^m$
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\nabla_{u,u} g(u(t),t)v \in \mathbf{R}^m$
        [Mv] = Time_Instance_Objective_uu_Apply(obj,v,u,t);

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^n$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\nabla_{z,z} R(z)v \in \mathbf{R}^n$
        [Mv] = Regularization_Objective_zz_Apply(obj,v,z);

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^m$
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        % * lambda: the adjoint state $\lambda(t) \in \mathbf{R}^m$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\lambda(t)^T f_{u,u}(u(t),z,t)v \in \mathbf{R}^m$
        [Mv] = Time_Instance_RHS_uu_Apply(obj,v,u,z,t,lambda);

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^n$
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        % * lambda: the adjoint state $\lambda(t) \in \mathbf{R}^m$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\lambda(t)^T f_{u,z}(u(t),z,t)v \in \mathbf{R}^m$
        [Mv] = Time_Instance_RHS_uz_Apply(obj,v,u,z,t,lambda);

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^m$
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        % * lambda: the adjoint state $\lambda(t) \in \mathbf{R}^m$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\lambda(t)^T f_{z,u}(u(t),z,t)v \in \mathbf{R}^n$
        [Mv] = Time_Instance_RHS_zu_Apply(obj,v,u,z,t,lambda);

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^n$
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        % * lambda: the adjoint state $\lambda(t) \in \mathbf{R}^m$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\lambda(t)^T f_{z,z}(u(t),z)v \in \mathbf{R}^n$
        [Mv] = Time_Instance_RHS_zz_Apply(obj,v,u,z,t,lambda);

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^n$
        % * z: the control $z \in \mathbf{R}^n$
        % * lambda: the adjoint state $\lambda(t) \in \mathbf{R}^m$
        %
        % Description:
        %
        % * Evaluate the hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\lambda(t)^T h_{z,z}(z)v \in \mathbf{R}^n$
        [Mv] = Initial_Condition_zz_Apply(obj,v,z,lambda);

    end

    methods (Access = public)

        %% Instantiation of base class pure virtual functions for gradient computation

        %%
        % Input:
        %
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Evaluate the discretized objective and its gradients.
        %
        % Output:
        %
        % * val: $J(u,z) \in \mathbf{R}$
        % * grad_u: $\nabla_u J(u,z) \in \mathbf{R}^{mN}$
        % * grad_z: $\nabla_z J(u,z) \in \mathbf{R}^n$
        function [val, grad_u, grad_z] = Objective(obj,u,z)
            val = 0;
            grad_u = 0*u;
            for k = 1:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
               [valk,gradk] = obj.Time_Instance_Objective(u(I),obj.t_mesh(k));
               val = val + obj.w(k)*valk;
               grad_u(I) = obj.w(k)*gradk;
            end
            [valk,grad_z] = obj.Regularization_Objective(z);
            val = val + valk;
        end

        %%
        % Input:
        %
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Solve the system of equation (discretized ODE) $c(u,z)=0$ for $u$.
        %
        % Output:
        %
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        function [u] = State_Solve(obj,z) % Input z, evaluate u=S(z)
            u = zeros(obj.m*obj.N,1);
            u(1:obj.m) = obj.Initial_Condition(z);
            for k = 2:obj.N
               Im = ((k-2)*obj.m+1):((k-1)*obj.m);
               I = ((k-1)*obj.m+1):(k*obj.m);
               u(I) = obj.State_Eq_Time_Step(u(Im),z,obj.t_mesh(k),obj.t_mesh(k-1));
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Solve the system of linear equations $c_u(u,z)^T x = v$ for $x \in \mathbf{R}^{mN}$
        %
        % Output:
        %
        % * Mv: $c_u(u,z)^{-T}v \in \mathbf{R}^{mN}$
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
            I = ((obj.N-1)*obj.m+1):(obj.N*obj.m);
            dt = obj.t_mesh(end)-obj.t_mesh(end-1);
            Mv(I,:) = obj.Linearized_Adjoint_Time_Step_Solve(v(I,:),u(I),z,obj.t_mesh(end),dt);

            for k = (obj.N-1):-1:2
                Im = (k*obj.m+1):((k+1)*obj.m);
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k+1)-obj.t_mesh(k);
                Mv(I,:) = obj.Linearized_Adjoint_Time_Step_Solve(v(I,:) + Mv(Im,:),u(I),z,obj.t_mesh(k),dt);
            end
            I = 1:(obj.m);
            Im = (obj.m+1):(2*obj.m);
            Mv(I,:) = Mv(Im,:) + v(I,:);
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Compute the control jacobian transpose matrix-vector product.
        %
        % Output:
        %
        % * Mv: $c_z(u,z)^Tv \in \mathbf{R}^{n}$
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            [~, h_z] = obj.Initial_Condition(z);
            Mv = -h_z'*v(1:obj.m,:);
            for k = 2:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                [~, ~,f_z] = obj.Time_Instance_RHS(u(I),z,obj.t_mesh(k));
                Mv = Mv - dt*f_z'*v(I,:);
            end
        end

        %% Instantiation of base class pure virtual functions for hessian-vector product computation

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Solve the system of linear equations $c_u(u,z) x = v$ for $x \in \mathbf{R}^{mN}$.
        %
        % Output:
        %
        % * Mv: $c_u(u,z)^{-1}v \in \mathbf{R}^{mN}$
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
            I = 1:obj.m;
            Mv(I,:) = v(I,:);
            for k = 2:obj.N
                Im = ((k-2)*obj.m+1):((k-1)*obj.m);
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                Mv(I,:) = obj.Linearized_Time_Step_Solve(v(I,:) + Mv(Im,:),u(I),z,obj.t_mesh(k),dt);
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{n}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Compute the control jacobian matrix-vector product.
        %
        % Output:
        %
        % * Mv: $c_z(u,z)v \in \mathbf{R}^{mN}$
        function [Mv] = c_z_Apply(obj,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
            I = 1:obj.m;
            [~, h_z] = obj.Initial_Condition(z);
            Mv(I,:) = -h_z*v;
            for k = 2:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                [~, ~,f_z] = obj.Time_Instance_RHS(u(I),z,obj.t_mesh(k));
                Mv(I,:) = Mv(I,:) - dt*f_z*v;
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        % * lambda: the discretized adjoint state $\lambda=(\lambda(t_1),\lambda(t_2),\dots,\lambda(t_N))^T \in \mathbf{R}^{mN}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{u,u} (\lambda^Tc(u,z))v$.
        %
        % Output:
        %
        % * Mv: $\lambda^Tc_{u,u}(u,z)v \in \mathbf{R}^{mN}$
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
            for k = 2:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                f_uu = obj.Time_Instance_RHS_uu_Apply(v(I,:),u(I),z,obj.t_mesh(k),lambda(I));
                Mv(I,:) = -dt*f_uu;
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{n}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        % * lambda: the discretized adjoint state $\lambda=(\lambda(t_1),\lambda(t_2),\dots,\lambda(t_N))^T \in \mathbf{R}^{mN}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{u,z} (\lambda^Tc(u,z))v$.
        %
        % Output:
        %
        % * Mv: $\lambda^Tc_{u,z}(u,z)v \in \mathbf{R}^{mN}$
        function [Mv] = c_uz_Apply(obj,v,u,z,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
            for k = 2:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                f_uz = obj.Time_Instance_RHS_uz_Apply(v,u(I),z,obj.t_mesh(k),lambda(I));
                Mv(I,:) = -dt*f_uz;
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        % * lambda: the discretized adjoint state $\lambda=(\lambda(t_1),\lambda(t_2),\dots,\lambda(t_N))^T \in \mathbf{R}^{mN}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{z,u} (\lambda^Tc(u,z))v$.
        %
        % Output:
        %
        % * Mv: $\lambda^Tc_{z,u}(u,z)v \in \mathbf{R}^{n}$
        function [Mv] = c_zu_Apply(obj,v,u,z,lambda)
            num_vecs = size(v,2);
            Mv = zeros(obj.n,num_vecs);
            for k = 2:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                f_zu = obj.Time_Instance_RHS_zu_Apply(v(I,:),u(I),z,obj.t_mesh(k),lambda(I));
                Mv = Mv - dt*f_zu;
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{n}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        % * lambda: the discretized adjoint state $\lambda=(\lambda(t_1),\lambda(t_2),\dots,\lambda(t_N))^T \in \mathbf{R}^{mN}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{z,z} (\lambda^Tc(u,z))v$.
        %
        % Output:
        %
        % * Mv: $\lambda^Tc_{z,z}(u,z)v \in \mathbf{R}^{n}$
        function [Mv] = c_zz_Apply(obj,v,u,z,lambda)
            Mv = -obj.Initial_Condition_zz_Apply(v,z,lambda(1:obj.m));
            for k = 2:obj.N
                I = ((k-1)*obj.m+1):(k*obj.m);
                dt = obj.t_mesh(k)-obj.t_mesh(k-1);
                f_zz = obj.Time_Instance_RHS_zz_Apply(v,u(I),z,obj.t_mesh(k),lambda(I));
                Mv = Mv - dt*f_zz;
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\nabla_{u,u}J(u,z)v \in \mathbf{R}^{mN}$
        function [Mv] = Objective_uu_Apply(obj,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
            for k = 1:obj.N
               I = ((k-1)*obj.m+1):(k*obj.m);
               Mv(I,:) = obj.w(k)*obj.Time_Instance_Objective_uu_Apply(v(I,:),u(I),obj.t_mesh(k));
            end
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{n}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\nabla_{u,z}J(u,z)v \in \mathbf{R}^{mN}$
        function [Mv] = Objective_uz_Apply(obj,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(obj.m*obj.N,num_vecs);
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{mN}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\nabla_{z,u}J(u,z)v \in \mathbf{R}^{n}$
        function [Mv] = Objective_zu_Apply(obj,v,u,z)
            num_vecs = size(v,2);
            Mv = zeros(obj.n,num_vecs);
        end

        %%
        % Input:
        %
        % * v: a direction $v \in \mathbf{R}^{n}$
        % * u: the discretized state $u=(u(t_1),u(t_2),\dots,u(t_N))^T \in \mathbf{R}^{mN}$
        % * z: the control $z \in \mathbf{R}^n$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        %
        % Output:
        %
        % * Mv: $\nabla_{z,z}J(u,z)v \in \mathbf{R}^{n}$
        function [Mv] = Objective_zz_Apply(obj,v,u,z)
            Mv = obj.Regularization_Objective_zz_Apply(v,z);
        end

    end

    %%
    % Input:
    %
    % * m: the dimension of the state $u(t)$
    % * n: the dimension of the control $z$
    % * $T$: the final time
    % * N: the number of nodes in the time mesh
    %
    % Description:
    %
    % * Constructor which sets dimensions and time mesh data structures.
    methods (Access = public)
        function obj = Constrained_ODE_Optimization(m,n,T,N)
            obj.m = m; % State dimension
            obj.n = n;
            obj.T = T; % Specify the final time
            obj.N = N; % Number of time nodes
            obj.t_mesh = linspace(0,T,N)';
            w = ones(N,1);
            w(2:end-1) = 2;
            obj.w = T*w/sum(w);
            obj.time_step_solver_options = optimoptions('fsolve','Display','none','SpecifyObjectiveGradient',true);
        end

        %% Finite difference tests

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Finite difference check for $f_u(u,z,t)$.
        %
        % Output:
        %
        % * diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Jacobian_u_Check(obj,u,z,t)
            [f, f_u] = obj.Time_Instance_RHS(u,z,t);
            v = randn(obj.m,1);
            v = v/norm(v);
            fv = f_u*v;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_fv = zeros(obj.m,p);
            diffs = zeros(p,1);
            for k = 1:p
                [fk] = obj.Time_Instance_RHS(u+h(k)*v,z,t);
                fd_fv(:,k) = (fk-f)/h(k);
                diffs(k) = norm(fd_fv(:,k)-fv)/norm(fv);
            end
            if obj.verbose
                disp('State jacobian finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
            disp(' ')
            end
        end

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Finite difference check for $f_z(u,z,t)$.
        %
        % Output:
        %
        % * diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Jacobian_z_Check(obj,u,z,t)
            [f, ~, f_z] = obj.Time_Instance_RHS(u,z,t);
            v = randn(obj.n,1);
            v = v/norm(v);
            fv = f_z*v;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_fv = zeros(obj.m,p);
            diffs = zeros(p,1);
            for k = 1:p
                [fk] = obj.Time_Instance_RHS(u,z+h(k)*v,t);
                fd_fv(:,k) = (fk-f)/h(k);
                diffs(k) = norm(fd_fv(:,k)-fv)/norm(fv);
            end
            if obj.verbose
                disp('Control jacobian finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Finite difference check for $f_{u,u}(u,z,t)$.
        %
        % Output:
        %
        % * diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_uu_Check(obj,u,z,t)
            v = randn(obj.m,1);
            v = v/norm(v);
            lambda = randn(obj.m,1);
            Mv = obj.Time_Instance_RHS_uu_Apply(v,u,z,t,lambda);
            [~, f_u] = obj.Time_Instance_RHS(u,z,t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(obj.m,p);
            diffs = zeros(p,1);
            for k = 1:p
                [~,f_uk] = obj.Time_Instance_RHS(u+h(k)*v,z,t);
                fd_Mv(:,k) = (f_uk'*lambda-f_u'*lambda)/h(k);
                diffs(k) = norm(fd_Mv(:,k)-Mv)/norm(Mv);
            end
            if obj.verbose
                disp('Hessian uu finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Finite difference check for $f_{u,z}(u,z,t)$.
        %
        % Output:
        %
        % * diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_uz_Check(obj,u,z,t)
            v = randn(obj.n,1);
            v = v/norm(v);
            lambda = randn(obj.m,1);
            Mv = obj.Time_Instance_RHS_uz_Apply(v,u,z,t,lambda);
            [~, f_u] = obj.Time_Instance_RHS(u,z,t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(obj.m,p);
            diffs = zeros(p,1);
            for k = 1:p
                [~,f_uk] = obj.Time_Instance_RHS(u,z+h(k)*v,t);
                fd_Mv(:,k) = (f_uk'*lambda-f_u'*lambda)/h(k);
                diffs(k) = norm(fd_Mv(:,k)-Mv)/norm(Mv);
            end
            if obj.verbose
                disp('Hessian uz finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Finite difference check for $f_{z,u}(u,z,t)$.
        %
        % Output:
        %
        % * diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_zu_Check(obj,u,z,t)
            v = randn(obj.m,1);
            v = v/norm(v);
            lambda = randn(obj.m,1);
            Mv = obj.Time_Instance_RHS_zu_Apply(v,u,z,t,lambda);
            [~, ~, f_z] = obj.Time_Instance_RHS(u,z,t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(obj.n,p);
            diffs = zeros(p,1);
            for k = 1:p
                [~,~,f_zk] = obj.Time_Instance_RHS(u+h(k)*v,z,t);
                fd_Mv(:,k) = (f_zk'*lambda-f_z'*lambda)/h(k);
                diffs(k) = norm(fd_Mv(:,k)-Mv)/norm(Mv);
            end
            if obj.verbose
                disp('Hessian zu finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end

        %%
        % Input:
        %
        % * u: the state $u(t) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * t: the time in the interval $[0,T]$
        %
        % Description:
        %
        % * Finite difference check for $f_{z,z}(u,z,t)$.
        %
        % Output:
        %
        % * diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_zz_Check(obj,u,z,t)
            v = randn(obj.n,1);
            v = v/norm(v);
            lambda = randn(obj.m,1);
            Mv = obj.Time_Instance_RHS_zz_Apply(v,u,z,t,lambda);
            [~, ~, f_z] = obj.Time_Instance_RHS(u,z,t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(obj.n,p);
            diffs = zeros(p,1);
            for k = 1:p
                [~,~,f_zk] = obj.Time_Instance_RHS(u,z+h(k)*v,t);
                fd_Mv(:,k) = (f_zk'*lambda-f_z'*lambda)/h(k);
                diffs(k) = norm(fd_Mv(:,k)-Mv)/norm(Mv);
            end
            if obj.verbose
                disp('Hessian zz finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end

    end

    %% Private functions
     methods (Access = protected)

        %%
        % Input:
        %
        % * ukm: the state $u(t_{k-1}) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * tk: next time step $t_k$
        % * tdm: previous time step $t_{k-1}$
        %
        % Description:
        %
        % * Advance in time by solving $u_k-u_{k-1} - dt f(u_k,z,t_k)=0$.
        %
        % Output:
        %
        % * uk: the state $u(t_k) \in \mathbf{R}^m$
         function [uk] = State_Eq_Time_Step(obj,ukm,z,tk,tkm)
               dt = tk-tkm;
               uk = fsolve(@(u)obj.Nonlinear_Step(u,ukm,z,tk,dt),ukm,obj.time_step_solver_options);
         end

         %%
         % Input:
         %
         % * v: the direction $v \in \mathbf{R}^m$
         % * u: the state $u \in \mathbf{R}^m$
         % * z: the control $z \in \mathbf{R}^n$
         % * tk: the time in the interval $[0,T]$
         % * dt: time step size
         %
         % Description:
         %
         % * Solve the system of linear equations $(I_m - dt*f_u(u,z,t_k))x=v$
         %
         % Output:
         %
         % * Mv: $(I_m - dt*f_u(u,z,t_k))^{-1}v$
         function [Mv] = Linearized_Time_Step_Solve(obj,v,u,z,tk,dt)
            [~, f_u] = obj.Time_Instance_RHS(u,z,tk);
            A = eye(obj.m) - dt*f_u;
            Mv = linsolve(A,v);
         end

         %%
         % Input:
         %
         % * v: the direction $v \in \mathbf{R}^m$
         % * u: the state $u \in \mathbf{R}^m$
         % * z: the control $z \in \mathbf{R}^n$
         % * tk: the time in the interval $[0,T]$
         % * dt: time step size
         %
         % Description:
         %
         % * Solve the system of linear equations $(I_m - dt*f_u(u,z,t_k)^T)x=v$
         %
         % Output:
         %
         % * Mv: $(I_m - dt*f_u(u,z,t_k)^T)^{-1}v$
         function [Mv] = Linearized_Adjoint_Time_Step_Solve(obj,v,u,z,tk,dt)
             [~, f_u] = obj.Time_Instance_RHS(u,z,tk);
             A = eye(obj.m) - dt*f_u';
             Mv = linsolve(A,v);
         end

        %%
        % Input:
        %
        % * uk: the state $u(t_k) \in \mathbf{R}^m$
        % * ukm: the state $u(t_{k-1}) \in \mathbf{R}^m$
        % * z: the control $z \in \mathbf{R}^n$
        % * tk: the time in the interval $[0,T]$
        % * dt: time difference $t_k-t_{k-1}$
        %
        % Description:
        %
        % * Evaluate the residual and jacobian for the system of nonlinear equations $u_k-u_{k-1} - dt f_u(u_k,z,t_k)=0$.
        %
        % Output:
        %
        % * f: value of the residual $u_k-u_{k-1} - dt f(u_k,z,t_k) \in
        % \mathbf{R}^m$
        % * J: value of the $m \times m$ state jacobian of the residual
         function [f,J] = Nonlinear_Step(obj,uk,ukm,z,tk,dt)
            [val, val_u] = obj.Time_Instance_RHS(uk,z,tk);
            f = uk - ukm - dt*val;
            J = eye(obj.m,obj.m);
            J = J - dt*val_u;
         end

     end

end
