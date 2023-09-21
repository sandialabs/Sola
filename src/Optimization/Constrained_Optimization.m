%%
% Solve the optimization problem
%
% $$\min_z J(S(z),z)$$
%
% where
% $$S(z)$ solves the constraint equation $c(u,z)=0$
%
% i.e. $$c(S(z),z)=0$ for all $$z$
%
% where 
%
% * $u \in \mathbf{R}^{n_u}$
% * $z \in \mathbf{R}^{n_z}$
% * $c(u,z) \in \mathbf{R}^{n_u}$

classdef Constrained_Optimization < handle
     
    %%
    % Member properties are the default optimizer settings for iterations and tolerances
    properties
        opt_tol;
        fun_tol;
        iteration_limit;
        step_tol;
        max_cg_iter;
        cg_tol;
        verbose;
    end
    
    methods (Abstract, Access = public)

         %% Pure virtual functions for gradient computation
         
         %%
         % Input:
         %
         % * u: the state $u \in \mathbf{R}^{n_u}$ 
         % * z: the control $z \in \mathbf{R}^{n_z}$
         %
         % Description:
         %
         % * Evaluate the objective and its gradients.
         %
         % Output:
         %
         % * $J(u,z) \in \mathbf{R}$
         % * $\nabla_u J(u,z)\in \mathbf{R}^{n_u}$
         % * $\nabla_z J(u,z)\in \mathbf{R}^{n_z}$
        [val, grad_u, grad_z] = Objective(obj,u,z); 
        
        %%
        % Input:
        %
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Solve the state equation, i.e. evaluate $S(z)$ where $c(S(z),z)=0$.
        %
        % Output:
        %
        % * u: $u=S(z)\in \mathbf{R}^{n_u}$
        [u] = State_Solve(obj,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Solve the linear system $c_u(u,z)^T x = v$ for $x$.
        % 
        % Output:
        %
        % * Mv: $c_u(u,z)^{-T}v \in \mathbf{R}^{n_u}$
        [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the constraint jacobian transpose matrix-vector product
        % 
        % Output:
        %
        % * Mv: $c_z(u,z)^{T}v \in \mathbf{R}^{n_z}$
        [Mv] = c_z_Transpose_Apply(obj,v,u,z); 
        
        %% Pure virtual functions for hessian-vector product computation
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Solve the linear system $c_u(u,z) x = v$ for $x$.
        % 
        % Output:
        %
        % * Mv: $c_u(u,z)^{-1}v \in \mathbf{R}^{n_u}$
        [Mv] = c_u_Inverse_Apply(obj,v,u,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_z}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the constraint jacobian matrix-vector product
        % 
        % Output:
        %
        % * Mv: $c_z(u,z)v \in \mathbf{R}^{n_u}$
        [Mv] = c_z_Apply(obj,v,u,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        % * lambda: the adjoint state $\lambda \in \mathbf{R}^{n_u}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{u,u} (\lambda^Tc(u,z))v$.
        % 
        % Output:
        %
        % * Mv: $\lambda^T c_{u,u}(u,z)v \in \mathbf{R}^{n_u}$
        [Mv] = c_uu_Apply(obj,v,u,z,lambda); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_z}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        % * lambda: the adjoint state $\lambda \in \mathbf{R}^{n_u}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{u,z} (\lambda^Tc(u,z))v$.
        % 
        % Output:
        %
        % * Mv: $\lambda^T c_{u,z}(u,z)v \in \mathbf{R}^{n_u}$
        [Mv] = c_uz_Apply(obj,v,u,z,lambda); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        % * lambda: the adjoint state $\lambda \in \mathbf{R}^{n_u}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{z,u} (\lambda^Tc(u,z))v$.
        % 
        % Output:
        %
        % * Mv: $\lambda^T c_{z,u}(u,z)v \in \mathbf{R}^{n_z}$
        [Mv] = c_zu_Apply(obj,v,u,z,lambda);
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_z}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        % * lambda: the adjoint state $\lambda \in \mathbf{R}^{n_u}$
        %
        % Description:
        %
        % * Compute the constraint hessian-vector product $\nabla_{z,z} (\lambda^Tc(u,z))v$.
        % 
        % Output:
        %
        % * Mv: $\lambda^T c_{z,z}(u,z)v \in \mathbf{R}^{n_z}$
        [Mv] = c_zz_Apply(obj,v,u,z,lambda);
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        % 
        % Output:
        %
        % * Mv: $\nabla_{u,u} J(u,z)v \in \mathbf{R}^{n_u}$
        [Mv] = Objective_uu_Apply(obj,v,u,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_z}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        % 
        % Output:
        %
        % * Mv: $\nabla_{u,z} J(u,z)v \in \mathbf{R}^{n_u}$
        [Mv] = Objective_uz_Apply(obj,v,u,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_u}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        % 
        % Output:
        %
        % * Mv: $\nabla_{z,u} J(u,z)v \in \mathbf{R}^{n_z}$
        [Mv] = Objective_zu_Apply(obj,v,u,z); 
        
        %%
        % Input:
        %
        % * v: a direction $v\in \mathbf{R}^{n_z}$
        % * u: the state $u \in \mathbf{R}^{n_u}$ 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the objective hessian-vector product.
        % 
        % Output:
        %
        % * Mv: $\nabla_{z,z} J(u,z)v \in \mathbf{R}^{n_z}$
        [Mv] = Objective_zz_Apply(obj,v,u,z); 
        
    end
    
    methods (Access = public)
        
        %%
        % Description:
        %
        % * Constructor sets the default optimization tolerances and iteration limits
        function obj = Constrained_Optimization( )
            obj.opt_tol = 10^-8;
            obj.fun_tol = 10^-6;
            obj.iteration_limit = 10^3;
            obj.step_tol = 10^-6;
            obj.max_cg_iter = 50;
            obj.cg_tol = 10^-4;
            obj.verbose = true;
        end
      
        %% Optimization functions
        %%
        % Input:
        % 
        % * z0: initial iterate control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Run trust region optimization with fminunc.
        %
        % Output:
        %
        % * u: the optimal state solution $u \in \mathbf{R}^{n_u}$
        % * z: the optimal control $z \in \mathbf{R}^{n_z}$
        function [u,z] = Optimize(obj,z0)
            HessMultFcn = @(hessian_data,v) obj.Jhat_hessVec(hessian_data,v);
            verb = 'iter-detailed';
            if obj.verbose == false
               verb = 'none'; 
            end
            options = optimoptions(@fminunc,'Display',verb,'Algorithm','trust-region','SpecifyObjectiveGradient',true,...
                                    'OptimalityTolerance',obj.opt_tol,'FunctionTolerance',obj.fun_tol,'MaxIterations',obj.iteration_limit,...
                                    'StepTolerance',obj.step_tol,'SubproblemAlgorithm','cg','MaxPCGIter',obj.max_cg_iter,...
                                    'TolPCG',obj.cg_tol,'HessianMultiplyFcn',HessMultFcn);
            z = fminunc(@(z)obj.Jhat(z),z0,options);
            u = obj.State_Solve(z);
        end
        
        %%
        % Input:
        % 
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Evaluate the reduced objective $\hat{J}(z)=J(S(z),z)$ and its gradient.
        %
        % Output:
        %
        % * val: $\hat{J}(z)$
        % * grad: $\nabla_z \hat{J}(z)$
        % * hessian_data: concatenation of state, control, and adjoint to pass to hessian-vector multiply function
        function [val,grad,hessian_data] = Jhat(obj,z)
            u = obj.State_Solve(z);
            [val, grad_u, grad_z] = obj.Objective(u,z);
            lambda = obj.c_u_Transpose_Inverse_Apply(-grad_u,u,z);
            grad = obj.c_z_Transpose_Apply(lambda,u,z);
            grad = grad + grad_z;
            hessian_data = [u;z;lambda];
        end
        
        %%
        % Input:
        % 
        % * hessian_data: output from Jhat function containing the state $u$, control $z$, and adjoint $\lambda$
        % * v: a direction $v \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Compute the reduced objective hessian-vector product.
        %
        % Output:
        %
        % * Hv: $\nabla_{z,z} \hat{J}(z)v$
        function [Hv] = Jhat_hessVec(obj,hessian_data,v)
            p = length(v);
            m = (length(hessian_data)-p)/2;
            u = hessian_data(1:m);
            z = hessian_data((m+1):(m+p));
            lambda = hessian_data((m+p+1):end);
            
            w = obj.c_z_Apply(v,u,z);
            mu = obj.c_u_Inverse_Apply(-w,u,z);
            yJ = obj.Objective_uu_Apply(mu,u,z) + obj.Objective_uz_Apply(v,u,z);
            yc = obj.c_uu_Apply(mu,u,z,lambda) + obj.c_uz_Apply(v,u,z,lambda);
            gamma = obj.c_u_Transpose_Inverse_Apply(-(yJ+yc),u,z);
            xJ = obj.Objective_zu_Apply(mu,u,z) + obj.Objective_zz_Apply(v,u,z);
            xc = obj.c_z_Transpose_Apply(gamma,u,z) + obj.c_zu_Apply(mu,u,z,lambda) + obj.c_zz_Apply(v,u,z,lambda);
            Hv = xJ + xc;
        end
        
        %% Finite difference test functions
        %%
        % Input:
        %
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Finite difference check for $\nabla_{z} \hat{J}(z)$.
        %
        % Output:
        % 
        % * diffs: vector of finite difference errors
        function [diffs] = Finite_Difference_Gradient_Check(obj,z)
            [val,grad] = obj.Jhat(z);
            n = length(grad);
            dz = randn(n,1);
            dz = dz/norm(dz);
            grad_dz = dz'*grad;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_grad = zeros(p,1);
            for k = 1:p
                valk = obj.Jhat(z+h(k)*dz);
                fd_grad(k) = (valk-val)/h(k);
            end
            
            diffs = abs(grad_dz-fd_grad)/abs(grad_dz);
            if obj.verbose
                disp('Gradient finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end
        
        %%
        % Input:
        %
        % * z: the control $z \in \mathbf{R}^{n_z}$
        %
        % Description:
        %
        % * Finite difference check for $\nabla_{z,z} \hat{J}(z)$.
        %
        % Output:
        % 
        % * diffs: vector of finite difference errors
        function [diffs] = Finite_Difference_Hessian_Check(obj,z)
            [~,grad,hessian_data] = obj.Jhat(z);
            n = length(grad);
            v = randn(n,1);
            v = v/norm(v);
            Hv = obj.Jhat_hessVec(hessian_data,v);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_hv = zeros(n,p);
            diffs = zeros(p,1);
            for k = 1:p
                [~,gradk] = obj.Jhat(z+h(k)*v);
                fd_hv(:,k) = (gradk-grad)/h(k);
                diffs(k) = norm(fd_hv(:,k)-Hv)/norm(Hv);
            end
            if obj.verbose
                disp('Hessian finite difference check')
                for k = 1:p
                    disp(['h = ',num2str(h(k)),' and error = ',num2str(diffs(k))])
                end
                disp(' ')
            end
        end
        
    end
end

