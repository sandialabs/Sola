classdef Adv_Diff_Constraint < Dynamic_Constraint


    properties
        A;
        M;
        x;
    end

    methods (Access = public)
        
        function [f, f_y, f_z] = Time_Instance_RHS(this,y,z,t)
            % Extract the control for the given time.
            w = this.Temporal_Weights(t);
            zt = reshape(z,this.m,this.N)*w;
            
            % Evaluate the RHS and its derivatives.
            f = linsolve(this.M,-this.A*y + this.M*zt);
            f_y = linsolve(this.M,-this.A);
            f_z = kron(w',eye(this.m));
        end
        
        function [h, h_z] = Initial_Condition(this,z)
            h = zeros(this.m,1);
            h_z = zeros(this.m,this.n);
        end
        
        function [Mv] = Time_Instance_RHS_yy_Apply(this,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(this.m,num_vecs);
        end
        
        function [Mv] = Time_Instance_RHS_yz_Apply(this,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(this.m,num_vecs);
        end
        
        function [Mv] = Time_Instance_RHS_zy_Apply(this,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(this.n,num_vecs);
        end
        
        function [Mv] = Time_Instance_RHS_zz_Apply(this,v,y,z,t,lambda)
            num_vecs = size(v,2);
            Mv = zeros(this.n,num_vecs);
        end
        
        function [Mv] = Initial_Condition_zz_Apply(this,v,z,lambda)
            num_vecs = size(v,2);
            Mv = zeros(this.n,num_vecs);
        end
        
    end
    
    methods (Access = public)
        function this = Adv_Diff_Constraint(m,n,T,N)
            this = this@Dynamic_Constraint(m,n,T,N);
            
            Pe = 1; % Peclet number

            % Spatial domain
            this.x = linspace(0,1,m)';
            dx = this.x(2)-this.x(1);

            % Mass matrix
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*dx*M;

            % Stiffness matrix (diffusion)
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/dx)*S;

            % Viscosity matrix (advection)
            V = diag(0*ones(1,m)) + (1/2)*diag(ones(1,m-1),1) + (-1/2)*diag(ones(1,m-1),-1);
            V(1,1) = -1/2;
            V(end,end) = 1/2;

            % Discretized PDE: M f(y,z) = -Ay + Mz, A = S + Pe V.
            A = S + Pe*V;

            % Store properties.
            this.A = A;
            this.M = M;
        end
        
        function [w] = Temporal_Weights(this,t)
            w = (this.t_mesh-t)/(this.t_mesh(2)-this.t_mesh(1));
            Im = intersect(find(w<=0),find(abs(w)<=1));
            Ip = intersect(find(w>0),find(abs(w)<=1));
            I = find(abs(w)>1);
            w(I) = 0;
            w(Im) = 1+w(Im);
            w(Ip) = 1-w(Ip);
        end

    end

end
