classdef Thermal_Constraint < Dynamic_Constraint
    
    
    properties
        x;
        M;
        S;
        forcing;
    end
    
    methods (Access = public)

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            D = this.Assembly(z);
            f_y = -linsolve(this.M,D);
            f = f_y*y + (10^3)*this.forcing(this.x,t);
            f_z = -linsolve(this.M,this.Assembly_z_Jacobian(y));
        end

        function [h, h_z] = Initial_Condition(this, z)
           h = ones(this.m,1);
           h_z = zeros(this.m,this.m);
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
           Mv = 0*v; 
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(this.m,size(v,2));
            for k = 1:size(v,2)
                D = this.Assembly(z);
                D_pert = this.Assembly(z+v(:,k));
                Mv(:,k) = -(D_pert'- D')*linsolve(this.M,lambda);
            end
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            Mv = zeros(this.m,size(v,2));
            for k = 1:size(v,2)
                D_diff = this.Assembly_z_Jacobian(y);
                D_diff_pert = this.Assembly_z_Jacobian(y+v(:,k));
                Mv(:,k) = -(D_diff_pert' - D_diff')*linsolve(this.M,lambda);
            end
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
           Mv = 0*v; 
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
           Mv = 0*v; 
        end

        
    end
    
    methods (Access = public)
        
        function [val] = Forcing_Function(this,x,t)
            val = this.forcing(x,t);
        end
        
        function [kappa] = Diffusion_Coeff(this,x,z)
            kappa = interp1(this.x,z,x);
        end
        
        function [D] = Assembly(this,z)
            h = this.x(2)-this.x(1);
            x1 = (0:(this.m-2))*h + (h/2)*(-1/sqrt(3) + 1);
            x2 = (0:(this.m-2))*h + (h/2)*(1/sqrt(3) + 1);
            diff_x = this.Diffusion_Coeff([x1;x2],z);
            
            s = sum(diff_x,1);
            D = diag(([0,s] + [s,0])*(1/h)/2) + (-1)*diag(s,1)*(1/h)/2 + (-1)*diag(s,-1)*(1/h)/2;
            
            %             phi_down_dot = -[1;1]/h;
            %             phi_up_dot = [1;1]/h;
            %             diff_x1 = diff_x(1:(this.m-1));
            %             diff_x2 = diff_x(this.m:end);
            %             D = zeros(this.m,this.m);
            %             for i = 1:this.m
            %                 if i > 1
            %                     perm = [diff_x1(i-1) ; diff_x2(i-1)];
            %                     D(i-1,i) = (h/2)*sum(phi_up_dot.*phi_down_dot.*perm);
            %                     D(i,i) = (h/2)*sum(phi_up_dot.*phi_up_dot.*perm);
            %                 end
            %                 if i < this.m
            %                     perm = [diff_x1(i) ; diff_x2(i)];
            %                     D(i,i) = D(i,i) + (h/2)*sum(phi_down_dot.*phi_down_dot.*perm);
            %                     D(i+1,i) = (h/2)*sum(phi_up_dot.*phi_down_dot.*perm);
            %                 end
            %             end
        end

        function [D_diff] = Assembly_z_Jacobian(this,u)
            h = this.x(2)-this.x(1);
            up = (u(2:end)-u(1:end-1))/h;
            d = [ -(1/2)*up(1) ; (1/2)*up(1:end-1) - (1/2)*up(2:end) ; (1/2)*up(end) ];
            D_diff = diag(d) + (-1/2)*diag(up,1) + (1/2)*diag(up,-1);

            %             D_diff = zeros(this.m,this.m);
            %             x1 = (h/2)*(-1/sqrt(3) + 1);
            %             x2 = (h/2)*(1/sqrt(3) + 1);
            %             phi_down = [x2;x1]/h;
            %             phi_up = [x1;x2]/h;
            %             phi_down_dot = -[1;1]/h;
            %             phi_up_dot = [1;1]/h;
            %             for i = 1:this.m
            %                 if i > 1
            %                     u_prime = u(i-1)*phi_down_dot + u(i)*phi_up_dot;
            %                     D_diff(i-1,i) = (h/2)*sum(phi_up.*phi_down_dot.*u_prime);
            %                     D_diff(i,i) = (h/2)*sum(phi_up.*phi_up_dot.*u_prime);
            %                 end
            %                 if i < this.m
            %                     u_prime = u(i)*phi_down_dot + u(i+1)*phi_up_dot;
            %                     D_diff(i,i) = D_diff(i,i) + (h/2)*sum(phi_down.*phi_down_dot.*u_prime);
            %                     D_diff(i+1,i) = (h/2)*sum(phi_down.*phi_up_dot.*u_prime);
            %                 end
            %             end
        end
        
        function this = Thermal_Constraint(m, n, T, N)
            this = this@Dynamic_Constraint(m, n, T, N);

            this.x = linspace(0,1,m)';
            
            h = this.x(2)-this.x(1);
            
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            this.M = M;
            
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            this.S = S;
          
            this.forcing = @(x,t) this.x.^2;
            
        end

    end
end

