classdef Mass_Spring_LoFi < Constrained_ODE_Optimization

    properties
        mass_spring_coupled;
        m1;
        m2;
        k1;
        k2;
        k3;
        f2;
        P_z;
    end

    methods

        function obj = Mass_Spring_LoFi(mass_spring_coupled)
            obj = obj@Constrained_ODE_Optimization(2,mass_spring_coupled.n,mass_spring_coupled.T,mass_spring_coupled.N);
            obj.mass_spring_coupled = mass_spring_coupled;

            obj.m1 = mass_spring_coupled.m1;
            obj.m2 = mass_spring_coupled.m2;
            obj.k1 = mass_spring_coupled.k1;
            obj.k2 = mass_spring_coupled.k2;
            obj.k3 = mass_spring_coupled.k3;
            obj.f2 = mass_spring_coupled.f2;
            obj.P_z = mass_spring_coupled.P_z;

        end

        function [val, grad_y] = Time_Instance_Objective(obj,y,t)
            [val,grad_y] = obj.mass_spring_coupled.Time_Instance_Objective(y,t);
            grad_y = grad_y(1:2);
        end

        function [val,grad_z] = Regularization_Objective(obj,z)
            [val,grad_z] = obj.mass_spring_coupled.Regularization_Objective(z);
        end

        % ODE system with four states y=(x_1,v_1)
        % x_1' = v_1
        % v_1' = (1/m_1)*( k_2*x_2 - (k_1+k_2)*x_1 + f_1(z) )
        function [f, f_y, f_z] = Time_Instance_RHS(obj,y,z,t)
            x1 = y(1);
            v1 = y(2);

            coeffs = obj.mass_spring_coupled.Temporal_Weights(t);
            f1 = (obj.P_z*z)'*coeffs;

            f = zeros(2,1);
            f(1) = v1;
            f(2) = (1/obj.m1)*( - (obj.k1+obj.k2)*x1 + f1 );

            f_y = zeros(2,2);
            f_y(1,2) = 1;
            f_y(2,1) = -(1/obj.m1)*(obj.k1+obj.k2);

            f_z = zeros(2,size(obj.P_z,2));
            f_z(2,:) = (1/obj.m1)*obj.P_z'*coeffs;
        end

        function [h, h_z] = Initial_Condition(obj,z)
            h = zeros(2,1);
            h_z = zeros(2,size(obj.P_z,2));
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(obj,v,y,t)
            Mv = zeros(size(v));
            Mv(1,:) = v(1,:);
        end

        function [Mv] = Regularization_Objective_zz_Apply(obj,v,z)
            Mv = obj.mass_spring_coupled.Regularization_Objective_zz_Apply(v,z);
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(obj,v,y,z,t,lambda)
            Mv = 0*v;
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(obj,v,y,z,t,lambda)
            Mv = zeros(2,size(v,2));
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(obj,v,y,z,t,lambda)
            Mv = zeros(length(z),size(v,2));
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(obj,v,y,z,t,lambda)
            Mv = 0*v;
        end

        function [Mv] = Initial_Condition_zz_Apply(obj,v,z,lambda)
            Mv = 0*v;
        end

    end



end
