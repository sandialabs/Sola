classdef HDSA_Sabl_MD_Continuation_Interface < HDSA_MD_Continuation_Interface
    
    properties
        con;
        z_current;
        u_current;
        m;
    end
    
    methods (Access = public)
        
        function [u] = State_Solve(this,z)
            if norm(z-this.z_current)==0
                u = this.u_current;
            else
                u = this.con.State_Solve(z);
                this.u_current = u;
                this.z_current = z;
            end
        end
        
        function [u_out] = Apply_Solution_Operator_Jacobian(this,z_in,z)
            u = this.State_Solve(z);
            tmp = this.con.c_z_Apply(z_in,u,z);
            u_out = -this.con.c_u_Inverse_Apply(tmp,u,z);
        end
        
    end
    
    methods
        function this = HDSA_Sabl_MD_Continuation_Interface(md_interface,con)
            this@HDSA_MD_Continuation_Interface(md_interface);
            this.con = con;
            this.z_current = md_interface.Load_Optimal_z();
            this.u_current = md_interface.Load_Optimal_u();
            this.m = length(this.u_current);
        end
        
        
    end
    
end