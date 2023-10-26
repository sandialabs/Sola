classdef HDSA_MD_Continuation_Interface < handle
    
    properties
        md_interface;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        [u] = State_Solve(this,z);
        
        [u_out] = Apply_Solution_Operator_Jacobian(this,z_in,z);
        
    end
    
    methods
        function this = HDSA_MD_Continuation_Interface(md_interface)
            this.md_interface = md_interface;
        end
        
        
    end
    
end