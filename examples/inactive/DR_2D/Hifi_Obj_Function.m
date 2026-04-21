%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Hifi_Obj_Function < handle

    properties
        diff_react_hifi
        obj
    end

    methods

        function val = Jhat(this, z)
            u = this.diff_react_hifi.State_Solve(this.diff_react_hifi.Map_z_to_Control_Fun(z));
            val = this.obj.J(u, z);
        end

        function this = Hifi_Obj_Function(diff_react_hifi, obj)
            this.diff_react_hifi = diff_react_hifi;
            this.obj = obj;
        end

    end
end
