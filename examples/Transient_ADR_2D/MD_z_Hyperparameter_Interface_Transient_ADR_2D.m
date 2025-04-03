classdef MD_z_Hyperparameter_Interface_Transient_ADR_2D < MD_z_Hyperparameter_Interface

    properties
        opt
        basis1
        basis2
        t
    end

    methods (Access = public)

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function [u] = State_Solve(this, z)
            tmp = this.opt.con.State_Solve(z(:,1));
            u1 = this.Map_Reduced_to_Full(tmp);
            u = zeros(length(u1),size(z,2));
            u(:,1) = u1;
            for k = 2:size(z,2)
               tmp = this.opt.con.State_Solve(z(:,k));
               u(:,k) = this.Map_Reduced_to_Full(tmp);
            end
        end

        function [u_full] = Map_Reduced_to_Full(this,u_red)
            u_reshape = reshape(u_red, this.opt.con.n_1+this.opt.con.n_2, this.opt.con.n_t);
            Y_rom_1 = this.basis1.Decompress(u_reshape(1:this.opt.con.n_1, :));
            Y_rom_2 = this.basis2.Decompress(u_reshape(this.opt.con.n_1 + 1:end, :));
            Y_rom = [Y_rom_1; Y_rom_2];
            u_full = Y_rom(:);
        end

        function this = MD_z_Hyperparameter_Interface_Transient_ADR_2D(num_state_solves,opt,basis1,basis2)
            this@MD_z_Hyperparameter_Interface('transient vector',num_state_solves);
            this.opt = opt;
            this.basis1 = basis1;
            this.basis2 = basis2;
            this.t = opt.con.t_mesh(1:(end-1));
        end

    end

end
