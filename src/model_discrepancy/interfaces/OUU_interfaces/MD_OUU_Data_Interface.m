classdef MD_OUU_Data_Interface < MD_Data_Interface

    properties
        n_u
        n_r
    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [us_opt] = Load_Optimal_us(this, s)

        [z_opt] = Load_Optimal_z(this)

    end

    methods

        function [u_mat] = Reshape_State_to_Mat(this, u_vec)
            u_mat = reshape(u_vec, this.n_r, this.n_u)';
        end

        function [u_vec] = Reshape_State_to_Vec(this, u_mat)
            u_mat = u_mat';
            u_vec = u_mat(:);
        end

    end

    %% Virtual functions for user implementation
    methods

        % Defaults to return empty array
        % Overload function to load high-fidelity data
        function [Z] = Load_Z_Data(this)
            Z = [];
            disp('Load_Z_Data must be implemented for any analyses except for optimal experimental design');
        end

        % Defaults to return empty array
        % Overload function to load high-fidelity data
        function [Ds] = Load_ds_Data(this, s)
            Ds = [];
            disp('Load_ds_Data must be implemented for any analyses except for optimal experimental design');
        end

        % Defaults to return all state elements
        % Overload function to extract component i from the state
        % Returns a vector of integers I index elements of component i
        function [I] = Separate_State_Components_Per_Sample(this, i)
            I = 1:this.n_u;
        end

        %%
        function [u_opt] = Load_Optimal_u(this)
            u1 = this.Load_Optimal_us(1);
            this.n_u = length(u1);
            u_opt = zeros(this.n_u, this.n_r);
            u_opt(:, 1) = u1;
            for s = 2:this.n_r
                u_opt(:, s) = this.Load_Optimal_us(s);
            end
            u_opt = this.Reshape_State_to_Vec(u_opt);
        end

        function [D] = Load_d_Data(this)
            D1 = this.Load_ds_Data(1);
            N = size(D1, 2);
            Dtmp = zeros(this.n_u, this.n_r, N);
            Dtmp(:, 1, :) = D1;
            for s = 2:this.n_r
                Ds = this.Load_ds_Data(s);
                Dtmp(:, s, :) = Ds;
            end
            D = zeros(this.n_u * this.n_r, N);
            for k = 1:N
                D(:, k) = this.Reshape_State_to_Vec(Dtmp(:, :, k));
            end
        end

        function [I] = Separate_State_Components(this, i)
            I_sample = this.Separate_State_Components_Per_Sample(i);
            m = length(I_sample);
            I = zeros(m, this.n_r);
            for s = 1:this.n_r
                I(:, s) = I_sample + this.n_u * (s - 1);
            end
            I =  this.Reshape_State_to_Vec(I);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_OUU_Data_Interface(n_r)
            this.n_r = n_r;
        end

    end

end
