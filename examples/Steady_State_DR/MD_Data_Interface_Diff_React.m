classdef MD_Data_Interface_Diff_React < MD_Data_Interface

    properties
        ensemble_id_k
        ensemble_id_i
        design_type
    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Optimization_Results.mat').u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat').z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            if isempty(this.ensemble_id_i)
                Z = load('Optimization_Results.mat').Z;
            elseif strcmp(this.design_type, 'OED')
                Z = load('OED_Ensemble_Results.mat', 'oed_Z_samps').oed_Z_samps{this.ensemble_id_k, this.ensemble_id_i};
            elseif strcmp(this.design_type, 'Random')
                Z = load('OED_Ensemble_Results.mat', 'rand_Z_samps').rand_Z_samps{this.ensemble_id_k, this.ensemble_id_i};
            elseif strcmp(this.design_type, 'SubRandom')
                Z = load('OED_Ensemble_Results.mat', 'subrand_Z_samps').subrand_Z_samps{this.ensemble_id_k, this.ensemble_id_i};
            end
        end

        function [D] = Load_d_Data(this)
            if isempty(this.ensemble_id_i)
                D = load('Optimization_Results.mat').D;
            elseif strcmp(this.design_type, 'OED')
                D = load('OED_Ensemble_Results.mat', 'oed_D_samps').oed_D_samps{this.ensemble_id_k, this.ensemble_id_i};
            elseif strcmp(this.design_type, 'Random')
                D = load('OED_Ensemble_Results.mat', 'rand_D_samps').rand_D_samps{this.ensemble_id_k, this.ensemble_id_i};
            elseif strcmp(this.design_type, 'SubRandom')
                D = load('OED_Ensemble_Results.mat', 'subrand_D_samps').subrand_D_samps{this.ensemble_id_k, this.ensemble_id_i};
            end
        end

        function this = MD_Data_Interface_Diff_React(varargin)
            switch nargin
                case 0
                    % do nothing
                case 2
                    [this.u_init, this.z_init] = deal(varargin{:});
                case 3
                    [this.ensemble_id_k, this.ensemble_id_i, this.design_type] = deal(varargin{:});
                otherwise
                    error("Please enter the correct number of inputs into the data interface.");
            end

            if isempty(this.z_opt)
                this.z_opt = this.z_init;
                this.u_opt = this.u_init;
            end

            if isempty(this.z_init)
                this.z_init = this.z_opt;
                this.u_init = this.u_opt;
            end
        end

    end

end
