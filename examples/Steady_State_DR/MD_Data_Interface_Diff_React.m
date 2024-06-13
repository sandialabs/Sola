classdef MD_Data_Interface_Diff_React < MD_Data_Interface

    properties
        ensemble_id_k
        ensemble_id_i
        design_type
    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            if isempty(this.u_opt)
                u_opt = load('Optimization_Results.mat').u_lofi;
            else
                u_opt = this.u_opt;
            end
        end

        function [z_opt] = Load_Optimal_z(this)
            if isempty(this.z_opt)
                z_opt = load('Optimization_Results.mat').z_lofi;
            else
                z_opt = this.z_opt;
            end
        end

        function [Z] = Load_Z_Data(this)
            if ~isempty(this.Z)
                Z = this.Z;
            elseif isempty(this.ensemble_id_i)
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
            if ~isempty(this.D)
                D = this.D;
            elseif isempty(this.ensemble_id_i)
                D = load('Optimization_Results.mat').D;
            elseif strcmp(this.design_type, 'OED')
                D = load('OED_Ensemble_Results.mat', 'oed_D_samps').oed_D_samps{this.ensemble_id_k, this.ensemble_id_i};
            elseif strcmp(this.design_type, 'Random')
                D = load('OED_Ensemble_Results.mat', 'rand_D_samps').rand_D_samps{this.ensemble_id_k, this.ensemble_id_i};
            elseif strcmp(this.design_type, 'SubRandom')
                D = load('OED_Ensemble_Results.mat', 'subrand_D_samps').subrand_D_samps{this.ensemble_id_k, this.ensemble_id_i};
            end
        end

        function Set_Z_and_D(this, Z, D)
            this.Z = Z;
            this.D = D;
        end

        function this = MD_Data_Interface_Diff_React(varargin)
            % This Operator is overloaded to allow for (u_opt, z_opt); (u_opt, z_opt, Z, D); or (e_k, e_i, type)
            % This is temporarily done to preserve backward-compatibility with OED Drivers
            % Also, automatically loads data so it doesn't have to be repeated in MD_Data_Interface
            switch nargin
                case 2
                    this.u_opt = varargin{1};
                    this.z_opt = varargin{2};
                case 3
                    this.ensemble_id_k = varargin{1};
                    this.ensemble_id_i = varargin{2};
                    this.design_type = varargin{3};
                case 4
                    this.u_opt = varargin{1};
                    this.z_opt = varargin{2};
                    this.Z = varargin{3};
                    this.D = varargin{4};
                otherwise
                    error("Please enter the correct number of inputs into MD_Data_Interface_Diff_React.");
            end
        end

    end

end
