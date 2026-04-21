%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Basis < handle
    % Dimension reduction for mapping between :math:`n_y`-dimensional
    % (full) states and :math:`n_y'`-dimensional (reduced) states.

    properties (Abstract)
        n_y                 % Ambient full dimension :math:`n_y`.
        r                   % Latent reduced dimension :math:`n_y'`.
    end

    methods (Abstract, Access = public)

        [states_compressed] = Compress(this, states)
        % Compress full-order states to the reduced-order space.
        %
        % Parameters
        % ----------
        % states
        %   State(s) :math:`\y\in\R^{n_y}` to compress.
        %
        % Returns
        % -------
        % states_compressed : vector(s)
        %   Compressed state(s) :math:`\hat{\y}\in\R^{r}`.

        [states] = Decompress(this, states_compressed)
        % Decompress reduced-order states to the full-order space.
        %
        % Parameters
        % ----------
        % states_compressed
        %   Compressed state(s) :math:`\hat{\y}\in\R^{r}` to decompress.
        %
        % Returns
        % -------
        % states : vector(s)
        %   Decompressed state(s) :math:`\y\in\R^{n_y}`.

    end

    methods (Access = public)

        function [states_projected] = Project(this, states)
            % Project full-order states to the span of the basis.
            %
            % This is compression followed by decompression:
            % ``states_projected = Decompress(Compress(states))``.
            %
            % Parameters
            % ----------
            % states
            %   State(s) :math:`\y\in\R^{n_y}` to project.
            %
            % Returns
            % -------
            % states_projected : vector(s)
            %   Projected state(s) :math:`\tilde{\y}\in\R^{n_y}`.
            states_compressed = this.Compress(states);
            states_projected = this.Decompress(states_compressed);
        end

    end

end
