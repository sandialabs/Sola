classdef Basis < handle
    % Dimension reduction for mapping between :math:`n_y`-dimensional
    % (full) states and :math:`n_y'`-dimensional (reduced) states.

    properties (Abstract)
        n_y                 % Ambient full dimension :math:`n_y`.
        r                   % Latent reduced dimension :math:`n_y'`.
    end

    methods (Access = public)

        function [states_compressed] = Compress(this, states)
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
            error('Compress() not implemented');
        end

        function [states] = Decompress(this, states_compressed)
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
            error('Decompress() not implemented');
        end

    end

end
