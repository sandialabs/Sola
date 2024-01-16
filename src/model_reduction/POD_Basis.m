classdef POD_Basis < Basis
    % Dimension reduction via proper orthogonal decomposition (POD).
    %
    % For a collection of state measurements
    %
    % .. math:: \Y = [~\y(t_1)~~\cdots~~\y(t_{n_t})~]\in\R^{n_y \times n_t},
    %
    % the POD approximation for the state :math:`\y` is given by
    %
    % .. math:: \y = \bar{\y} + \V\hat{\y} = \sum_{i=1}^{r}\hat{y}_i \v_i
    %
    % where :math:`\V` consists of the first :math:`r` left singular vectors of
    % :math:`\Y - \bar{\y}\1\trp`,
    %
    % .. math:: \V = [~\v_1~~\cdots~~\v_r~] = \boldsymbol{\Phi}_{:,1:r},
    %  \qquad \texttt{svd}(\Y - \bar{\y}\1\trp) = \boldsymbol{\Phi\Sigma\Psi}\trp.
    %
    % The reference snapshot :math:`\bar{\y}` is either the mean snapshot
    % :math:`\bar{\y} = \frac{1}{n_t}\sum_{j=1}^{n_t}\y_j` or zero (default).

    properties
        r                   % Number of basis vectors :math:`r`.
        economize           % If ``true``, store only the first ``r`` singular vectors.
    end

    properties (SetAccess = protected)
        ybar
        maxdim
        singular_vectors
        singular_values
        W
    end

    properties (Dependent)
        n_y                 % Dimension :math:`n_y` of the state :math:`\y` to approximate.
        V                   % Basis matrix :math:`\V\in\R^{n_y \times r}`.
    end

    methods

        %% Constructor

        function this = POD_Basis(Y, shift, W)
            % Parameters
            % ----------
            % Y
            %   Snapshot matrix
            %   :math:`\Y = [~\y(t_1)~~\cdots~~\y(t_{n_t})~]\in\R^{m\times k}`.
            % shift
            %   (Optional) if ``true``, set the reference snapshot to the mean
            %   :math:`\bar{\y} = \frac{1}{n_t}\sum_{j=1}^{n_t}\y_j`.
            %   If ``false`` (default), use :math:`\bar{\y} = \0`.
            % W
            %   (Optional) Weight matrix :math:`\W\in\R^{m\times m}`.
            %   If provided, the basis :math:`\V` is weighted to satisfy
            %   :math:`\V\trp\W\V = \I_{r}`.
            %   If provided as a vector :math:`\w\in\R^{m\times 1}`, use
            %   :math:`\W = \text{diag}(\w)`.
            arguments
                Y (:, :) double
                shift = false
                W (:, :) {mustBeNumeric} = []
            end

            this.economize = false;
            this.maxdim = size(Y, 2);       % Maximum number of basis vectors.
            this.r = this.maxdim;

            % Shift snapshots by the mean if desired.
            if shift
                this.ybar = mean(Y, 2);
                Y = Y - this.ybar;
            else
                this.ybar = zeros(size(Y, 1), 1);
            end

            % Weighting.
            if size(W, 1) == 0
                W = ones(size(Y, 1), 1);
            end
            if size(W, 1) == 1 || size(W, 2) == 1
                this.W = diag(W);
                Wsqrt = diag(sqrt(W));
                Winvsqrt = diag(1 ./ sqrt(W));
            else
                this.W = W;
                Wsqrt = sqrtm(W);
                Winvsqrt = linsolve(Wsqrt, eye(size(W, 1)));
            end

            % Do the SVD and store the singular values / vectors.
            [Phi, svdvals, ~] = svd(Wsqrt * Y, "econ");
            this.singular_vectors = Winvsqrt * Phi;
            this.singular_values = diag(svdvals);
        end

        %% Getters and setters.

        function [state_dimension] = get.n_y(this)
            % Dimension :math:`n_y` of the state :math:`\y` being approximated.
            state_dimension = size(this.singular_vectors, 1);
        end

        function [basis] = get.V(this)
            % Basis matrix :math:`\V\in\R^{n_y \times r}`.
            basis = this.singular_vectors(:, 1:this.r);
        end

        function set.r(this, r)
            % Set the reduced dimension.
            %
            % Parameters
            % ----------
            % r
            %   Reduced dimension, an integer between 1 and :math:`\min\{n_t, n_y\}`.

            if r > this.maxdim
                error('invalid reduced dimension r');
            end

            % Economize if desired.
            if this.economize && r < this.r
                this.singular_vectors = this.singular_vectors(:, 1:r);
                this.maxdim = r;
            end
            this.r = r;
        end

        function [dimension] = Set_Reduced_Dimension_From_Cumulative_Energy(this, energy)
            % Set the ``r`` property based on the cumulative energy criteria:
            % choose :math:`r` to be the smallest integer such that
            %
            % .. math:: \frac{\sum_{i=1}^{r}\sigma_i^2}{\sum_{j=1}^{n_t}\sigma_j^2} \ge \kappa
            %
            % for some user-specified energy level :math:`\kappa`.
            %
            % Parameters
            % ----------
            % energy
            %   Cumulative energy level :math:`\kappa`.
            %
            % Returns
            % -------
            % dimension : uint8
            %   Reduced dimension :math:`r`.
            svdvals2 = this.singular_values.^2;
            cumulative_energies = cumsum(svdvals2) / sum(svdvals2);
            dimension = sum(cumulative_energies < energy) + 1;
            this.r = dimension;
        end

        function [dimension] = Set_Reduced_Dimension_From_Residual_Energy(this, energy)
            % Set the ``r`` property based on the residual energy criteria:
            % choose :math:`r` to be the smallest integer such that
            %
            % .. math:: 1 - \frac{\sum_{i=1}^{r}\sigma_i^2}{\sum_{j=1}^{n_t}\sigma_j^2} < \epsilon
            %
            % for some user-specified tolerance :math:`\epsilon`.
            %
            % Parameters
            % ----------
            % energy
            %   Residual energy tolerance :math:`\epsilon`.
            %
            % Returns
            % -------
            % dimension : uint8
            %   Reduced dimension :math:`r`.
            dimension = this.Set_Reduced_Dimension_From_Cumulative_Energy(1 - energy);
            % svdvals2 = this.singular_values .^ 2;
            % residual_energies = 1 - (cumsum(svdvals) / sum(svdvals));
            % dimension = sum(residual_energies > energy) + 1;
            % this.r = dimension;
        end

        %% Dimensionality reduction: compression and decompression

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
            states_compressed = (this.V' * this.W) * (states - this.ybar);
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
            states = (this.V * states_compressed) + this.ybar;
        end

    end
end
