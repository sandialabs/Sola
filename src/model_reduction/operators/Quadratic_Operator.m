classdef Quadratic_Operator < OpInf_Operator
    % Quadratic state operator
    % :math:`\mathcal{H}(\y,\q) = \H[\q\otimes\q],~~\H\in\R^{n_y \times n_y^2}`.
    %
    % Internally, this is represented efficiently with a matrix
    % :math:`\tilde{\H}\in\R^{n_y \times (n_y(n_y + 1)/2)}` such that
    %
    % .. math:: \H[\y\otimes \y] = \tilde{\H}[\y\tilde{\otimes}\y]
    %
    % for all :math:`\y\in\R^{n_y}`.

    properties (Access = protected)
        mask
        prejac
    end

    methods (Access = public)

        function this = Quadratic_Operator(H)
            arguments
                H (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator(H);
        end

        function Set_Entries(this, H)
            % Set the entries of the operator.
            %
            % Parameters
            % ----------
            % H
            %   Operator entries, either the full :math:`\H\in\R^{n_y \times n_y^2}`
            %   or the compressed :math:`\tilde{\H}\in\R^{n_y \times (n_y (n_y + 1)/2)}`.

            % Extract dimensions.
            n_y = size(H, 1);
            r2 = Quadratic_Operator.Column_Dimension(n_y);
            ncols = size(H, 2);

            % Get the compressed Kronecker mask.
            msk = zeros(r2, 2);
            count = 1;
            for i = 1:n_y
                for j = 1:i
                    msk(count, :) = [i j];
                    count = count + 1;
                end
            end
            this.mask = msk;

            % Get both the full and compressed representations.
            if ncols == n_y^2
                Hc = this.compress_quadratic(H);
            elseif ncols == r2
                Hc = H;
                H = this.expand_quadratic(Hc);
            end

            % Compute the pre-Jacobian tensor.
            Ht = reshape(H, n_y, n_y, n_y);
            this.prejac = reshape(Ht + permute(Ht, [1, 3, 2]), n_y, n_y^2);

            Set_Entries@OpInf_Operator(this, Hc);
        end

        function [out] = Apply(this, y, ~)
            % Apply the operator to the given state(s) and input(s):
            % :math:`\mathcal{H}(\y,\q) = \H[\y \otimes \y]` or
            % :math:`\mathcal{H}(\Y,\Q) = \H[\Y \odot \Y]`
            % where :math:`\odot` is the Khatri--Rao product, i.e.,
            % the column-wise Kronecker product.
            %
            % Parameters
            % ----------
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}`,
            %   or an :math:`n_y \times n_t` collection of these.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`
            %   or an :math:`n_q \times n_t` collection of these.
            %
            % Returns
            % -------
            % out : :math:`n_y`-vector or :math:`n_y\times n_t` matrix
            %   Application of the operator, :math:`\H[\y\otimes\y]`
            %   or :math:`\H[\Y\odot\Y]`.
            out = this.entries * (y(this.mask(:, 1), :) .* y(this.mask(:, 2), :));
        end

        function [jac] = Jacobian_y(this, y, ~)
            % Construct the state Jacobian of the operator.
            %
            % Parameters
            % ----------
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            %
            % Returns
            % -------
            % jac : :math:`n_y \times n_y` matrixr
            %   State Jacobian
            %   :math:`\H[(\I\otimes\y) + (\y\otimes\I)]\in\R^{n_y \times n_y}`.
            jac = this.prejac * kron(eye(this.n_y), y);
        end

        function [Mv] = Hessian_yy_Apply(this, v, ~, ~, lambda)
            % Compute the action of the :math:`y,y` Hessian of the operator.
            %
            % Parameters
            % ----------
            % v
            %   Search direction :math:`\v\in\R^{n_y}`.
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            % lambda
            %   Adjoint state :math:`\lambda(t)\in\R^{n_y}` at time :math:`t`.
            %
            % Returns
            % -------
            % Mv : :math:`n_y`-vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\mathcal{H}_{y,y}(\y,\q)\v
            %   = \mathcal{H}_{y}(\v,\q)\trp\bflambda`.
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y, num_vecs);
            for j = 1:num_vecs
                Mv(:, j) = this.Jacobian_y(v(:, j))' * lambda;
            end
        end

    end

    methods (Static, Access = public)

        function [d] = Column_Dimension(n_y, ~)
            % Column dimension of the operator entries.
            %
            % Parameters
            % ----------
            % n_y
            %   Dimension :math:`n_y` of the ODE state :math:`\y(t)`.
            % n_q
            %   Dimension :math:`n_q` of the input :math:`\q(t)`.
            %
            % Returns
            % -------
            % d : uint8
            %   Operator column dimension :math:`d = n_y (n_y + 1) / 2`.
            d = n_y * (n_y + 1) / 2;
        end

        function [block] = Datablock(Y, ~)
            % Construct the data matrix block corresponding to the operator.
            %
            % Since
            %
            % .. math::
            %  \min_{\H}\sum_{j=1}^{n_t}\left\| \H(\y_j, \q_j) \right\|_{2}^{2}
            %  = \min_{\H}\left\|
            %  \H[~\y_1\otimes\y_1~~\cdots~~\y_{n_t}\otimes\y_{n_t}~]
            %  \right\|_{F}^{2}
            %  = \min_{\H}\left\| \H[\Y\odot\Y]\right\|_{F}^{2},
            %
            % the data block would ordinarily be
            % :math:`\D = \Y\odot\Y
            % = [~\y_1\otimes\y_1~~\cdots~~\y_{n_t}\otimes\y_{n_t}~]
            % \in \R^{n_y^2 \times n_t}`.
            %
            % However, to preserve the symmetry of :math:`\H` and reduce computation,
            % a compressed Kronecker product :math:`\tilde{\otimes}` is used, so that
            % :math:`\D = \Y\tilde{\odot}\Y
            % = [~\y_1\tilde{\otimes}\y_1~~\cdots~~\y_{n_t}\tilde{\otimes}\y_{n_t}~]
            % \in \R^{(n_y(n_y + 1)/2) \times n_t}`.
            %
            % Parameters
            % ----------
            % Y
            %   State data :math:`\Y\in\R^{n_y \times n_t}`.
            %   Each column is a single state vector :math:`\y(t_j)`.
            % Q
            %   Input data :math:`\Q\in\R^{n_q \times n_t}`.
            %   Each column is a single input vector :math:`\q(t_j)`.
            %
            % Returns
            % -------
            % block : :math:`(n_y(n_y + 1)/2) \times n_t` matrix
            %     Data matrix block
            %     :math:`\D = \Y\tilde{\odot}\Y
            %     \in\R^{(n_y(n_y + 1)/2) \times n_t}`.
            block = Quadratic_Operator.kron_compressed(Y);
        end

    end

    methods (Static, Access = protected)

        function [y2] = kron_compressed(y)
            % Compressed Kronecker product.
            %
            % Compute the unique terms of the Kronecker product
            % :math:`\y \otimes \y` given by the convention
            %
            % .. math::
            %  \y \tilde{\otimes} \y = \left[\begin{array}{c}
            %  y_1 \y \\ y_2 \y_{:2} \\ y_3 \y_{:3} \\ \vdots
            %  \end{array}\right] \in \R^{n_y (n_y + 1) / 2}.
            %
            % Parameters
            % ----------
            % y
            %   Vector :math:`\y\in\R^{n_y}`
            %   to compute the compressed Kronecker product of,
            %   or a matrix :math:`\Y\in\R^{n_y \times n_t}` to compute
            %   the column-wise compressed Kronecker product of.
            %
            % Returns
            % -------
            % y2 : vector or matrix
            %   Compressed Kronecker product
            %   :math:`\y\tilde{\otimes}\y\in\R^{n_y (n_y + 1) / 2}` or
            %   :math:`\Y\tilde{\odot}\Y\in\R^{(n_y (n_y + 1)/2) \times n_t}`.
            n_y = size(y, 1);
            y2 = zeros(Quadratic_Operator.Column_Dimension(n_y), size(y, 2));
            index = 1;
            for i = 1:n_y
                newindex = index + i;
                y2(index:newindex - 1, :) = y(i, :) .* y(1:i, :);
                index = newindex;
            end
        end

        function [Hc] = compress_quadratic(H)
            % Convert the matricized quadratic tensor :math:`\H\in\R^{r \times r^2}`
            % to the compressed repesentation :math:`\tilde{\H}\in\R^{r \times (r(r+1)/2)}`
            % such that for all :math:`\y\in\R^{r}`,
            %
            % .. math:: \H[\y\otimes \y] = \tilde{\H}[\y\tilde{\otimes}\y]
            %
            % where :math:`\tilde{\otimes}` denotes the compressed Kronecker product
            % (see ``kron_compressed()``).
            %
            % Parameters
            % ----------
            % H
            %   Matricized quadratic tensor :math:`\H\in\R^{r \times r^2}`.
            %
            % Returns
            % -------
            % Hc : matrix
            %   Compressed matricized quadratic tensor
            %   :math:`\tilde{\H}\in\R^{r \times (r(r+1)/2)}`.

            % Extract dimensions.
            r = size(H, 1);
            r2 = size(H, 2);
            if r2 ~= r^2
                error('invalid size(H) for compress_quadratic(H), must be r x r^2');
            end
            s = Quadratic_Operator.Column_Dimension(r);

            % Fill the rows.
            Hc = zeros(r, s);
            fj = 1;
            for i = 0:(r - 1)
                for j = 0:i
                    if i == j
                        Hc(:, fj) = H(:, (i * r) + j + 1);
                    else
                        Hc(:, fj) = H(:, (i * r) + j + 1) + H(:, (j * r) + i + 1);
                    end
                    fj = fj + 1;
                end
            end
        end

        function [H] = expand_quadratic(Hc)
            % Convert the compressed matricized quadratic tensor
            % :math:`\tilde{\H}\in\R^{r \times (r(r+1)/2)}` to the full representation
            % :math:`\H\in\R^{r \times r^2}` such that for all :math:`\y\in\R^{r}`,
            %
            % .. math:: \H[\y\otimes \y] = \tilde{\H}[\y\tilde{\otimes}\y]
            %
            % where :math:`\tilde{\otimes}` denotes the compressed Kronecker product
            % (see ``kron_compressed()``).
            %
            % Parameters
            % ----------
            % Hc : matrix
            %   Compressed matricized quadratic tensor
            %   :math:`\tilde{\H}\in\R^{r \times (r(r+1)/2)}`.
            %
            % Returns
            % -------
            % H
            %   Matricized quadratic tensor :math:`\H\in\R^{r \times r^2}`.

            % Extract dimensions.
            r = size(Hc, 1);
            s = size(Hc, 2);
            if s ~= Quadratic_Operator.Column_Dimension(r)
                error('invalid size(Hc) for expand_quadratic(Hc), must be r x (r(r+1)/2)');
            end

            % Fill the rows.
            H = zeros(r, r^2);
            fj = 1;
            for i = 0:(r - 1)
                for j = 0:i
                    if i == j
                        H(:, (i * r) + j + 1) = Hc(:, fj);
                    else
                        fill = Hc(:, fj) / 2;
                        H(:, (i * r) + j + 1) = fill;
                        H(:, (j * r) + i + 1) = fill;
                    end
                    fj = fj + 1;
                end
            end
        end

    end
end
