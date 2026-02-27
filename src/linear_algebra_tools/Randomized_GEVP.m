classdef Randomized_GEVP < handle

    properties
        vec
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [vec_out] = Apply_Operator(this, vec_in)

        [vec_out] = Apply_Weighting_Operator(this, vec_in)

        [vec_out] = Apply_Weighting_Operator_Inverse(this, vec_in)

        [samples] = Generate_Random_Samples(this, num_samples)

    end

    methods

        function this = Randomized_GEVP(vec)
            this.vec = 0 * vec;
        end

        function [evecs, evals] = Compute_GEVP(this, num_evals, oversampling)

            kpp = num_evals + oversampling;
            tmp1 = this.Generate_Random_Samples(kpp);
            tmp2 = this.Apply_Operator(tmp1);
            Y = this.Apply_Weighting_Operator_Inverse(tmp2);

            Q = this.CholQR(Y, 'weighting');

            AQ = this.Apply_Operator(Q);
            T = Q' * AQ;

            T = 0.5 * (T + T');

            try
                R_T = chol(T);
                M = AQ * linsolve(R_T, eye(size(R_T, 1)));
                [~, WQ, R] = this.CholQR(M, 'weighting_inverse');
                [U_M, Sigma_M, ~] = svd(R);

                scale = sign(U_M(1, :));
                U_M = U_M .* (ones(size(U_M, 1), 1) * scale);

                evecs = WQ * U_M(:, 1:num_evals);
                evals = diag(Sigma_M(1:num_evals, 1:num_evals)).^2;
            catch
                [S, Lambda] = eig(T, 'vector');
                [~, I] = sort(Lambda, 'descend');
                S = S(:, I);
                Lambda = Lambda(I);
                evecs = Q * S(:, 1:num_evals);
                evals = Lambda(1:num_evals);
            end
        end

        function [Q, WQ, R] = CholQR(this, Z, type)

            R_Z = chol(Z' * Z);
            Q_Z = Z * linsolve(R_Z, eye(size(R_Z, 1)));
            % The commented line below is superior to the two lines above.
            % However, the approach above is preferable for parallel
            % implmentations and hence is there for comparison.
            % [Q_Z, R_Z] = qr(Z, "econ");

            if strcmp(type, 'weighting')
                W_Q_Z = this.Apply_Weighting_Operator(Q_Z);
            elseif strcmp(type, 'weighting_inverse')
                W_Q_Z = this.Apply_Weighting_Operator_Inverse(Q_Z);
            else
                disp('Error specifying type in CholQR');
            end
            C = Q_Z' * W_Q_Z;
            R_C = chol(C);
            R_C_inv = linsolve(R_C, eye(size(R_C, 1)));
            Q = Q_Z * R_C_inv;
            WQ = W_Q_Z * R_C_inv;
            R = R_C * R_Z;
        end

    end

end
