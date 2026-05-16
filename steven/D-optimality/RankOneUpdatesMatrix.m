classdef RankOneUpdatesMatrix < handle

    properties
        M
        updates
        inv_updates
        denoms
        t
    end

    methods (Access = public)
        function this = RankOneUpdatesMatrix(M)
            this.M = M;
            this.updates = [];
            this.inv_updates = [];
            this.denoms = [];
            this.t = 0;
        end

        function this = Add_Update(this, u)
            v = this.Inverse_Apply(u);
            this.t = this.t + 1;
            this.updates(:, this.t) = u;
            this.inv_updates(:, this.t) = v;
            this.denoms(this.t) = 1 + u' * this.M * v;
        end

        function [Ax] = Apply(this, x)
            Ax = x;
            for i = 1:this.t
                u = this.updates(:, i);
                c = u' * this.M * x;
                Ax = Ax + c * u;
            end
        end

        function [Bx] = Inverse_Apply(this, x)
            Bx = x;
            for i = 1:this.t
                v = this.inv_updates(:, i);
                u = this.updates(:, i);
                c = u' * this.M * Bx / this.denoms(i);
                Bx = Bx - c * v;
            end
        end
    end
end