classdef SSA_System < handle

    properties
        k1;
        k2;
        k3;
        k4;
        k5;
        k6;
        T;
        N;
        vol;
        nA;
        con;
    end

    methods
        function this = SSA_System(con)
            this.k1 = con.k1;
            this.k2 = con.k2;
            this.k3 = con.k3;
            this.k4 = con.k4;
            this.k5 = con.k5;
            this.k6 = con.k6;
            this.T = con.T;
            this.N = con.n_t;
            this.vol = con.vol;
            this.nA = con.nA;
            this.con = con;
        end

        function [Y] = SSA_Mean(this,z,num_samples)
            Y_tmp = this.SSA_Solve(z);
            Y_samps = zeros(length(Y_tmp),num_samples);
            Y_samps(:,1) = Y_tmp;
            for k = 2:num_samples
                Y_samps(:,k) = this.SSA_Solve(z);
            end
            Y = mean(Y_samps,2);
        end

        function [Y] = SSA_Solve(this,z)

            X = this.con.h(z)*this.nA*this.vol/this.con.state_scale;

            V = zeros(9,6);
            V(:,1) = [-1  -1  0   0   0   0   1   0   0];
            V(:,2) = [0   0   -1  -1  0   0   0   1   0];
            V(:,3) = [-1  0   -1  0   1   0   0   0   0];
            V(:,4) = [0  -1   0   -1  0   1   0   0   0];
            V(:,5) = [0   0   0   0   -1  -1  0   0   1];
            V(:,6) = [0   0   0   0   0   0   -1  -1  1];

            %%%%%%%%%% Parameters and Initial Conditions %%%%%%%%%
            c = zeros(6,1);
            c(1) = this.k1/(this.nA*this.vol); 
            c(2) = this.k2/(this.nA*this.vol); 
            c(3) = this.k3/(this.nA*this.vol); 
            c(4) = this.k4/(this.nA*this.vol); 
            c(5) = this.k5/(this.nA*this.vol); 
            c(6) = this.k6/(this.nA*this.vol); 

            t = 0;
            count = 1;
            tvals(1) = 0;
            Xvals(:,1) = X;
            while t < this.T
                a(1) = c(1)*X(1)*X(2);
                a(2) = c(2)*X(3)*X(4);
                a(3) = c(3)*X(1)*X(3);
                a(4) = c(4)*X(2)*X(4);
                a(5) = c(5)*X(5)*X(6);
                a(6) = c(6)*X(7)*X(8);
                asum = sum(a);
                j = min(find(rand<cumsum(a/asum)));
                tau = log(1/rand)/asum;
                X = X + V(:,j);

                count = count + 1;
                t = t + tau;
                tvals(count) = t;
                Xvals(:,count) = X;
            end

            t = linspace(0,this.T,this.N)';
            Y = interp1(tvals,Xvals',t)';
            Y = Y(:)*this.con.state_scale/(this.nA*this.vol);
        end

    end
end