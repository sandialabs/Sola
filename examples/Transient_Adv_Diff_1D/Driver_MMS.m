clear
close all
clc
addpath(genpath('../../src'))

plot_solution = false;

m = 200;
N = 51;
T = 1;
n = N*m;
adv_diff = Adv_Diff(m,n,T,N);
x = adv_diff.x;
t = adv_diff.t_mesh;
z0 = zeros(m,N);
for k = 1:N
    tk = t(k);
    z0(:,k) = cos(2*pi*x) + 4*pi^2*tk*cos(2*pi*x) - 2*pi*tk*sin(2*pi*x);
end
z0 = z0(:);

u = adv_diff.State_Solve(z0);
u_reshape = reshape(u,m,N);

u_true = zeros(m,N);
for k = 1:N
    tk = t(k);
    u_true(:,k) = tk*cos(2*pi*x);
end

if plot_solution
    figure,
    for k = 1:N
        plot(x,u_reshape(:,k),'-',x,u_true(:,k),'--','LineWidth',3)
        pause(.05)
    end
end

disp(norm(u-u_true(:))/norm(u_true(:)))