%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
rng(1324)

n_u = 200;
d = 3;
plot_sample = false;

h = 1/(n_u-1);
kappa_u = .2;

C_d = 4 * d;
beta_u = kappa_u^2/C_d;

[M, S] = Assemble_Mass_and_Stiffness(n_u);
M = full(M);
S = full(S);
[V,Gamma] = eig(S,M);

if plot_sample
    Lambda = (C_d/4) * beta_u * Gamma + eye(n_u);
    omega = randn(n_u,1);
    u = V * sqrt(inv(Lambda)) * omega;
    x = linspace(0,1,n_u)';
    figure,
    plot(x,u)
end

v = abs(diag(Gamma));
if d == 1
    S = v;
elseif d == 2
    [a,b] = ndgrid(v,v);
    S = a(:) + b(:);
elseif d == 3
    [a,b,c] = ndgrid(v,v,v);
    S = a(:) + b(:) + c(:);
end
Gamma_vec = sort(S,'ascend');

C = 1/6;
coeff_var = (h^2 * C * Gamma_vec) ./ (1 + beta_u * Gamma_vec).^2;
Lambda_vec = 1 + beta_u * Gamma_vec;

v = 1./Lambda_vec;
m = length(v);

total_var = sum(v.^2);
percent_var = cumsum(v.^2)/sum(v.^2);

var_error = (h^2 * C / beta_u) * sum( (Lambda_vec-1)./(Lambda_vec.^2) );
total_var_rel_error = var_error/total_var;

disp('Total variance relative error:')
disp(total_var_rel_error)


[~,T] = min(abs(percent_var-.9));

leftColor  = [0 0.4470 0.7410];      % color for v
compColor  = [0.8500 0.3250 0.0980]; % color for sqrt(coeff_var)
rightColor = [0,0,0];                % color for percent_var

indexSets = {1:T, (T+1):m};
firstPlotLastTick = [];
sharedExponent = [];

for k = 1:numel(indexSets)
    I = indexSets{k};

figure('Units','pixels','Position',[100 100 900 650])
hold on

ax = gca;

yyaxis left
ax.YColor = leftColor;
plot(I, v(I), '-', 'Color', leftColor, 'LineWidth', 5)
plot(I, sqrt(coeff_var(I)), '--', 'Color', compColor, 'LineWidth', 5)
xlim([I(1)-length(I)/100, I(end)])
ylabel('Coefficient')

yyaxis right
ax.YColor = rightColor;
plot(I, percent_var(I), '-', 'Color', rightColor, 'LineWidth', 5)
xlim([I(1)-length(I)/100, I(end)])
ylabel('Percent of Coefficient Variance')
ylim([0 1])

xlabel('Coefficient Index')
legend({'$u$ coefficient','$u-\overline{u}$ coefficient bound','Percent of Variance'}, ...
    'Location', 'east', 'Interpreter', 'latex')

set(gca, 'FontSize', 24)

ax = gca;
ax.XAxis.Exponent = 0;  % suppress automatic x10^n label

if k == 1
    xt = xticks;
    firstPlotLastTick = xt(end);

    maxTick = max(abs(xt));
    if maxTick > 0
        sharedExponent = floor(log10(maxTick));
    else
        sharedExponent = 0;
    end
else
    xt = xticks;
    xt = xt(xt > firstPlotLastTick);
    xt = [firstPlotLastTick, xt];
    xticks(xt)
end

% Use the same exponent for both plots
xt = xticks;
scaledTicks = xt / 10^sharedExponent;
xticklabels(arrayfun(@(x) sprintf('$%.3g\\times10^{%d}$', x, sharedExponent), ...
    scaledTicks, 'UniformOutput', false))
ax.TickLabelInterpreter = 'latex';

% Force consistent diagonal labels
ax.XTickLabelRotation = 35;

% Reserve fixed space for long rotated tick labels and x-label
ax.Units = 'normalized';
ax.Position = [0.13, 0.30, 0.72, 0.62];

% Move x-axis label down slightly so it is never cut off
xl = xlabel('Coefficient Index');
xl.Units = 'normalized';
xl.Position(2) = -0.20;

print(gcf, '-depsc', sprintf('coefficient_error_d_%d_k_%d.eps', d, k))

end
