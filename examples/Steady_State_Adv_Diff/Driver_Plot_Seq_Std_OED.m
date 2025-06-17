clear;
close all;
clc;

load Truth_Results.mat;
load Std_OED_Results.mat;
load Seq_OED_Results.mat;

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
N = length(Jhat_std_oed);

figure;
hold on;
xlim([0 N]);
yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
plot(0:N, [Jhat_lofi; Jhat_std_oed], ".-", "Color", "#1F618D", "DisplayName", "Standard OED");
plot(0:N, [Jhat_lofi; Jhat_seq_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
xlabel("Evaluations ($N$)", "Interpreter", "latex");
ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
legend("location", "east", "Interpreter", "latex");
title("Optimization Objective over Evals");

m = length(z_lofi);
x = linspace(0, 1, m)';
% set(h, {"color"}, num2cell(copper(50*(1:N),:), 2))
figure;
plot(x, std_oed_Z{end});
title('Std OED');
figure;
plot(x, seq_oed_Z{end});
title('Seq OED');
% figure;
% plot(x, z_hifi);
% title('High-fidelity solution');
% figure;
% plot(x, reshape(cell2mat(seq_oed_mean_z), m, N));
% title('Seq OED Updates');

% for k = 1:N
%     figure(1);
%     plot(x, std_oed_Z{k});
%     title('Std OED');
%     figure(2);
%     plot(x, seq_oed_Z{k});
%     title('Seq OED');
%     pause();
% end
