clear;
close all;
clc;

load Truth_Results.mat;
load Std_OED_Results.mat;
load Seq_OED_Results.mat;

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
