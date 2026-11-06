%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runs every continuation-improvements experiment in turn.  Each experiment
% is also a standalone script that you can run on its own:
%
%   Driver_CI_Mean_To_Sample      mean-to-sample vs nominal-to-sample
%   Driver_CI_Adaptive_Mean       time-step heuristic, posterior mean (fast)
%   Driver_CI_Adaptive_Samples    time-step heuristic, mean-to-sample
%   Driver_CI_Mean_Newton_Polish  final Newton steps, posterior mean (fast)
%
% They share CI_Setup (settings, posterior sampling, mean solves) and
% CI_Sample_Reference (cached sample reference).  Each saves a ci-*-results.mat
% file and has a matching Plot_CI_* script, which can be rerun on its own.

Driver_CI_Mean_To_Sample;
Driver_CI_Adaptive_Mean;
Driver_CI_Adaptive_Samples;
Driver_CI_Mean_Newton_Polish;
