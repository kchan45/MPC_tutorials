%% An implementation of economic MPC for the plasma thermal effect/dose delivery problem
% created with influences from code by Angelo D. Bonzanini and Joel A. Paulson
% (c) Kimberly Chan, Mesbah Lab 2022


close all; clear; clc;
%% add package paths
addpath(genpath('./config'), genpath('./utils'))

%% user inputs
economic_flag = true;   % economic objective flag (CEM setpoint tracking for plasma case studies)
robust_flag = false;    % robust MPC implementation flag (multistage MPC)
offset_flag = false;    % offset-free MPC (reference tracking)
save_file = [];     % filename to save data for future use; leave empty if no save is desired
font_size = 15;     % font size for plotting
line_width = 3;     % line width for plotting

%% load problem data
prob_info = economic_load_problem_info();

%% set up controller
c = EconomicMPC(prob_info);
c = c.get_mpc();
[c, res, feas] = c.solve_mpc(); % do test run

%% initialize simulation
Nsim = 200;
sim = Simulation(Nsim, prob_info);
% run simulation
[sim, sim_data] = sim.run_closed_loop(c, [], [], economic_flag, robust_flag, offset_flag);
disp('sim_data fields:')
disp(fieldnames(sim_data))

%% plot data

sim_ref = sim_data.Yrefsim;
CEMsim = sim_data.CEMsim(:,1:end-1);
mask = CEMsim < sim_ref;
ctime = sim_data.ctime(mask);
fprintf('Total Runtime: %.2f s\n', sum(ctime))
fprintf('Average Runtime: %.4f s\n', mean(ctime))

CEMplot = CEMsim(mask);
ref_plot = sim_ref(mask);
Tplot = sim_data.Xhat(1,1:end-1) + prob_info.xssp(1);
Tplot = Tplot(mask);
t = linspace(0,length(CEMplot),length(CEMplot))*prob_info.ts;

figure(1)
subplot(111)
hold on
plot(t, ref_plot, 'k--', 'LineWidth', line_width, 'DisplayName', 'Reference')
plot(t, CEMplot, 'b-', 'LineWidth', line_width, 'DisplayName', 'CEM')
xlabel('Time (s)')
ylabel('CEM (min)')
legend()
set(gca, 'Fontsize', font_size)

figure(2)
subplot(111)
hold on
T_max = (prob_info.x_max(1)+prob_info.xss(1))*ones(size(Tplot));
plot(t, T_max, 'r--', 'LineWidth', line_width, 'DisplayName', 'Constraint')
plot(t, Tplot, 'b-', 'LineWidth', line_width, 'DisplayName', 'Surface Temperature')
xlabel('Time (s)')
ylabel('Surface Temperature ($^\circ$C)', 'Interpreter', 'latex')
legend()
set(gca, 'Fontsize', font_size)

figure(3)
subplot(111)
hold on
labels = {'Power', 'Flow Rate'};
for i = 1:length(labels)
    u_plot = sim_data.Usim(i,:)+prob_info.ussp(i);
    stairs(t, u_plot(mask), '-', 'LineWidth', line_width, 'DisplayName', labels{i})
end
xlabel('Time (s)')
ylabel('Inputs')
legend()
set(gca, 'Fontsize', font_size)

%% save data
if ~isempty(save_file)
    save(save_file, 'sim_data', 'prob_info')
end

%% remove package paths
rmpath(genpath('./config'), genpath('./utils'))