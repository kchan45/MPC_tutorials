%% An implementation of offset-free MPC for temperature tracking using a simple linear plasma model
% created with influences from code by Joel A. Paulson
% (c) Kimberly J. Chan, Mesbah Lab 2022


close all; clear; clc;
%% add package paths
addpath(genpath('./config'), genpath('./utils'))

%% user inputs
economic_flag = false;  % economic objective flag (CEM setpoint tracking for plasma case studies)
robust_flag = false;    % robust MPC implementation flag (multistage MPC)
offset_flag = true;     % offset-free MPC (reference tracking)
offset_style = 'limon'; % offset-free MPC style ('limon' or 'rawlings')
save_file = [];     % filename to save data for future use; leave empty if no save is desired
font_size = 15;     % font size for plotting
line_width = 3;     % line width for plotting

%% load problem data
prob_info = offsetfree_load_problem_info();

%% set up observer
ekf = EKF();
ekf = ekf.get_observer(prob_info);

%% set up controller
c = OffsetFreeMPC(prob_info);
c = c.get_mpc();
[c, res, feas] = c.solve_mpc(); % do test run

%% initialize simulation
Nsim = ceil(5*60/0.5);
sim = Simulation(Nsim, prob_info);
% run simulation
[sim, sim_data] = sim.run_closed_loop(c, [], ekf, economic_flag, robust_flag, offset_flag);
disp('sim_data fields:')
disp(fieldnames(sim_data))

%% plot data
ctime = sim_data.ctime;
fprintf('Total Runtime: %.2f s\n', sum(ctime))
fprintf('Average Runtime: %.4f s\n', mean(ctime))

sim_ref = sim_data.Yrefsim + prob_info.xss(1); % reference was saved as deviation variable using the PREDICTION model as basis
Tplot = sim_data.Xsim(1,1:end-1) + prob_info.xssp(1); % states are taken as deviation variables using the PLANT model as basis
t = linspace(0,length(sim_ref),length(sim_ref))*prob_info.ts;

figure()
subplot(111)
hold on
plot(t, sim_ref, 'k--', 'LineWidth', line_width, 'DisplayName', 'Reference')
T_max = (prob_info.x_max(1)+prob_info.xss(1))*ones(size(Tplot));
plot(t, T_max, 'r--', 'LineWidth', line_width, 'DisplayName', 'Constraint')
plot(t, Tplot, 'b-', 'LineWidth', line_width, 'DisplayName', 'Surface Temperature')
xlabel('Time (s)')
ylabel('Surface Temperature ($^\circ$C)', 'Interpreter', 'latex')
legend()
set(gca, 'Fontsize', font_size)

figure()
subplot(111)
hold on
labels = {'Power', 'Flow Rate'};
for i = 1:length(labels)
    u_plot = sim_data.Usim(i,:)+prob_info.ussp(i);
    stairs(t, u_plot, '-', 'LineWidth', line_width, 'DisplayName', labels{i})
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