clear;
clc;

%% Test VP_spike_classify
F_source = [30,110,150,200];
F_target = [50,100,115,150];
tau_q = 5;
sigma = @(x) x;
[F_add, F_rm, F_mv, d] = VP_spike_classify(F_source,F_target,tau_q,sigma);
source_move_times = F_source(F_mv(1,:));
target_move_times = F_target(F_mv(2,:));
moves = "";
for j=1:length(target_move_times)
    moves = moves + source_move_times(j) + " -> " + target_move_times(j) + ", ";
end
disp("Testing VP_spike_classify with params:" + newline...
    + "source spike times: " + mat2str(F_source) + newline...
    + "target spike times: " + mat2str(F_target) + newline...
    + "spike movement cost (tau_q): " + tau_q + newline...
    + "sigma movment function (sigma): identity" + newline + newline...
    + "Results:" + newline...
    + "Times of spikes to add from target: " + mat2str(F_target(F_add)) + newline...
    + "Times of spikes to remove from source: " + mat2str(F_source(F_rm)) + newline...
    + "Times of spike movements (source -> target): " + moves + newline...
    + "Victor-Purpura distance: " + d);
   
%% Test IF_get_lambdas
% variables for test
R           = 2e4;              %[Ohm] resistance
C           = 1e-6;             %[F] capacitance
tau_m       = R*C;              %[sec] membrane time constant
tau_s       = tau_m/4;          %[sec] synaptic time constant
alpha   = tau_m/tau_s;
kappa   = alpha^(alpha/(alpha-1))/(alpha-1); % normalize K to max of 1V
K       = @(T) (T > 0).*kappa.*(exp(-T/tau_m) - exp(-T/tau_s)); % Voltage Kernel function
N = 4; % number of neurons
input_times = [0.08,0.1,0.15,0.2,0.25,0.27,0.35,0.45]; % spike times
input_neurons = [1,1,2,2,2,3,4,3];
t = 0:1e-4:0.5;

lambdas = IF_get_lambdas(N, input_times, input_neurons, t, K);
% print test params
disp(newline + "Testing IF_get_lambdas with params: " + newline...
    + "Number of neurons: " + N + newline...
    + "Input spike times: " + mat2str(input_times) + newline...
    + "Input neurons: " + mat2str(input_neurons) + newline...
    + "Using standard tempotron voltage kernel function")
% graph results
hold on;
sgtitle("Testing Function: IF\_get\_lambdas" + newline...
    + "voltage of each neuron according to time")
figure(1);
for i=1:N
    subplot(N,1,i);
    voltage_mv = lambdas(i,:)*1e3;
    t_ms = t*1e3;
    plot(t_ms,voltage_mv);
    title("neuron number " + i);
    xlabel("time (ms)");
    ylabel("voltage (mV)");
end
hold off;

%% Test IF_sim
% variables for test
R           = 2e4;              %[Ohm] resistance
C           = 1e-6;             %[F] capacitance
tau_m       = R*C;              %[sec] membrane time constant
tau_s       = tau_m/4;          %[sec] synaptic time constant
alpha   = tau_m/tau_s;
kappa   = alpha^(alpha/(alpha-1))/(alpha-1); % normalize K to max of 1V
K       = @(T) (T > 0).*kappa.*(exp(-T/tau_m) - exp(-T/tau_s)); % Voltage Kernel function
input_times = [0.08,0.1,0.15,0.2,0.25,0.27,0.35,0.45]; % spike times
input_neurons = [1,1,2,2,2,3,4,3];
t = 0:1e-4:0.5;
W = [0.03,0.02,0.03,0.03]; % pre synpatic weights
theta = 30e-3; % threshold potential

% print test params
disp(newline + "Testing IF_sim with params: " + newline...
    + "Input spike times: " + mat2str(input_times) + newline...
    + "Input neurons: " + mat2str(input_neurons) + newline...
    + "Using standard tempotron voltage kernel function" + newline...
    + "Neuron weights: " + mat2str(W) + newline...
    + "Membrane time constant (tau_m): " + tau_m + newline...
    + "Threshold potential (theta): " + theta)

[V, spk_times] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta);
figure(2);
t_ms = t*1e3;
V_mv = V*1e3;
plot(t_ms,V_mv)
title("Testing Function: IF\_sim")
subtitle("integration of presynpatic spikes according to the SRM_0 model"...
    + newline + "spike times (ms): " + mat2str(spk_times*1e3))
xlabel("time (ms)");
ylabel("voltage (mV)");
yline(theta*1e3,"--");
legend("Voltage (mV)", "Threshold (mV)");

%% Test chronotron_learn
% variables for test
R           = 2e4;              %[Ohm] resistance
C           = 1e-6;             %[F] capacitance
tau_m       = R*C;              %[sec] membrane time constant
tau_s       = tau_m/4;          %[sec] synaptic time constant
alpha   = tau_m/tau_s;
kappa   = alpha^(alpha/(alpha-1))/(alpha-1); % normalize K to max of 1V
K       = @(T) (T > 0).*kappa.*(exp(-T/tau_m) - exp(-T/tau_s)); % Voltage Kernel function
input_times = [0.08,0.1,0.15,0.2,0.25,0.27,0.35,0.45]; % spike times
input_neurons = [1,1,2,2,2,3,4,3];
t = 0:1e-4:0.5;
theta = 30e-3;
W = [0.03,0.02,0.03,0.03]; % pre-synpatic weights
y0 = [10*1e-3,270*1e-3,360*1e-3,370*1e-3];
tau_q = tau_m/10000000; % spike movement cost
eta = 1e-3; % learning rate
gamma_r = tau_q; 

deltaW = chronotron_learn(input_times, input_neurons, y0, t, W, K, tau_m, theta, tau_q, eta, gamma_r);
disp(newline + "Testing chronotron_learn: " + newline...
    + "Resulting deltaW: " + mat2str(deltaW));

