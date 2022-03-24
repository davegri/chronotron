tic
clear;
clc;

% load data
load('train_data.mat');

% parameters
R           = 2e4;              %[Ohm] resistance
C           = 1e-6;             %[F] capacitance
tau_m       = R*C;              %[sec] membrane time constant
tau_s       = tau_m/4;          %[sec] synaptic time constant
tau_q = tau_m;                  % spike movement cost
theta = 30e-3;                  % [V] threshold potential
eta = 1e-4;                     % learning rate
gamma_r = tau_q; 
epochs = 10;
t_init = 0;
t_fin = 0.5;
delta_t = 1e-4;
t = t_init:delta_t:t_fin;

% We round the time vector to the nearest tick size in order to make sure
% Comparisons between times work correctly
t = arrayfun(@(x) roundn(x, log10(delta_t)), t);

% define params for voltage kernel function
alpha = tau_m/tau_s;
kappa = alpha^(alpha/(alpha-1))/(alpha-1); % normalize K to max of 1V
K = @(T) (T > 0).*kappa.*(exp(-T/tau_m) - exp(-T/tau_s)); % Voltage Kernel function
W = theta/2*rand(1,N); % pre-synpatic weights
sigma = @(x) (x.^2)/2; % spike distance cost function

% round time data to nearest delta_t
for i=1:length(Samples)
    Samples(i).times = arrayfun(@(x) roundn(x, log10(delta_t)), Samples(i).times);
    Samples(i).y0 = arrayfun(@(x) roundn(x, log10(delta_t)), Samples(i).y0);
end

% Train chronotron on samples
for i=1:epochs
    random_indx = randperm(length(Samples)); % random permutation of samples
    for j = 1:length(random_indx)
        sample = Samples(random_indx(j));
        % calculate weight change for sample
        deltaW = chronotron_learn(sample.times, sample.neurons,...
            sample.y0, t, W, K, tau_m, theta, tau_q, eta, gamma_r);
        % update weights
        W = W + deltaW;
    end
end

% plot chrontron learning results for 3 random samples
sgtitle("Voltage of chronotron relative to time - after learning" + newline + "3 random samples")
random_idx = randperm(length(Samples)); % random permutation of samples
for i=1:3
    % calculate voltage and spikes for random sample
    sample = Samples(random_idx(i));
    [V, spk_times] = IF_sim(sample.times, sample.neurons, t, W, K, tau_m, theta);
    
    % plot random sample voltage as a function of time
    subplot(3,1,i);
    plot(t*1e3,V*1e3);
    yline(theta*1e3,"--");
    xlabel("time (ms)");
    ylabel("voltage (mV)");
    title("sample number " + random_idx(i));
    

    % calculate VP distance for sample and display as subtitle
    [~, ~, ~, d] = VP_spike_classify(spk_times, sample.y0, tau_q, sigma);
    subtitle("Victor-Purpura Distance: " + d);
    
    % plot markers for teacher spikes and chronotron spikes
    hasTeacherSpks = ~isempty(sample.y0);
    hasSpks = ~isempty(spk_times);
    hold on;
    % we check this condition to make sure the legend entries match the
    % spikes, in the situation where some type of spike dosen't exist
    if hasTeacherSpks && ~hasSpks
        plot(sample.y0*1e3, (theta*1e3-5)*ones(1,length(sample.y0)),'bl*');
        legend("Voltage (mV)", "Threshold (mV)", "teacher spikes", "actual spikes");
    else
        plot(spk_times*1e3, theta*1e3*ones(1,length(spk_times)),'r*');
        plot(sample.y0*1e3, (theta*1e3-5)*ones(1,length(sample.y0)),'bl*');
        legend("Voltage (mV)", "Threshold (mV)", "actual spikes", "teacher spikes");
    end
    hold off;
end

% calculate VP-Distances for all samples
vp_distances = zeros(1,length(Samples));
for i=1:length(Samples)
   [~, spk_times] = IF_sim(Samples(i).times, Samples(i).neurons, t, W, K, tau_m, theta);
   [~, ~, ~, d] = VP_spike_classify(spk_times, Samples(i).y0, tau_q, sigma);
   vp_distances(i) = d;
end

% mean and std for VP distances
mean_vp_dist = mean(vp_distances);
std_vp_dist = std(vp_distances);

disp("Mean Victor-Purpura distance for samples is: " + mean_vp_dist...
    + newline + "With STD: " + std_vp_dist);
toc