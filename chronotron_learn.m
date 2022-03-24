% this function calculates the chronotron E-learning rule for a given example

% Arguments
% input_times: the spike times for the presynpatic neurons
% input_neurons: The indexes of the pre-synpatic neurons that fired for each spike time
% y0: teacher spike times
% t: time vector
% W: presynpatic weights
% K: voltage kernel function
% tau_m: membranal time constant
% theta: threshold potential
% tau_q: victor-purpura distance constant
% eta: learning rate
% gamma_r: chronotron E-learning parameter

% Returns
% deltaW: the synpatic weights update

function [deltaW] = chronotron_learn(input_times, input_neurons, y0, t, W, K, tau_m, theta, tau_q, eta, gamma_r)
    sigma = @(x) (x.^2)/2;
    [~, spk_times] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta);
    [F_add, F_rm, F_mv] = VP_spike_classify(spk_times, y0, tau_q, sigma);
    lambdas = IF_get_lambdas(length(W), input_times, input_neurons, t, K);

    add_spk_times = y0(F_add);
    rm_spk_times = spk_times(F_rm);
    if ~isempty(F_mv)
        mv_origin_spk_times = spk_times(F_mv(1,:));
        mv_target_spk_times = y0(F_mv(2,:));
        lambdas_mv_origin_spk_idxs = ismember(t, mv_origin_spk_times);
    end
    % get the indexes in the lambdas that match the spike times
    lambdas_add_spk_idxs = ismember(t, add_spk_times);
    lambdas_rm_spk_idxs = ismember(t, rm_spk_times);
    
    deltaW = zeros(1,length(W));
    % loop over existing weights
    for i=1:length(W)
        mv_spk_volt = 0;
        % calculate each part of chronotron learning rule
        add_spk_volt = sum(lambdas(i,lambdas_add_spk_idxs));
        rm_spk_volt = sum(lambdas(i,lambdas_rm_spk_idxs));
        if ~isempty(F_mv)
            mv_origin_volt = lambdas(i,lambdas_mv_origin_spk_idxs);
            mv_spk_time_diff = mv_origin_spk_times - mv_target_spk_times;
            mv_spk_volt = sum(mv_spk_time_diff.*mv_origin_volt);
        end
        % chronotron learning rule
        deltaW(i) = eta*(add_spk_volt - rm_spk_volt + (gamma_r/tau_q^2)*mv_spk_volt);
    end
end







