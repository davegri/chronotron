% This function calculates the voltage and spike times of a post-synpatic
% neuron given it's presynaptic inputs and their associated weights.
% According to the SRM_0 model.

% Arguments
% input_times: the spike times for the presynpatic neurons
% input_neurons: The indexes of the pre-synpatic neurons that fired for each spike time
% t: time vector
% W: synpatic weights
% K: voltage kernel function
% tau_m: membranal time constant
% theta: threshold potential

% Returns
% V: voltage of postsynaptic neuron by time
% spk_times: spike times of postsynpatic neuron

function [V, spk_times] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta)
    spk_times = [];
    % calculate voltages of presynaptic neurons by time
    lambdas = IF_get_lambdas(length(W), input_times, input_neurons, t, K);
    
    % multiply and sum voltages by pre synpatic weights
    V = W * lambdas;
    
    % find spikes and calculate offsets to reset voltage
    offset = zeros(1,length(t));
    t_spk_idx = find(V>theta,1);
    while(~isempty(t_spk_idx))
        spk_times(end + 1) = t(t_spk_idx); % save spike time
        % calculate and subtract the offset for all future times from current spike
        offset(t_spk_idx:end) = exp((t(t_spk_idx)-t(t_spk_idx:end))/tau_m);
        V(t_spk_idx:end) = V(t_spk_idx:end) - theta*offset(t_spk_idx:end);
        % find index of next spike
        t_spk_idx = find(V(t_spk_idx+1:end)>theta,1) + t_spk_idx;
    end
end

