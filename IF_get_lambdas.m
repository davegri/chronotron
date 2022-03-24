% this function integrates the total voltage over time for each neuron
% according to the integrate and fire model. 

% Arguments:
% N: the number of neurons
% input_times: the spike times for all neurons
% input_neurons: The indexes of the neurons that fired for each spike time
% t: time vector
% K: voltage kernel function

% Returns:
% lambdas: an N*length(t) sized matrix, where the (i,j)th place corresponds
% to the voltage of the i'th neuron at time t(j).
function [lambdas] = IF_get_lambdas(N, input_times, input_neurons, t, K)
    lambdas = zeros(N, length(t)); % init all voltages at 0
    for t_idx = 1:length(t)
        % loop over all spike times
        for input_time_idx=1:length(input_times)
            input_time = input_times(input_time_idx); % spike time
            input_neuron = input_neurons(input_time_idx); % neuron for spike time
            % sum the voltage contribution from this spike at this time
            lambdas(input_neuron, t_idx) = lambdas(input_neuron,t_idx) + K(t(t_idx) - input_time);
        end
    end
end