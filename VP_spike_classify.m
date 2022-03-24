% This function calculates the spikes needed to be added, removed, or moved
% in order to transform a source spike train into a target spike train
% according to the Victor-Purpura Distance, as well as the distance itself.

% Arguments
% F_source: vector of source spike train times
% F_target: vector of target spike train times
% tau_q: constant that determines the cost of moving spikes (large tau_q is
% a lower cost, small tau q is a larger cost)
% sigma: a function that transforms the cost of moving spikes (should be
% monotonically increasing)

% Returns
% F_add: vector of the indexes of spikes that need to be added from F_target
% F_rm: vector of the indexes of spikes that need to be removed from F_source
% F_mv: Matrix with two rows, each column corresponds to a movement of
% spikes from source to target, where the first number in the column is the
% time in the source, and the second is the time in the target
% d: The actual Victor-Purpura distance (non-negative double)

function [F_add, F_rm, F_mv, d] = VP_spike_classify(F_source, F_target, tau_q, sigma)
    n_source = length(F_source);
    n_target = length(F_target);

    D = zeros(n_source + 1, n_target + 1);
    S_add = cell(n_source + 1, n_target + 1);
    S_rm = cell(n_source + 1, n_target + 1);

    S_mv = cell(n_source + 1, n_target + 1);
    S_mv(:,:) = {zeros(0,2)};

    D(1,:) = 0:n_target;
    for j=1:n_target + 1
       S_add{1,j} = 1:j-1;
    end

    for i=2:n_source + 1
        D(i,1) = i - 1;
        S_rm{i,1} = [S_rm{i-1,1}, i-1];
        for j= 2:n_target + 1
            a1 = D(i-1,j) + 1;
            a2 = D(i,j-1) + 1;
            a3 = D(i-1,j-1) + sigma(abs(F_source(i-1)-F_target(j-1))/tau_q);
            if a1 <= a2 && a1 <= a3
                D(i,j) = a1;
                S_add{i,j} = S_add{i-1,j};
                S_rm{i,j} = [S_rm{i-1,j}, i-1];
                S_mv{i,j} = S_mv{i-1,j};
            elseif a2 <= a3
                D(i,j) = a2;
                S_add{i,j} = [S_add{i,j-1}, j-1];
                S_rm{i,j} = S_rm{i,j-1};
                S_mv{i,j} = S_mv{i,j-1};
            else
                D(i,j) = a3;
                S_add{i,j} = S_add{i-1,j-1};
                S_rm{i,j} = S_rm{i-1,j-1};
                S_mv{i,j} = [S_mv{i-1,j-1}, [i-1;j-1]];
            end 
        end
    end
    F_add = S_add{end,end};
    F_rm = S_rm{end,end};
    F_mv = S_mv{end,end};
    d = D(end,end);
end