% Clear everything
clc; clear all; close all;

%% Params

% Time horizon considered (in years) for reliability calculation
T_HORIZON = 800;
% Time mission considered in years
T_M = 200;
% Time for calculating conditional reliability in years
T_COND = 20;
% Number of Monte Carlo iterations
N_MC_SIM = 5000;
% Number of system components
N_COMP = 2;
% Failure rate of each component [years^(-1)]
LAMBDA = 0.01;
% For reproducibility
rng(2561);
% System state matrix
system_state_matrix = zeros(N_MC_SIM, T_HORIZON + 1);
% Fix equal to 1 element k, 1 of the state matrix for each row since R(0)=1
system_state_matrix(1:N_MC_SIM, 1) = 1;
% All possible existing combinations of system state
ALL_POSSIBLE_STATES = dec2bin(0:(2^N_COMP-1)) - '0';
% System type
SYS_TYPE = 'PARALLEL';

%% Set system type

% Working and not working states: working statets are those with 1
if strcmp(SYS_TYPE, 'SERIES')
    SYS_STATES = reshape([0 0 0 1], [2^N_COMP, 1]);
elseif strcmp(SYS_TYPE, 'PARALLEL')
    SYS_STATES = reshape([0 1 1 1], [2^N_COMP, 1]);
end

%% Monte Carlo iterations

% Monte Carlo iterations outer loop
for k = 1:N_MC_SIM
    
    % Components state matrix: Ncomp x (T_horizon + 1)
    components_state_matrix = zeros(N_COMP, T_HORIZON + 1);
    % Fix equal to 1 element (k, 1) of the state matrix for each row k
    % since reliability is 1 at beginning of the period
    components_state_matrix(1:N_COMP, 1) = 1;

    % Sample tf for each component
    for i = 1:N_COMP
        % Sample tf (tf is in years)
        tf = expinv(rand, 1/LAMBDA);
        % From year 1 to ceil(tf) component has worked, set 1.
        % Note that index is shifted! Year 1 is at index 2.
        components_state_matrix(i, 2:min(ceil(tf), T_HORIZON + 1)) = 1;
        % From year ceil(tf) to year T_m component has not worked, set 0
        components_state_matrix(i, min(ceil(tf) + 1, T_HORIZON + 1):(T_HORIZON + 1)) = 0;
    end

    % Determine system history.
    for i = 2:(T_HORIZON + 1)
        
        % State system at time i is equal to column i of
        % components_state_matrix
        current_system_state = components_state_matrix(1:N_COMP, i);
        
        % Find if current_system_state is a working state
        for j = 1:2^N_COMP
            % This finds the corresponding state from all possible states
            index_sum = sum(reshape(current_system_state, [1, N_COMP]) == ALL_POSSIBLE_STATES(j, 1:N_COMP));
            % Once found the corresponding state we exit the loop
            if index_sum == N_COMP
                break
            end
        end
    
        % Corresponding state is the j-th element of SYS_STATES
        system_status = SYS_STATES(j);
        
        % If system_status is 0 system is not working
        if ~ system_status
            system_state_matrix(k, i:(T_HORIZON+1)) = 0; 
            break;
        end
        
        % Else the system is still working
        system_state_matrix(k, i) = 1;
    
    end   
end

%% Calculate quantities of interest

% System reliability R(t)
R = sum(system_state_matrix)/N_MC_SIM;

% MTTF: integrating with unit increment
MTTF = trapz(R);

% Reliability at T_M
R_TM = R(T_M + 1);

% Conditional reliability at T = 20 years -> R(t + T)/R(T)
cond_R = R((T_COND*2 + 1):(T_HORIZON + 1))/R(T_COND + 1);
cond_t_axis = T_COND:(T_HORIZON - T_COND);

% Residual reliability
residual_reliability = mean([R(floor(MTTF) + 1) R(ceil(MTTF) + 1)]);

%% Print results

% Print MTTF
sprintf('MTTF %f years', MTTF)

% Reliability at mission time T_M
sprintf('Reliability at mission time %f', R_TM)

% Residual reliability
sprintf('Residual reliability %f', residual_reliability)

% Estimated probability that system survives 20 years
sprintf('Estimated reliability at 20 years %f', R(21))

%% Plot results

% Plot of reliability
plot(0:T_HORIZON, R, 'LineWidth', 1.5, 'color', 'b')
xlabel('Time [years]'); ylabel('R(t)');
xlim([0 , T_HORIZON])
ylim([0, 1])
grid
hold on

% Real reliability
if strcmp(SYS_TYPE, 'SERIES')
    R_true = exp(-2*LAMBDA*(0:T_HORIZON));
    plot(0:T_HORIZON, R_true, 'LineWidth', 1.5, 'color', 'g')
    plot(cond_t_axis, R_true((T_COND*2 + 1):(T_HORIZON + 1))/R_true(T_COND + 1), 'LineWidth', 1.5, 'color', 'cyan');
elseif strcmp(SYS_TYPE, 'PARALLEL')
    R_true = 2*exp(-LAMBDA*(0:T_HORIZON)) - exp(-2*LAMBDA*(0:T_HORIZON));
    plot(0:T_HORIZON, R_true, 'LineWidth', 1.5, 'color', 'g')
    plot(cond_t_axis, R_true((T_COND*2 + 1):(T_HORIZON + 1))/R_true(T_COND + 1), 'LineWidth', 1.5, 'color', 'cyan');
end

% Plot conditional reliability P(tf > t + T | tf > T)
plot(cond_t_axis, cond_R, 'LineWidth', 1.5, 'color', 'r')

% Plot a 1 indicator for reference %refline(0, 1)
line([T_COND T_COND], [0 1], 'col', 'k', 'LineWidth', 0.5, 'DisplayName', 'Marker of year 20')

% Legend
if strcmp(SYS_TYPE, 'SERIES')
    title('Series system Reliability')
elseif strcmp(SYS_TYPE, 'PARALLEL')
    title('Parallel system Reliability')
end

legend('MC Reliability', 'True Reliability', 'True conditional reliability | 20 years', 'MC conditional reliability | 20 years','Marker of year 20')
