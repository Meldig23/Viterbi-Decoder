%% Define parameters for a two-stage trellis diagram with 4-PAM
syms z;
transfer_function = [0.6, -1, 0.8] / sqrt(2); % F(z) coefficients for channel (scaled by sqrt(2))
No_symbols = 1e5; % Total number of symbols (including termination symbols)
mu = 2; % Channel memory

% Parameters for SNR and 4-PAM symbol settings
N = 1 * mu;
symbols = [-3 -1 1 3]; % 4-PAM symbol set
initialstate = [-1 -1]; % Symbols for initial state and final termination
M = length(symbols);
Eb = mean(symbols.^2); % Average energy per symbol for given constellation
SNR_db = 0:2:16; % SNR values in dB for analysis

% Generate random sequence of 4-PAM symbols
tx_symbols = symbols(randi([1, M], 1, No_symbols));
tx_symbols = [initialstate tx_symbols initialstate]; % Add termination symbols

%% Define states: two previous symbols form the state
state1 = repmat(symbols, 1, M); % First state (previous symbol)
state2 = repelem(symbols, M); % Second state (current symbol)
statesdiag = [state2; state1]; % State combinations (2 symbols)
states = flip(statesdiag, 2);
numStates = M^2;
theta = linspace(0, 2 * pi, numStates + 1);
statePositions = [cos(theta(1:end-1)); sin(theta(1:end-1))]';

%% Plot State Diagram
figure;
hold on;
title('State Diagram for 4-PAM, Two-Stage Trellis');
axis equal;

% Plot nodes (states)
for i = 1:numStates
    plot(statePositions(i, 1), statePositions(i, 2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'cyan');
    text(statePositions(i, 1), statePositions(i, 2), num2str([statesdiag(1, i), statesdiag(2, i)]), 'HorizontalAlignment', 'center', 'FontSize', 8);
end

% Plot transitions
colors = lines(M); % Different colors for each input symbol
for i = 1:numStates
    for j = 1:M
        currentInput = symbols(j);
        nextState = [currentInput; statesdiag(1, i)];
        nextIndex = find(ismember(statesdiag', nextState', 'rows'));
        if ~isempty(nextIndex)
            dx = statePositions(nextIndex, 1) - statePositions(i, 1);
            dy = statePositions(nextIndex, 2) - statePositions(i, 2);
            quiver(statePositions(i, 1), statePositions(i, 2), dx, dy, 0, 'Color', colors(j, :), 'MaxHeadSize', 0.3);
            f_value = transfer_function(1) * symbols(j) + transfer_function(2) * statesdiag(1, i) + transfer_function(3) * states(2, i);
            label = sprintf('%d/%.2f', currentInput, f_value);
            text(statePositions(i, 1) + 0.5*dx, statePositions(i, 2) + 0.5*dy, label, 'FontSize', 8, 'Color', colors(j, :));
        end
    end
end
hold off;

%% Pass symbols through ISI channel
channel_output = zeros(1, length(tx_symbols) - 2);
for i = 1:length(tx_symbols) - 2
    channel_output(i) = sum(tx_symbols(i + 2:-1:i) .* transfer_function);
end

% Alternative ISI channel output using convolution
c1 = conv(tx_symbols, transfer_function, "valid");
No_symbols = No_symbols +2;

%% Define state-space and mapping for Viterbi decoding
states = combinations(symbols, symbols).Variables;
arr = zeros(symbols(end) - symbols(1) + 1, symbols(end) - symbols(1) + 1);
counter = 1;
for stage = symbols
    for j = symbols
        arr(stage + 1 - symbols(1), j + 1 - symbols(1)) = counter;
        counter = counter + 1;
    end
end

% Compute channel output for each state transition
outputs = zeros(length(states), length(states));
for stage = 1:length(outputs)
    for j = 1:length(outputs)
        l1 = symbols(ceil(j / length(symbols)));
        if l1 ~= 0
            outputs(stage, j) = sum([l1, states(stage, :)] .* transfer_function);
        end
    end
end

%% Viterbi Algorithm
figure;
matrix = inf(length(states), length(states));
weight_of_survivor = zeros(1, length(states));
survivor_paths = zeros(length(states), No_symbols + 1);
SNR = 10.^(SNR_db / 10);
N0 = Eb ./ SNR;
SEP = zeros(1, length(N0));
N_array = mu * [1, 2, 4, 5, 10];

for N = N_array
    for a = 1:length(N0)
        noise = sqrt(N0(a) / 2) * randn(1, No_symbols);
        rx_signal = channel_output + noise;
        weight_of_survivor = inf(1, length(states));
        survivor_paths = zeros(length(states), No_symbols + 1);
        matrix = inf(length(states), length(states));
        decoded_symbols = zeros(1, No_symbols);
        c = 1;

        % Main loop for Viterbi algorithm
        for stage = 1:No_symbols + 1
            possible_states = 1:length(states);
            weight_of_survivor = inf(1, length(states));
            for j = possible_states
                [val, previous_index] = min(matrix(:, j));
                weight_of_survivor(j) = val;
                survivor_paths(j, stage) = previous_index;
            end

            [~, index] = min(weight_of_survivor);
            if stage == 1
                weight_of_survivor(1) = 0;
            end

            if stage ~= No_symbols + 1
                matrix = inf(length(states), length(states));
                for j = possible_states
                    for input = symbols
                        k1 = next_state_index(j, input, states, arr, symbols);
                        k2 = outputs(j, k1);
                        matrix(j, k1) = (rx_signal(stage) - k2)^2;
                    end
                end
                matrix = matrix + weight_of_survivor';
            end
        end

        SEP(a) = sum(tx_symbols(3:end - 2) ~= decoded_symbols(1:end - 2)) / (No_symbols);
    end
    hold on;
    semilogy(SNR_db, SEP, 'LineWidth', 2);
end

legend(('N = ') + string(N_array / mu) + ('mu'));
xlabel("SNR(dB)", "FontWeight", "bold");
ylabel("Symbol Error Probability (SEP)", "FontWeight", "bold");
title('SEP vs SNR for different values of N');

%% Helper Functions
function cal = next_state_index(current_state_index, input, states, arr, symbols)
    cal = arr(input + 1 - symbols(1), states(current_state_index, 1) + 1 - symbols(1));
end
