%% Assignment: Expected Matching Accuracy with Partial Elimination
% Description:
% You are asked to simulate and analyze a word-definition matching task.
% For a given number of words (n), you will compare the expected number of correct matches:
%   (1) when guessing randomly,
%   (2) when k wrong options are eliminated *per word*,
%   (3) when y wrong word-definition *pairs* are eliminated globally.
% The goal is to empirically and/or analytically calculate how accuracy improves
% with increasing elimination power, and how this scales with n.
% You will write functions and run experiments that compare these strategies for
% different values of n (e.g., 3 ≤ n ≤ 8) and different levels of elimination.
% This project combines simulation (via Monte Carlo), combinatorics (via permutation filtering),
% and statistics (via averaging over trials).
% Recommended: Use helper functions to structure your code cleanly.
%% Setup
% You may use built-in functions such as perms, randperm, and logical indexing.
% Avoid toolboxes unless specifically instructed.
% TODO (1): Initialize parameters
maxN = 6;           % Try values up to maxN
trials = 5000;     % Number of random permutations to average over
k_list = 1:3;       % Number of eliminated options per word (0, 1, 2)
max_y_needed = maxN * max(k_list);
y_list = 1:max_y_needed; % Number of globally eliminated word-definition pairs
% Store results in matrices:
E_random = zeros(1, maxN);
E_perword = zeros(length(k_list), maxN);
E_global = zeros(length(y_list), maxN);
% For the scaled graph (new additions)
E_scaled_k_elim = zeros(length(k_list), maxN); % Expected correct matches for scaled elimination
y_from_k = zeros(length(k_list), maxN); % Store the calculated y for each k and n
%% Part 1: Random Guessing Expectation
% TODO (2): For each n = 1 to maxN, simulate random assignments and compute expected # correct
for n = 1:maxN
    matchCount = 0;
    for t = 1:trials
        guess = randperm(n);
        matchCount = matchCount + sum(guess == 1:n);
    end
    E_random(n) = matchCount / trials;
end
%% Part 2.1: Eliminate k Wrong Options per Word
% TODO (3): For each value of n and k, simulate with each word having k wrong choices removed
for n = 1:maxN
    for ki = 1:length(k_list)
        k = k_list(ki);
        % Skip invalid elimination levels
        if k >= n
            E_perword(ki, n) = NaN;
            continue
        end
        matchCount = 0;
        validTrials = 0;
        for t = 1:trials
            allowed = true(n); % allowed(i,j) = true if word i can be matched with def j
            for i = 1:n
                choices = setdiff(1:n, i); % exclude the correct match for word i
                % Ensure k does not exceed available wrong choices
                if k > length(choices)
                    k_actual = length(choices);
                else
                    k_actual = k;
                end
                
                % Randomly select k_actual wrong options to forbid
                forbidden_indices = randperm(length(choices), k_actual);
                forbidden_defs = choices(forbidden_indices);
                allowed(i, forbidden_defs) = false;
            end
            
            perms = allValidPerms(allowed);
            if isempty(perms)
                continue
            end
            chosen = perms(randi(size(perms,1)), :);
            matchCount = matchCount + sum(chosen == 1:n);
            validTrials = validTrials + 1;
        end
        if validTrials > 0
            E_perword(ki, n) = matchCount / validTrials;
        else
            E_perword(ki, n) = NaN;
        end
    end
end
%% Part 2.2: Eliminate k Options per Word (Wrongly) (No specific fix needed, seems intentional)
% This part simulates eliminating k options per word without guaranteeing the correct match is kept.
% The logic seems to align with this "wrong" interpretation.
E_perword_wrong = zeros(length(k_list), maxN);
for n = 1:maxN
    for ki = 1:length(k_list)
        k = k_list(ki);
        if k >= n
            E_perword_wrong(ki, n) = NaN;
            continue
        end
        matchCount = 0;
        validTrials = 0;
        for t = 1:trials
            allowed = true(n);
            for i = 1:n
                if k >= n
                    forbidden = randperm(n, n-1); % At least one option left
                else
                    forbidden = randperm(n, k);
                end
                allowed(i, forbidden) = false;
            end
            perms = allValidPerms(allowed);
            if isempty(perms)
                continue
            end
            chosen = perms(randi(size(perms,1)), :);
            matchCount = matchCount + sum(chosen == 1:n);
            validTrials = validTrials + 1;
        end
        if validTrials > 0
            E_perword_wrong(ki, n) = matchCount / validTrials;
        else
            E_perword_wrong(ki, n) = NaN;
        end
    end
end
%% Part 3.1: Eliminate y Incorrect Global Pairs
% TODO (4): For each value of n and y, simulate by forbidding y random pairs from all n^2
for n = 1:maxN
    for yi = 1:length(y_list)
        y = y_list(yi);
        maxForbiddenPairs = n^2 - n;
        if y > maxForbiddenPairs
            E_global(yi, n) = NaN;
            continue
        end
        matchCount = 0;
        validTrials = 0;
        for t = 1:trials
            allowed = true(n);
            
            % Generate linear indices for all possible pairs
            allIdx = 1:n^2;
            % Generate linear indices for correct pairs (where word i matches def i)
            correctPairsLinearIdx = sub2ind([n n], 1:n, 1:n);
            % Get linear indices of all *incorrect* pairs
            incorrectPairsLinearIdx = setdiff(allIdx, correctPairsLinearIdx);
            
            % Ensure y does not exceed available incorrect pairs
            if y > length(incorrectPairsLinearIdx)
                y_actual = length(incorrectPairsLinearIdx);
            else
                y_actual = y;
            end
            % Randomly select y_actual incorrect pairs to forbid
            forbiddenLinearIdx = incorrectPairsLinearIdx(randperm(length(incorrectPairsLinearIdx), y_actual));
            
            % Convert linear indices back to subscripts (row, col) and set allowed to false
            [r, c] = ind2sub([n, n], forbiddenLinearIdx);
            for pair_idx = 1:length(r)
                allowed(r(pair_idx), c(pair_idx)) = false;
            end
            perms = allValidPerms(allowed);
            if isempty(perms)
                continue
            end
            chosen = perms(randi(size(perms,1)), :);
            matchCount = matchCount + sum(chosen == 1:n);
            validTrials = validTrials + 1;
        end
        if validTrials > 0
            E_global(yi, n) = matchCount / validTrials;
        else
            E_global(yi, n) = NaN;
        end
    end
end
%% Part 3.2: Eliminate y Global Pairs
% This part simulates eliminating y options globally without guaranteeing correct matches are kept.
E_global_wrong = zeros(length(y_list), maxN);
for n = 1:maxN
    for yi = 1:length(y_list)
        y = y_list(yi);
        if y > n^2 % If y is greater than all possible pairs
            E_global_wrong(yi, n) = NaN;
            continue
        end
        matchCount = 0;
        validTrials = 0;
        for t = 1:trials
            allowed = true(n);
            
            % Ensure y does not exceed total possible pairs
            if y > n^2
                y_actual = n^2;
            else
                y_actual = y;
            end
            forbiddenLinearIdx = randperm(n^2, y_actual); % May forbid correct pairs
            
            [r, c] = ind2sub([n, n], forbiddenLinearIdx);
            for pair_idx = 1:length(r)
                allowed(r(pair_idx), c(pair_idx)) = false;
            end
            
            perms = allValidPerms(allowed);
            if isempty(perms)
                continue
            end
            chosen = perms(randi(size(perms,1)), :);
            matchCount = matchCount + sum(chosen == 1:n);
            validTrials = validTrials + 1;
        end
        if validTrials > 0
            E_global_wrong(yi, n) = matchCount / validTrials;
        else
            E_global_wrong(yi, n) = NaN;
        end
    end
end
%% Part 4: Scaled Global Elimination to Match Per-Word Elimination
% A third graph that is (3) scaled to (2) as in the amount of interval of increasing pair eliminations
% matches the actual amount of eliminations being performed by k-per word eliminations.
for n = 1:maxN
    for ki = 1:length(k_list)
        k = k_list(ki);
        
        % The formula for the relation:
        % When k wrong options are eliminated per word for n words,
        % the total number of *wrong* pairs eliminated is n * k.
        % This is because for each of the 'n' words, 'k' distinct wrong definitions are removed.
        % So, y = n * k.
        y_scaled = n * k; 
        y_from_k(ki, n) = y_scaled; % Store for display/checking
        
        % Check for invalid y_scaled (cannot exceed n^2 - n incorrect pairs)
        maxIncorrectPairs = n^2 - n;
        if y_scaled > maxIncorrectPairs || k >= n
            E_scaled_k_elim(ki, n) = NaN;
            continue
        end
        
        matchCount = 0;
        validTrials = 0;
        for t = 1:trials
            allowed = true(n);
            
            allIdx = 1:n^2;
            correctPairsLinearIdx = sub2ind([n n], 1:n, 1:n);
            incorrectPairsLinearIdx = setdiff(allIdx, correctPairsLinearIdx);
            
            % Ensure y_scaled does not exceed available incorrect pairs
            if y_scaled > length(incorrectPairsLinearIdx)
                current_y_actual = length(incorrectPairsLinearIdx);
            else
                current_y_actual = y_scaled;
            end
            forbiddenLinearIdx = incorrectPairsLinearIdx(randperm(length(incorrectPairsLinearIdx), current_y_actual));
            
            [r, c] = ind2sub([n, n], forbiddenLinearIdx);
            for pair_idx = 1:length(r)
                allowed(r(pair_idx), c(pair_idx)) = false;
            end
            
            perms = allValidPerms(allowed);
            if isempty(perms)
                continue
            end
            chosen = perms(randi(size(perms,1)), :);
            matchCount = matchCount + sum(chosen == 1:n);
            validTrials = validTrials + 1;
        end
        if validTrials > 0
            E_scaled_k_elim(ki, n) = matchCount / validTrials;
        else
            E_scaled_k_elim(ki, n) = NaN;
        end
    end
end
%% Visualization
% TODO (5): Plot expected correct guesses for each method
% Per-Word Elimination Plot
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
for ki = 1:length(k_list)
    lineData = E_perword(ki,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5, 'DisplayName', sprintf('Elim %d/word', k_list(ki)));
    
    for n_idx = 1:maxN % Iterate over each data point
        if ~isnan(lineData(n_idx))
            baseline_at_n = E_random(n_idx);
            actual_at_n = lineData(n_idx);
            max_possible_gain = n_idx - baseline_at_n;
            
            if max_possible_gain > 0
                imp = 100 * (actual_at_n - baseline_at_n) / max_possible_gain;
            else
                imp = 0;
            end
            txt = sprintf(' %+0.1f%%', imp);
            text(n_idx + 0.2, actual_at_n, txt, 'FontSize', 8);
        end
    end
end
text(maxN + 0.2, E_random(maxN), '  0%', 'FontSize', 8, 'Color', 'k');
xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n (Number of Words)');
ylabel('Expected Correct Matches');
title('Effect of Per-Word Elimination on Accuracy (Guaranteed Correct)');
legend('Location', 'southwest');
grid on;
hold off;

% Global Pair Elimination Plot
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
for yi = 1:length(y_list)
    lineData = E_global(yi,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5, 'DisplayName', sprintf('Elim %d pairs', y_list(yi)));
    
    for n_idx = 1:maxN % Iterate over each data point
        if ~isnan(lineData(n_idx))
            baseline_at_n = E_random(n_idx);
            actual_at_n = lineData(n_idx);
            max_possible_gain = n_idx - baseline_at_n;
            
            if max_possible_gain > 0
                imp = 100 * (actual_at_n - baseline_at_n) / max_possible_gain;
            else
                imp = 0;
            end
            txt = sprintf(' %+0.1f%%', imp);
            text(n_idx + 0.2, actual_at_n, txt, 'FontSize', 8);
        end
    end
end
text(maxN + 0.2, E_random(maxN), '  0%', 'FontSize', 8, 'Color', 'k');
xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n (Number of Words)');
ylabel('Expected Correct Matches');
title('Effect of Global Pair Elimination on Accuracy (Guaranteed Correct)');
legend('Location', 'southwest');
grid on;
hold off;

% Per-Word Elimination Plot (Possibly Incorrect)
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
for ki = 1:length(k_list)
    lineData = E_perword_wrong(ki,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5, 'DisplayName', sprintf('Elim %d/word', k_list(ki)));
    
    for n_idx = 1:maxN % Iterate over each data point
        if ~isnan(lineData(n_idx))
            baseline_at_n = E_random(n_idx);
            actual_at_n = lineData(n_idx);
            max_possible_gain = n_idx - baseline_at_n;
            
            if max_possible_gain > 0
                imp = 100 * (actual_at_n - baseline_at_n) / max_possible_gain;
            else
                imp = 0;
            end
            txt = sprintf(' %+0.1f%%', imp);
            text(n_idx + 0.2, actual_at_n, txt, 'FontSize', 8);
        end
    end
end
text(maxN + 0.2, E_random(maxN), '  0%', 'FontSize', 8, 'Color', 'k');
xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n (Number of Words)');
ylabel('Expected Correct Matches');
title('Effect of Per-Word Elimination (May Forbid Correct)');
legend('Location', 'southwest');
grid on;
hold off;

% Global Pair Elimination (Possibly Incorrect)
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
for yi = 1:length(y_list)
    lineData = E_global_wrong(yi,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5, 'DisplayName', sprintf('Elim %d pairs', y_list(yi)));
    
    for n_idx = 1:maxN % Iterate over each data point
        if ~isnan(lineData(n_idx))
            baseline_at_n = E_random(n_idx);
            actual_at_n = lineData(n_idx);
            max_possible_gain = n_idx - baseline_at_n;
            
            if max_possible_gain > 0
                imp = 100 * (actual_at_n - baseline_at_n) / max_possible_gain;
            else
                imp = 0;
            end
            txt = sprintf(' %+0.1f%%', imp);
            text(n_idx + 0.2, actual_at_n, txt, 'FontSize', 8);
        end
    end
end
text(maxN + 0.2, E_random(maxN), '  0%', 'FontSize', 8, 'Color', 'k');
xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n (Number of Words)');
ylabel('Expected Correct Matches');
title('Effect of Global Pair Elimination (May Forbid Correct)');
legend('Location', 'southwest');
grid on;
hold off;

% Dynamic Global Elimination
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
for ki = 1:length(k_list)
    lineData = E_scaled_k_elim(ki,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5, 'DisplayName', sprintf('Scaled Elim (k=%d, y=n*k)', k_list(ki)));
    
    for n_idx = 1:maxN % Iterate over each data point
        if ~isnan(lineData(n_idx))
            baseline_at_n = E_random(n_idx);
            actual_at_n = lineData(n_idx);
            max_possible_gain = n_idx - baseline_at_n;
            
            if max_possible_gain > 0
                imp = 100 * (actual_at_n - baseline_at_n) / max_possible_gain;
            else
                imp = 0;
            end
            txt = sprintf(' %+0.1f%%', imp);
            text(n_idx + 0.2, actual_at_n, txt, 'FontSize', 8);
        end
    end
end
text(maxN + 0.2, E_random(maxN), '  0%', 'FontSize', 8, 'Color', 'k');
xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n (Number of Words)');
ylabel('Expected Correct Matches');
title('Dynamic Global Elimination');
legend('Location', 'southwest');
grid on;
hold off;

% Global Elimination with CONSTANT Absolute Eliminations (k-scaled Global Pair Elimination)
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2, 'DisplayName', 'Random');
n_values_for_y_selection = [2, 3, 4]; 
constant_y_values_to_plot = zeros(1, length(n_values_for_y_selection));
for i = 1:length(n_values_for_y_selection)
    n_for_y = n_values_for_y_selection(i);
    constant_y_values_to_plot(i) = n_for_y * (n_for_y - 1);
end
constant_y_values_to_plot = unique(sort(constant_y_values_to_plot));
constant_y_values_to_plot = constant_y_values_to_plot(constant_y_values_to_plot > 0);
for y_val_to_plot_idx = 1:length(constant_y_values_to_plot)
    y_plot_val = constant_y_values_to_plot(y_val_to_plot_idx);
    
    yi_idx = find(y_list == y_plot_val, 1);
    
    if ~isempty(yi_idx) && yi_idx <= size(E_global, 1)
        lineData = E_global(yi_idx, :);
        plot(1:maxN, lineData, 'LineWidth', 1.5, 'DisplayName', sprintf('Constant Elim %d pairs', y_plot_val));
        
        for n_idx = 1:maxN % Iterate over each data point
            if ~isnan(lineData(n_idx))
                baseline_at_n = E_random(n_idx);
                actual_at_n = lineData(n_idx);
                max_possible_gain = n_idx - baseline_at_n;
                
                if max_possible_gain > 0
                    imp = 100 * (actual_at_n - baseline_at_n) / max_possible_gain;
                else
                    imp = 0;
                end
                txt = sprintf(' %+0.1f%%', imp);
                text(n_idx + 0.2, actual_at_n, txt, 'FontSize', 8);
            end
        end
    else
        warning('Cannot plot for constant y=%d as it is not in y_list or E_global data.', y_plot_val);
    end
end
text(maxN + 0.2, E_random(maxN), '  0%', 'FontSize', 8, 'Color', 'k');
xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n (Number of Words)');
ylabel('Expected Correct Matches');
title('k-scaled Global Pair Elimination');
legend('Location', 'southwest');
grid on;
hold off;
%% Helper Function: Generate All Valid Permutations
% TODO (6): Implement function to return all valid permutations based on allowed matrix
function validPerms = allValidPerms(allowed)
    n = size(allowed, 1);
    if n == 0
        validPerms = [];
        return;
    end
    
    allP = perms(1:n);
    numPerms = size(allP, 1);
    
    mask = true(numPerms, 1);
    
    for i = 1:n
        current_column_indices = allP(:,i);
        linear_indices = sub2ind(size(allowed), repmat(i, numPerms, 1), current_column_indices);
        
        mask = mask & allowed(linear_indices);
    end
    validPerms = allP(mask,:);
end