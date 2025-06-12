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
trials = 10000;     % Number of random permutations to average over
k_list = 0:2;       % Number of eliminated options per word
y_list = 0:6;       % Number of globally eliminated word-definition pairs

% Store results in matrices:
E_random = zeros(1, maxN);
E_perword = zeros(length(k_list), maxN);
E_global = zeros(length(y_list), maxN);

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
            E_perword(ki, n) = NaN;  % or 0, or skip entirely
            continue
        end

        matchCount = 0;
        validTrials = 0;

        for t = 1:trials
            allowed = true(n); % allowed(i,j) = true if word i can be matched with def j
            for i = 1:n
                choices = setdiff(1:n, i); % exclude the correct match
                forbidden = choices(randperm(length(choices), k));
                % forbidden = randperm(n, k);
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
            E_perword(ki, n) = matchCount / validTrials;
        else
            E_perword(ki, n) = NaN;
        end
    end
end

%% Part 2.2: Eliminate k Options per Word
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
                forbidden = randperm(n, k);  % no guarantee correct match is kept
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

        maxForbiddenPairs = n^2 - n; % can't forbid correct pairs

        if y > maxForbiddenPairs
            E_global(yi, n) = NaN;
            continue
        end

        matchCount = 0;
        validTrials = 0;

        for t = 1:trials
            allowed = true(n);
            % forbiddenPairs = randperm(n^2, y);

            [W, D] = ndgrid(1:n, 1:n);
            allPairs = [W(:), D(:)];
            correctPairs = sub2ind([n n], 1:n, 1:n);
            allIdx = 1:n^2;
            allowedIdx = setdiff(allIdx, correctPairs);
            forbiddenPairs = allowedIdx(randperm(length(allowedIdx), y));
            
            for idx = forbiddenPairs
                i = mod(idx-1, n) + 1;
                j = floor((idx-1)/n) + 1;
                allowed(i,j) = false;
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
E_global_wrong = zeros(length(y_list), maxN);

for n = 1:maxN
    for yi = 1:length(y_list)
        y = y_list(yi);
        if y > n^2
            E_global_wrong(yi, n) = NaN;
            continue
        end

        matchCount = 0;
        validTrials = 0;

        for t = 1:trials
            allowed = true(n);
            forbiddenIdx = randperm(n^2, y);
            for idx = forbiddenIdx
                i = mod(idx-1, n) + 1;
                j = floor((idx-1)/n) + 1;
                allowed(i,j) = false;
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

%% Visualization
% TODO (5): Plot expected correct guesses for each method

% Per-Word Elimination Plot
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2);

for ki = 1:length(k_list)
    lineData = E_perword(ki,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5);
    
    % Calculate relative improvement (%) at last n
    imp = 100 * (lineData(end) - E_random(end)) / E_random(end);
    txt = sprintf(' %+0.1f%%', imp);
    text(maxN + 0.2, lineData(end), txt, 'FontSize', 8);
end

% Annotate random baseline (improvement is zero)
text(maxN + 0.2, E_random(end), '  0%', 'FontSize', 8, 'Color', 'k');

xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n');
ylabel('Expected correct guesses');
title('Effect of Per-Word Elimination on Accuracy');
legend([{ 'Random' }, compose('Elim %d/word', k_list)], 'Location', 'southwest');
hold off;


% Global Pair Elimination Plot
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2);

for yi = 1:length(y_list)
    lineData = E_global(yi,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5);
    
    % Calculate relative improvement (%) at last n
    imp = 100 * (lineData(end) - E_random(end)) / E_random(end);
    txt = sprintf(' %+0.1f%%', imp);
    text(maxN + 0.2, lineData(end), txt, 'FontSize', 8);
end

% Annotate random baseline (improvement is zero)
text(maxN + 0.2, E_random(end), '  0%', 'FontSize', 8, 'Color', 'k');

xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n');
ylabel('Expected correct guesses');
title('Effect of Global Pair Elimination on Accuracy');
legend([{ 'Random' }, compose('Elim %d pairs', y_list)], 'Location', 'southwest');
hold off;

% Per-Word Elimination (Incorrect) Plot
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2);

for ki = 1:length(k_list)
    lineData = E_perword_wrong(ki,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5);

    imp = 100 * (lineData(end) - E_random(end)) / E_random(end);
    txt = sprintf(' %+0.1f%%', imp);
    text(maxN + 0.2, lineData(end), txt, 'FontSize', 8);
end

text(maxN + 0.2, E_random(end), '  0%', 'FontSize', 8, 'Color', 'k');

xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n');
ylabel('Expected correct guesses');
title('Effect of Incorrect Per-Word Elimination');
legend([{ 'Random' }, compose('Elim %d/word', k_list)], 'Location', 'southwest');
hold off;

% Global Pair Elimination (Incorrect) Plot
figure;
hold on;
plot(1:maxN, E_random, 'k--', 'LineWidth', 2);

for yi = 1:length(y_list)
    lineData = E_global_wrong(yi,:);
    plot(1:maxN, lineData, 'LineWidth', 1.5);

    imp = 100 * (lineData(end) - E_random(end)) / E_random(end);
    txt = sprintf(' %+0.1f%%', imp);
    text(maxN + 0.2, lineData(end), txt, 'FontSize', 8);
end

text(maxN + 0.2, E_random(end), '  0%', 'FontSize', 8, 'Color', 'k');

xticks(1:maxN);
xlim([1 maxN + 1.5]);
xlabel('n');
ylabel('Expected correct guesses');
title('Effect of Incorrect Global Pair Elimination');
legend([{ 'Random' }, compose('Elim %d pairs', y_list)], 'Location', 'southwest');
hold off;

%% Helper Function: Generate All Valid Permutations
% TODO (6): Implement function to return all valid permutations based on allowed matrix
function validPerms = allValidPerms(allowed)
    n = size(allowed, 1);
    allP = perms(1:n);
    mask = true(size(allP,1),1);
    for i = 1:n
        mask = mask & allowed(sub2ind(size(allowed), i*ones(size(allP,1),1), allP(:,i)));
    end
    validPerms = allP(mask,:);
end
