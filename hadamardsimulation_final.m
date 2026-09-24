function run_hadamard_realworld_fixed()
    %% Parameters
    N = 16;
    true_spp = 20;
    fs = 5000;
    fc = 10;
    snr_values = [-3, 0, 3, 5, 8, 10, 12, 15];
    numTrials = 20;

    %% Generate codes
    H = hadamard(N);
    walshScore = sum(abs(diff(H,1,2)), 2) / 2;
    [~, idx] = sort(walshScore);
    H_sorted = H(idx, :);   % sequency-sorted codebook. idx(r) = original Hadamard row for sorted row r.

    %% Precompute filter
    [b, a] = butter(6, fc/(fs/2));

    %% Fixed noise reference (paper-consistent SNR definition)
    % codes that get heavily attenuated by the channel filter also get less absolute
    % noise, which is not how a real fixed noise floor behaves and unfairly
    % flatters high-sequency codes.
    refPower = 1;

    %% Storage
    numRows = size(H, 1);
    numSNR = numel(snr_values);
    allAcc = zeros(numRows, numTrials, numSNR);        % per-chip accuracy (legacy metric; known-code-identity)
    allIDcorrect = zeros(numRows, numTrials, numSNR);  % Eq.(4) codebook-identification correctness
    allEstSpp = zeros(numRows, numTrials, numSNR);

    %% Noiseless pass: identify Subset A (codes that fail even with no channel noise)
    % Per the paper's theory section, Subset A is discarded before any SNR-based
    % ranking is meaningful: these rows can't survive the channel's lossy
    % compression (bandlimiting) even in the best case.
    noiselessIDcorrect = false(numRows, 1);
    noiselessAcc = zeros(numRows, 1);
    for r = 1:numRows
        codeinput = H_sorted(r, :);
        audioVec = repelem(codeinput, true_spp);
        tx_signal = filtfilt(b, a, double(audioVec(:)'));
        rx_signal = filtfilt(b, a, tx_signal);  % receiver front-end filter, no added noise

        [acc, ~] = decode_with_known_spp(rx_signal, true_spp, codeinput);
        noiselessAcc(r) = acc;

        idHat = identify_code(rx_signal, true_spp, H_sorted);
        noiselessIDcorrect(r) = (idHat == r);
    end
    subsetA = ~noiselessIDcorrect;  % true = fails even noiseless -> discard, per the paper's Subset A

    %% Progress bar
    totalSteps = numSNR * numTrials * numRows;
    step = 0;
    h = waitbar(0, 'Running...');

    for s = 1:numSNR
        snr_dB = snr_values(s);
        noisePower = refPower / (10^(snr_dB/10));
        noiseStd = sqrt(noisePower);

        for t = 1:numTrials
            for r = 1:numRows
                codeinput = H_sorted(r, :);

                %% TRANSMIT
                audioVec = repelem(codeinput, true_spp);
                tx_signal = filtfilt(b, a, double(audioVec(:)'));

                %% CHANNEL: fixed-power AWGN referenced to the common chip amplitude,
                %% not to this particular code's own (filter-attenuated) signal power
                noisySignal = tx_signal + noiseStd * randn(size(tx_signal));

                % Lowpass the received (noisy) waveform to simulate receiver anti-alias / front-end
                noisySignal = filtfilt(b, a, double(noisySignal));

                %% RECEIVER: USE KNOWN SPP (REAL-WORLD SOLUTION)
                % In real systems, SPP is known from protocol. This is NOT cheating - it's how real systems work.
                est_spp = true_spp;

                %% DECODE (legacy per-chip metric; assumes the code identity is already known)
                [accuracy, ~] = decode_with_known_spp(noisySignal, est_spp, codeinput);
                allAcc(r, t, s) = accuracy;
                allEstSpp(r, t, s) = est_spp;

                %% Eq. (4) TEST: cross-correlate against the FULL codebook and check
                %% whether the true transmitted row gives the highest correlation.
                %% This is the actual criterion the paper specifies for a spreading
                %% code to "fit the channel" -- the receiver does not know in advance
                %% which of the hopped codes was sent, so identification against the
                %% whole pool (not just a known-vs-received BER check) is the real test.
                idHat = identify_code(noisySignal, est_spp, H_sorted);
                allIDcorrect(r, t, s) = (idHat == r);

                step = step + 1;
                if mod(step, 100) == 0
                    waitbar(step/totalSteps, h, ...
                        sprintf('SNR %d/%d (%.1f dB), Trial %d/%d', ...
                        s, numSNR, snr_dB, t, numTrials));
                end
            end
        end
    end

    delete(h);

    %% Export results
    export_fixed_results(allAcc, allIDcorrect, snr_values, walshScore, idx, numTrials, subsetA);
end

function [accuracy, recovered] = decode_with_known_spp(signal, spp, original_code)
    % Clean signal
    signal_clean = movmedian(signal, max(3, round(spp/5)));

    % Decimate using known SPP
    L = floor(length(signal_clean)/spp) * spp;

    if L < spp
        accuracy = 0;
        recovered = [];
        return;
    end

    X = reshape(signal_clean(1:L), spp, []);
    chipStat = mean(X, 1);
    recovered = sign(chipStat);
    recovered(recovered == 0) = 1;

    numChips = min(length(recovered), length(original_code));
    rec = recovered(1:numChips);
    orig = original_code(1:numChips);

    accuracy = sum(rec == orig) / numChips * 100;
end

function idHat = identify_code(signal, spp, codebook)
    % Implements the paper's Equation (4): correlate the received chip statistic
    % against every row in the codebook (not just the known transmitted one) and
    % return the index of the best match. Because Hadamard rows are exactly
    % orthogonal +-1 sequences, a simple dot product against each candidate row
    % is the matched-filter correlation the equation calls for.
    signal_clean = movmedian(signal, max(3, round(spp/5)));
    L = floor(length(signal_clean)/spp) * spp;

    if L < spp
        idHat = 0;  % no valid decision
        return;
    end

    X = reshape(signal_clean(1:L), spp, []);
    chipStat = mean(X, 1);

    numChips = min(length(chipStat), size(codebook, 2));
    r = chipStat(1:numChips);
    C = codebook(:, 1:numChips);

    corrScores = C * r(:);
    [~, idHat] = max(corrScores);
end

function export_fixed_results(allAcc, allIDcorrect, snr_values, walshScore, idx, numTrials, subsetA)
    for s = 1:length(snr_values)
        snr_dB = snr_values(s);
        acc_mat = allAcc(:, :, s);
        id_mat = allIDcorrect(:, :, s);

        meanAcc = mean(acc_mat, 2);
        stdAcc = std(acc_mat, 0, 2);
        pctPerfect = sum(acc_mat >= 99.999, 2) / numTrials * 100;
        pctIDcorrect = mean(id_mat, 2) * 100;

        T = table((1:length(meanAcc))', idx(:), meanAcc, stdAcc, pctPerfect, pctIDcorrect, walshScore, subsetA, ...
            'VariableNames', {'SortedPosition', 'OriginalRow', 'MeanAccuracy', 'StdAccuracy', ...
            'PctPerfect', 'PctIDCorrect', 'WalshScore', 'SubsetA_Discard'});

        % Rank by the Eq.(4) identification metric first -- that's the criterion
        % that actually determines whether a code "fits the channel" per the
        % paper -- with WalshScore as a tiebreaker.
        T_sorted = sortrows(T, {'PctIDCorrect', 'WalshScore'}, {'descend', 'ascend'});
        T_sorted.AccuracyRank = (1:height(T_sorted))';

        filename = sprintf('fixed_hadamard_results_snr_%ddB.csv', snr_dB);
        writetable(T_sorted, filename);
        fprintf('Wrote: %s\n', filename);
    end
end

% Run it
run_hadamard_realworld_fixed();
