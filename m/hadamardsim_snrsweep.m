% run_hadamard_trials_with_progress_snr_sweep.m
N = 16;
H = hadamard(N);
samplesPerChip = 20;
fs = 5000;
fc = 10;

snr_values = [-3, 0, 3, 5, 8, 10, 12, 15];

numTrials = 20;
numRows = size(H,1);
numSNR = numel(snr_values);

allAcc = zeros(numRows, numTrials, numSNR);
allEstSpp = zeros(numRows, numTrials, numSNR);

totalSteps = numSNR * numTrials * numRows;
step = 0;

% waitbar
h = waitbar(0, 'Starting...', 'Name', 'Hadamard Trials Progress (SNR Sweep)');
setappdata(h, 'canceling', 0);
uicontrol('Parent', h, 'Style', 'pushbutton', 'String', 'Cancel', ...
    'Units', 'normalized', 'Position', [0.82 0.02 0.15 0.07], ...
    'Callback', @(~,~) setappdata(h,'canceling',1));

try
    for s = 1:numSNR
        snr_dB = snr_values(s);
        for t = 1:numTrials
            sequency = sum(abs(diff(H,1,2)), 2) / 2;
            [~, idx] = sort(sequency);
            H_sorted = H(idx, :);

            accuracies = zeros(numRows,1);
            estSpp = zeros(numRows,1);

            for r = 1:numRows
                % Check cancel
                if getappdata(h, 'canceling')
                    warning('User canceled operation.');
                    error('CANCELED');
                end

                codeinput = H_sorted(r,:);
                analogSignal0 = repelem(codeinput, samplesPerChip);
                audioVec = analogSignal0(:);
                if max(abs(audioVec))>0
                    audioVec = audioVec / max(abs(audioVec));
                end
                tmpName = sprintf('temp_lossy_s%d_t%d_r%d.wav', s, t, r);
                audiowrite(tmpName, audioVec, fs);
                [audioLossy, fs_read] = audioread(tmpName);
                if fs_read ~= fs
                    audioLossy = resample(audioLossy, fs, fs_read);
                end
                audioLossy = audioLossy(:)'; delete(tmpName);
                analogSignal = audioLossy;

                filteredSignal = lowpass(analogSignal, fc, fs);
                noisySignal = awgn(filteredSignal, snr_dB, 'measured');

                y = noisySignal(:)';
                y_smooth = movmedian(y, max(3, round(length(y)/1000)));
                digital_est = sign(y_smooth); digital_est(digital_est==0)=1;
                d = [true, diff(digital_est) ~= 0];
                idx_ch = find(d); idx_ch = [idx_ch, length(digital_est)+1];
                runLengths = diff(idx_ch); runValues = digital_est(idx_ch(1:end-1));
                if isempty(runLengths)
                    est = 1;
                else
                    uniqueLens = unique(runLengths);
                    counts = histc(runLengths, uniqueLens);
                    [~, order] = sort(counts, 'descend');
                    K = min(3, numel(order));
                    estCandidates = uniqueLens(order(1:K));
                    est = max(1, round(median(estCandidates)));
                    minRun = max(1, round(0.25 * est));
                    if any(runLengths < minRun)
                        mergedVals = runValues; mergedLens = runLengths; i = 1;
                        while i <= length(mergedLens)
                            if mergedLens(i) < minRun && length(mergedLens) > 1
                                if i == 1
                                    mergedLens(2) = mergedLens(2) + mergedLens(1);
                                    mergedVals(1) = []; mergedLens(1) = [];
                                else
                                    mergedLens(i-1) = mergedLens(i-1) + mergedLens(i);
                                    mergedVals(i) = []; mergedLens(i) = []; i = i-1;
                                end
                            end
                            i = i + 1;
                        end
                        runLengths = mergedLens;
                        uniqueLens = unique(runLengths);
                        counts = histc(runLengths, uniqueLens);
                        [~, order] = sort(counts, 'descend'); K = min(3, numel(order));
                        est = max(1, round(median(uniqueLens(order(1:K)))));
                    end
                end

                L = floor(length(y)/est)*est;
                X = reshape(y(1:L), est, []);
                chipStat = mean(X,1);
                recoveredChips = sign(chipStat); recoveredChips(recoveredChips==0)=1;
                numChips = min(length(recoveredChips), length(codeinput));
                rec = recoveredChips(1:numChips); orig = codeinput(1:numChips);
                accuracies(r) = sum(rec==orig)/numChips*100;
                estSpp(r) = est;

                % update progress
                step = step + 1;
                waitbar(step/totalSteps, h, sprintf('SNR %d/%d (%.1f dB), Trial %d/%d, Row %d/%d (%.1f%%)', ...
                    s, numSNR, snr_dB, t, numTrials, r, numRows, step/totalSteps*100));
            end

            acc_by_orig = zeros(numRows,1);
            acc_by_orig(idx) = accuracies;
            allAcc(:,t,s) = acc_by_orig;
            allEstSpp(:,t,s) = estSpp(idx);
        end
    end
catch ME
    if ~strcmp(ME.message,'CANCELED')
        rethrow(ME);
    end
end

if ishandle(h)
    delete(h);
end

% For each SNR compute statistics and export CSVs
walshScore = sum(abs(diff(H,1,2)), 2) / 2;

for s = 1:numSNR
    snr_dB = snr_values(s);
    acc_mat = allAcc(:,:,s);

    meanAcc = mean(acc_mat, 2);
    stdAcc = std(acc_mat, 0, 2);
    pctPerfect = sum(acc_mat >= 99.999, 2) / numTrials * 100;

    T = table((1:numRows)', meanAcc, stdAcc, pctPerfect, walshScore, ...
        'VariableNames', {'OriginalRow','MeanAccuracy','StdAccuracy','PctPerfect','WalshScore'});

    T_sorted = sortrows(T, {'WalshScore','MeanAccuracy'}, {'ascend','descend'});
    T_sorted.WalshRank = (1:height(T_sorted))';

    % Display top rows briefly
    fprintf('SNR = %d dB: Top 5 rows by sequency/accuracy\n', snr_dB);
    disp(T_sorted(1:5,:));

    % Reorder for plotting (optional)
    orderByWalshThenAcc = T_sorted.OriginalRow;
    meanAcc_w = meanAcc(orderByWalshThenAcc);
    stdAcc_w  = stdAcc(orderByWalshThenAcc);
    pctPerfect_w = pctPerfect(orderByWalshThenAcc);

    figure;
    subplot(1,2,1);
    errorbar(1:numRows, meanAcc_w, stdAcc_w, 'o-');
    xlabel('Walsh Order (sorted by sequency; most accurate within group first)');
    ylabel('Accuracy (%)');
    title(sprintf('Mean ± STD (Sequency-Sorted)  SNR=%d dB', snr_dB));
    ylim([0 100]);

    subplot(1,2,2);
    bar(1:numRows, pctPerfect_w);
    xlabel('Walsh Order (sorted by sequency)');
    ylabel('% Trials ~100%');
    title(sprintf('Consistency (Sequency-Sorted)  SNR=%d dB', snr_dB));
    ylim([0 100]);



    outSummarySorted = sprintf('hadamard_row_stats_by_sequency_snr_%ddB.csv', snr_dB);
    writetable(T_sorted, outSummarySorted);

    accTable = array2table(acc_mat, 'VariableNames', strcat('Trial', string(1:numTrials)));
    accTable = addvars(accTable, (1:numRows)', 'Before', 1, 'NewVariableNames', 'OriginalRow');
    outAcc = sprintf('hadamard_accuracies_by_trial_snr_%ddB.csv', snr_dB);
    writetable(accTable, outAcc);



    fprintf('Wrote: %s\nWrote: %s\n', outSummarySorted, outAcc);
end