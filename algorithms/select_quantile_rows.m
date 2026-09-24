function [selTbl, finalIdx] = select_quantile_rows(T, nSelect, topCandidates, rankCol)
%SELECT_QUANTILE_ROWS Pick rows spread evenly across WalshScore quantiles.
%   Selects NSELECT rows from the TOPCANDIDATES highest-ranked rows of table
%   T (ranked by RANKCOL, descending), spacing the selection evenly across
%   WalshScore quantiles so the chosen codes span the sequency range rather
%   than clustering at one end.
%
%   [SELTBL, FINALIDX] = SELECT_QUANTILE_ROWS(T) selects 16 rows from the
%   top 64 candidates ranked by T.MeanAccuracy.
%
%   [SELTBL, FINALIDX] = SELECT_QUANTILE_ROWS(T, NSELECT, TOPCANDIDATES, RANKCOL)
%   overrides the defaults. Pass RANKCOL = 'PctIDCorrect' to rank by the
%   Eq.(4) codebook-identification metric (from the corrected
%   hadfixed9_2 pipeline) instead of the legacy per-chip MeanAccuracy.
%
%   T must contain columns: OriginalRow, WalshScore, and RANKCOL.
%   FINALIDX are row indices into the ORIGINAL table T (not into SELTBL),
%   sorted so that FINALIDX(i) corresponds to SELTBL row i.
%
%   Example:
%       T = readtable('fixed_hadamard_results_snr_10dB.csv');
%       [selTbl, finalIdx] = select_quantile_rows(T, 16, 64, 'PctIDCorrect');
%
%   See also SELECT_GREEDY_ROWS.

    if nargin < 2 || isempty(nSelect)
        nSelect = 16;
    end
    if nargin < 3 || isempty(topCandidates)
        topCandidates = 64;
    end
    if nargin < 4 || isempty(rankCol)
        rankCol = 'MeanAccuracy';
    end

    requiredCols = {'OriginalRow', 'WalshScore', rankCol};
    if ~all(ismember(requiredCols, T.Properties.VariableNames))
        error('select_quantile_rows:missingColumns', ...
            'Table T must contain columns: %s', strjoin(requiredCols, ', '));
    end

    ws = T.WalshScore;
    if ~isnumeric(ws); ws = str2double(string(ws)); end
    rankVal = T.(rankCol);
    if ~isnumeric(rankVal); rankVal = str2double(string(rankVal)); end

    % Select top candidates by the ranking column (descending)
    [~, orderAcc] = sort(rankVal, 'descend', 'MissingPlacement', 'last');
    candIdx = orderAcc(1:min(topCandidates, numel(orderAcc)));
    ws_cand = ws(candIdx);
    nAvail = numel(ws_cand);

    if nAvail <= nSelect
        selIdxLocal = (1:nAvail)';
    else
        % Quantile / even-spacing selection on WalshScore among candidates
        X = ws_cand(:);
        q = linspace(0, 1, nSelect + 1);
        qcenters = (q(1:end-1) + q(2:end)) / 2;   % bin-center quantiles avoid exact 0/1
        targetVals = quantile(X, qcenters);

        selIdxLocal = zeros(nSelect, 1);
        used = false(nAvail, 1);
        for b = 1:nSelect
            [~, orderIdx] = sort(abs(X - targetVals(b)), 'ascend');
            pick = orderIdx(find(~used(orderIdx), 1));
            if isempty(pick)
                pick = find(~used, 1); % fallback
            end
            selIdxLocal(b) = pick;
            used(pick) = true;
        end

        % Defensive fill if duplicates or fewer than required were found
        selIdxLocal = unique(selIdxLocal, 'stable');
        if numel(selIdxLocal) < nSelect
            remaining = setdiff((1:nAvail)', selIdxLocal);
            fillIdx = round(linspace(1, numel(remaining), nSelect - numel(selIdxLocal)));
            selIdxLocal = [selIdxLocal; remaining(fillIdx(:))];
        end
    end

    % Map local candidate indices back to table row indices, then sort by
    % WalshScore so FINALIDX and SELTBL rows correspond 1:1 in the same order.
    finalIdx = candIdx(selIdxLocal);
    [~, sortOrder] = sort(ws(finalIdx), 'ascend');
    finalIdx = finalIdx(sortOrder);
    selTbl = T(finalIdx, :);
end
