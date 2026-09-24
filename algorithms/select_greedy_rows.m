function [selTbl, finalIdx] = select_greedy_rows(T, nSelect, topCandidates, rankCol)
%SELECT_GREEDY_ROWS Pick rows via greedy farthest-point sampling on WalshScore.
%   Selects NSELECT rows from the TOPCANDIDATES highest-ranked rows of table
%   T (ranked by RANKCOL, descending), starting from the top-ranked candidate
%   and repeatedly adding whichever remaining candidate is farthest (in
%   WalshScore) from everything already chosen. This tends to spread
%   selections more aggressively toward the extremes than quantile binning.
%
%   [SELTBL, FINALIDX] = SELECT_GREEDY_ROWS(T) selects 16 rows from the top
%   64 candidates ranked by T.MeanAccuracy.
%
%   [SELTBL, FINALIDX] = SELECT_GREEDY_ROWS(T, NSELECT, TOPCANDIDATES, RANKCOL)
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
%       [selTbl, finalIdx] = select_greedy_rows(T, 16, 64, 'PctIDCorrect');
%
%   See also SELECT_QUANTILE_ROWS.

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
        error('select_greedy_rows:missingColumns', ...
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
        % Greedy farthest-point sampling on WalshScore among candidates
        selLocal = zeros(nSelect, 1);
        selLocal(1) = 1; % start with the highest-ranked candidate
        distToSel = abs(ws_cand - ws_cand(selLocal(1)));
        for ii = 2:nSelect
            [~, nextI] = max(distToSel);
            selLocal(ii) = nextI;
            distToSel = min(distToSel, abs(ws_cand - ws_cand(nextI)));
        end
        selIdxLocal = unique(selLocal, 'stable');
        if numel(selIdxLocal) < nSelect
            remaining = setdiff((1:nAvail)', selIdxLocal);
            while numel(selIdxLocal) < nSelect
                dmin = min(abs(ws_cand(remaining) - ws_cand(selIdxLocal)'), [], 2);
                [~, mi] = max(dmin);
                selIdxLocal(end+1) = remaining(mi); %#ok<AGROW>
                remaining(mi) = [];
            end
        end
    end

    % Map local candidate indices back to table row indices, then sort by
    % WalshScore so FINALIDX and SELTBL rows correspond 1:1 in the same order.
    finalIdx = candIdx(selIdxLocal);
    [~, sortOrder] = sort(ws(finalIdx), 'ascend');
    finalIdx = finalIdx(sortOrder);
    selTbl = T(finalIdx, :);
end
