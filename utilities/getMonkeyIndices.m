function idx = getMonkeyIndices(ids)
%GETMONKEYINDICES Per-monkey row indices for a list of session/unit IDs.
%
% idx = getMonkeyIndices(ids) accepts either:
%   - a table with a 'ses_ID' column (e.g., a subset from createDatSubset), or
%   - a cell array of session/unit ID strings (e.g., unique(dat.ses_ID)).
%
% Returns a struct with logical masks .An, .Ch, .Mi, the same length as
% the (unique'd) ID list, so they can be used directly to index any
% array/table with one row per session or unit, e.g.:
%
%   idx = getMonkeyIndices(dat);
%   LSF_slope_An = fits(idx.An, 3, 1);
%
% Monkey identity is matched by substring on the ID itself (e.g., 'An' in
% 'AnDR_20210614_Unit1_BNP'), rather than hardcoded numeric ranges.
% Hardcoded ranges (e.g., An = 1:55; Ch = 56:68; Mi = 69:155;) differ
% across data subsets (behavioral vs. neural vs. pupil) because each
% subset excludes a different number of sessions/units, so they must be
% recomputed for every subset — deriving them from the IDs directly
% avoids having to keep several hardcoded triples in sync.
%
% NOTE: Monkey Mi's session/unit IDs use the raw code "MM" (e.g.
% "MMDR_20250211_02_Unit1_BNP"), not "Mi" — confirmed against
% mergedTable_proc.ses_ID. "Mi" never appears in any real ID. The struct
% field below is still named .Mi (matching the paper's monkey label and
% every other script's variable naming), but it is matched via "MM".
%
% idx.names = {'An','Ch','Mi'} is included for convenience when looping
% over monkeys.

if istable(ids)
    ids = unique(ids.ses_ID);
end

idx.An = contains(ids, 'An');
idx.Ch = contains(ids, 'Ch');
idx.Mi = contains(ids, 'MM');

assert(~any((idx.An & idx.Ch) | (idx.An & idx.Mi) | (idx.Ch & idx.Mi)), ...
    'getMonkeyIndices:ambiguousID', ...
    'Some session/unit IDs matched more than one monkey substring.');

idx.names = {'An', 'Ch', 'Mi'};

end
