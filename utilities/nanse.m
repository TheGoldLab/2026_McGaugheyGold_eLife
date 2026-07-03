function se = nanse(x,dim)
% nanse(x,dim) returns the std error of the matix x 
%   along the dimension dim

% 1/15/97 mns wrote it

% if ~isempty(x)
%     if nargin < 2
%         dim = [];
%     end
%     %se = nanstd(x,[],dim) ./ sqrt(sum(~isnan(x),dim));
% else
%   se = nan;
% end
%


% 1/15/2026 kdm edited for MATLAB2025a

if ~isempty(x)
    if nargin < 2 || isempty(dim)
        % Find first non-singleton dimension (mimics old behavior)
        dim = find(size(x) > 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end
    se = std(x, [], dim, 'omitnan') ./ sqrt(sum(~isnan(x), dim));
else
    se = nan;
end