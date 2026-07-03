
function [mergedTable_proc_sub] = createDatSubset(mergedTable_proc, subset_code)
% createDatSubset Create data subsets based on available data types
%
% Extracts sessions with specific combinations of behavioral (B), neural (N),
% and pupil (P) data from the master data table. Data type is encoded in the
% session ID suffix (e.g., "SessionName_BNP" has all three data types).
%
% Inputs:
%   mergedTable_proc - Master data table containing all sessions
%   subset_code      - String specifying desired data subset:
%                      'B'   = Behavioral data only
%                      'BN'  = Behavioral + Neural data
%                      'BNP' = Behavioral + Neural + Pupil data
%                      'BP'  = Behavioral + Pupil data
%                      'N'   = Neural data (any session with neural recordings)
%                      'NP'  = Neural + Pupil data
%
% Output:
%   mergedTable_proc_sub - Filtered data table containing only sessions
%                          that match the specified data type criteria
%
% Example usage:
%   % Get sessions with both neural and pupil data
%   [dat_NP] = createDatSubset(mergedTable_proc, 'NP');
%
%   % Get sessions with behavioral data only
%   [dat_B] = createDatSubset(mergedTable_proc, 'B');
%
% Note: Session IDs must follow the naming convention "SessionName_CODE"
% where CODE contains letters B, N, and/or P indicating available data types.

    dat = mergedTable_proc;
    
    % Get unique session/unit names:
    
    unique_names = unique(dat.ses_ID);
    num_names = length(unique_names);
    
    % Filter by data type:
    
    mergedTable_proc_B = [];
    mergedTable_proc_BN = [];
    mergedTable_proc_BNP = [];
    mergedTable_proc_BP = [];
    mergedTable_proc_NP = [];
    mergedTable_proc_N = [];
    
    for ss = 1:num_names
    
        this_sesID = unique_names{ss};
        Lses = strcmp(dat.ses_ID, unique_names{ss});
        session_data = dat(Lses, 1:end);
    
        % Extract code for data types present
        % BNP (behavior, neural, pupil)
    
        sesID_trim_idx = find(this_sesID == '_', 1, 'last');
        sesID_data_code = this_sesID(sesID_trim_idx+1:end);
    
        if contains(sesID_data_code,'B'); % Behavioral data 
    
            mergedTable_proc_B = [mergedTable_proc_B;session_data];
    
            if contains(sesID_data_code,'BN'); % Behavioral and neural data
    
                mergedTable_proc_BN = [mergedTable_proc_BN;session_data];
    
                if contains(sesID_data_code,'BNP'); % Behavioral, neural, and pupil data
    
                    mergedTable_proc_BNP = [mergedTable_proc_BNP;session_data];
                end
            end
    
            if contains(sesID_data_code,'B') && contains(sesID_data_code,'P'); % Behavioral and pupil data
    
                mergedTable_proc_BP = [mergedTable_proc_BP;session_data];
            end 
        end

        if contains(sesID_data_code,'N'); % For neural and pupil comparison

                mergedTable_proc_N = [mergedTable_proc_N;session_data];

            if contains(sesID_data_code,'P'); % For neural and pupil comparison
    
                mergedTable_proc_NP = [mergedTable_proc_NP;session_data];
            end
        end
    end

    % Get desired output
    
    if strcmp(subset_code,'B') 
        mergedTable_proc_sub = mergedTable_proc_B;
    
    elseif strcmp(subset_code,'BN') 
        mergedTable_proc_sub = mergedTable_proc_BN;
    
    elseif strcmp(subset_code,'BNP') 
        mergedTable_proc_sub = mergedTable_proc_BNP;
    
    elseif strcmp(subset_code,'BP') 
        mergedTable_proc_sub = mergedTable_proc_BP;
    
    elseif strcmp(subset_code,'NP') 
        mergedTable_proc_sub = mergedTable_proc_NP;

    elseif strcmp(subset_code,'N') 
        mergedTable_proc_sub = mergedTable_proc_N;
    end

end




