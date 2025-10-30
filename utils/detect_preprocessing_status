function is_preprocessed = detect_preprocessing_status(DataMap)
% DETECT_PREPROCESSING_STATUS - Auto-detect if EEG data has been preprocessed
%
% Detection methods:
% 1. Check for common preprocessed data characteristics
% 2. Look for preprocessing report in data structure
% 3. Check channel names against common noisy channels

    is_preprocessed = false;
    
    if isempty(DataMap) || ~isvalid(DataMap)
        return;
    end
    
    channels = keys(DataMap);
    if isempty(channels)
        return;
    end
    
    % Method 1: Check if common noisy channels are already removed
    common_noisy_channels = {'FT9', 'FT10', 'TP9', 'TP10', 'F9', 'F10', 'P9', 'P10'};
    available_channels = upper(string(channels));
    noisy_present = any(ismember(upper(common_noisy_channels), available_channels));
    
    % If no noisy channels are present, likely preprocessed
    if ~noisy_present
        is_preprocessed = true;
        fprintf('   -> No common noisy channels detected - assuming preprocessed data\n');
        return;
    end
    
    % Method 2: Check data characteristics (simplified)
    % Sample a channel to check data properties
    sample_ch = channels{1};
    sample_data = DataMap(sample_ch);
    
    if isfield(sample_data, 'group1') && ~isempty(sample_data.group1)
        % Check first subject's data
        test_data = sample_data.group1{1};
        
        % Preprocessed data often has specific characteristics:
        % - Normalized amplitudes
        % - Specific frequency content
        % - Cleaner signal properties
        
        % Simple check: look for extreme outliers (raw data often has them)
        data_std = std(test_data);
        data_mean = mean(test_data);
        z_scores = abs(test_data - data_mean) / data_std;
        extreme_outliers = sum(z_scores > 5) / numel(test_data);
        
        if extreme_outliers < 0.01 % Less than 1% extreme outliers
            is_preprocessed = true;
            fprintf('   -> Data appears clean - assuming preprocessed\n');
        end
    end
    
    % Method 3: Check for preprocessing metadata
    % This would require the data to be saved with preprocessing report
    try
        % Try to load the original file and check for preprocessing info
        [filepath, name, ext] = fileparts(DataMap.source_file); % Assuming source file is stored
        possible_report_file = fullfile(filepath, 'preprocessing_report.mat');
        if exist(possible_report_file, 'file')
            is_preprocessed = true;
            fprintf('   -> Preprocessing report found - data is preprocessed\n');
        end
    catch
        % If we can't check, continue with other methods
    end
end
