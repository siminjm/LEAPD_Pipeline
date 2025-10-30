function additional_noisy = detect_additional_noisy_channels(DataMap)
% DETECT_ADDITIONAL_NOISY_CHANNELS - Automatically detect noisy channels in raw data

    additional_noisy = {};
    
    if isempty(DataMap) || ~isvalid(DataMap)
        return;
    end
    
    channels = keys(DataMap);
    fprintf('Scanning %d channels for noise...\n', numel(channels));
    
    noise_scores = containers.Map();
    
    for i = 1:numel(channels)
        ch = channels{i};
        chData = DataMap(ch);
        
        % Combine all subjects' data for this channel
        all_ch_data = [];
        if isfield(chData, 'group1') && ~isempty(chData.group1)
            for subj = 1:numel(chData.group1)
                all_ch_data = [all_ch_data; chData.group1{subj}(:)]; %#ok<AGROW>
            end
        end
        if isfield(chData, 'group2') && ~isempty(chData.group2)
            for subj = 1:numel(chData.group2)
                all_ch_data = [all_ch_data; chData.group2{subj}(:)]; %#ok<AGROW>
            end
        end
        
        if isempty(all_ch_data)
            continue;
        end
        
        % Calculate noise metrics
        % 1. Standard deviation (too high or too low)
        data_std = std(all_ch_data);
        
        % 2. Kurtosis (outliers)
        data_kurt = kurtosis(all_ch_data);
        
        % 3. Line noise ratio (60Hz power)
        Fs = 500; % Assuming standard sampling rate
        [pxx, f] = pwelch(all_ch_data, [], [], [], Fs);
        line_noise_band = (f >= 58 & f <= 62);
        broadband = (f >= 1 & f <= 45);
        line_noise_ratio = mean(pxx(line_noise_band)) / mean(pxx(broadband));
        
        % Combined noise score
        noise_score = data_kurt * line_noise_ratio * abs(log(data_std));
        noise_scores(ch) = noise_score;
    end
    
    % Identify outliers in noise scores
    if ~isempty(noise_scores)
        scores = cell2mat(values(noise_scores));
        score_names = keys(noise_scores);
        
        z_scores = abs(scores - mean(scores)) / std(scores);
        noisy_idx = find(z_scores > 3); % 3 standard deviations
        
        if ~isempty(noisy_idx)
            additional_noisy = score_names(noisy_idx);
            fprintf('Automatically detected noisy channels: %s\n', strjoin(additional_noisy, ', '));
        else
            fprintf('No additional noisy channels detected automatically\n');
        end
    end
end
