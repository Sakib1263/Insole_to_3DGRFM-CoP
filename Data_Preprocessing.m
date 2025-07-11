%% Prepare IMU2VGRF Data for Walk or Run
clear;
clc;
% Configurations
common_sampling_rate = 256;  % Resample each stream to this rate
overlapping_percentage = 0;  % in percent, vary between 50 to 90
augmentation_factor = (100-overlapping_percentage)/100;
final_segment_length = 1024;
% Global Constants
weight_array = [75.2, 65.1, 56.0, 80.8, 61.5, 81.0, 61.6, 62.4, 69.0];
ch_IMU_min = Inf;
ch_CoM_min = Inf;
ch_Force_FP_RL_min = Inf;
ch_Force_FP_LL_min = Inf;
ch_Insole_RL_min = Inf;
ch_Insole_LL_min = Inf;
ch_EMG_min = Inf;
trial_mode = "Run";
% Directory Information
rootdir = 'Data';
rootdir_info = dir(rootdir);
foldernames = {rootdir_info.name};
foldernames = foldernames(3:end);
% Global Arrays
data_IMU_segmented_all = zeros(10000,final_segment_length,50);
data_CoM_segmented_all = zeros(10000,final_segment_length,50);
data_Force_FP_RL_segmented_all = zeros(10000,final_segment_length,50);
data_Force_FP_LL_segmented_all = zeros(10000,final_segment_length,50);
data_Insole_RL_segmented_all = zeros(10000,final_segment_length,50);
data_Insole_LL_segmented_all = zeros(10000,final_segment_length,50);
data_EMG_segmented_all = zeros(10000,final_segment_length,50);
subinfo_all = zeros(10000,1);
speedinfo_all = zeros(10000,1);
% Loop through all Folders
% j=1:length(foldernames) % Assuming same folder names are present in both image and mask directories, as it is here
global_counter = 0;
for j=1:length(foldernames)
    foldername = cell2mat(foldernames(j));
    subjdir = fullfile(rootdir,foldername);
    subjdir_info = dir(subjdir);
    subjdir_foldernames = {subjdir_info.name};
    subjdir_foldernames = subjdir_foldernames(3:end);
    subinfo = str2double(extractAfter(foldername,4));
    weight_temp = weight_array(j);
    for jj=1:length(subjdir_foldernames)
        foldername = cell2mat(subjdir_foldernames(jj));
        if foldername == trial_mode
            trial_subjdir = fullfile(subjdir,foldername);
            trial_subjdir_info = dir(trial_subjdir);
            trial_subjdir_filenames = {trial_subjdir_info.name};
            trial_subjdir_filenames = trial_subjdir_filenames(3:end);
            for jjj=1:length(trial_subjdir_filenames)
                filename = cell2mat(trial_subjdir_filenames(jjj));
                fullfile_dir = fullfile(trial_subjdir,filename);
                disp(fullfile_dir);
                if trial_mode == "Walk"
                    speed_info = str2double(extractBefore(extractAfter(filename,12),'.mat'));
                elseif trial_mode == "Run"
                    speed_info = str2double(extractBefore(extractAfter(filename,11),'.mat'));
                end
                % Import Data
                data = load(fullfile_dir).Datastr;
                % Extract IMU information
                IMU_SF = data.IMU.IMUFrameRate; % IMU sampling frequency
                data_IMU = data.IMU.IMUData; % IMU data
                label_IMU = data.IMU.IMUDataLabel; % IMU labels
                data_CoM = data.IMU.CoM; % Center of Mass data
                timesamples_IMU = data.IMU.IMUmvnData.time';
                % Extract Force Plate Force information
                FP_GRF_SF = data.Force.FrameRate; % Force Plate (FP) sampling frequency
                data_Force_FP_RL = data.Force.RightForceData; % FP Right Leg Force Data
                data_Force_FP_LL = data.Force.LeftForceData; % FP Left Leg Force Data
                label_FP = data.Force.DataLabel; % FP labels
                % Extract Insole information
                Insole_GRF_SF = data.Insole.ResampleRate; % Insole sampling frequency
                data_Insole_RL = data.Insole.FilteredRight; % Insole Right Leg Force Data
                data_Insole_LL = data.Insole.FilteredLeft; % Insole Left Leg Force Data
                label_Insole = data.Insole.DataLabel; % Insole labels
                % Extract EMG information (Not required in current study)
                EMG_SF = data.EMG.FrameRate; % EMG sampling frequency
                data_EMG = data.EMG.Channels; % EMG Data
                label_EMG = data.EMG.DataLabel; % EMG labels
                % Get duration for each channel and determine the least
                data_IMU_duration = floor(length(data_IMU)/IMU_SF);
                data_Force_FP_RL_duration = floor(length(data_Force_FP_RL)/FP_GRF_SF);
                data_Force_FP_LL_duration = floor(length(data_Force_FP_LL)/FP_GRF_SF);
                data_Insole_RL_duration = floor(length(data_Insole_RL)/Insole_GRF_SF);
                data_Insole_LL_duration = floor(length(data_Insole_LL)/Insole_GRF_SF);
                data_EMG_duration = floor(length(data_EMG)/EMG_SF);
                minimum_duration = min([data_IMU_duration,data_Force_FP_RL_duration,data_Force_FP_LL_duration,data_Insole_RL_duration,data_Insole_LL_duration,data_EMG_duration]);
                % Cut extra points in the end before synchronization
                % data_IMU = data_IMU(end-minimum_duration*IMU_SF+1:end,:);
                % data_CoM = data_CoM(end-minimum_duration*IMU_SF+1:end,:);
                % data_Force_FP_RL = data_Force_FP_RL(end-minimum_duration*FP_GRF_SF+1:end,:);
                % data_Force_FP_LL = data_Force_FP_LL(end-minimum_duration*FP_GRF_SF+1:end,:);
                % data_Insole_RL = data_Insole_RL(end-minimum_duration*Insole_GRF_SF+1:end,:);
                % data_Insole_LL = data_Insole_LL(end-minimum_duration*Insole_GRF_SF+1:end,:);
                % data_EMG = data_EMG(end-minimum_duration*EMG_SF+1:end,:);
                data_IMU = data_IMU(1:minimum_duration*IMU_SF,:);
                data_CoM = data_CoM(1:minimum_duration*IMU_SF,:);
                data_Force_FP_RL = data_Force_FP_RL(1:minimum_duration*FP_GRF_SF,:);
                data_Force_FP_LL = data_Force_FP_LL(1:minimum_duration*FP_GRF_SF,:);
                data_Insole_RL = data_Insole_RL(1:minimum_duration*Insole_GRF_SF,:);
                data_Insole_LL = data_Insole_LL(1:minimum_duration*Insole_GRF_SF,:);
                data_EMG = data_EMG(1:minimum_duration*EMG_SF,:);
                % Resample data to a uniform rate
                data_IMU = resample(data_IMU,common_sampling_rate,IMU_SF);
                data_CoM = resample(data_CoM,common_sampling_rate,IMU_SF);
                data_Force_FP_RL = resample(data_Force_FP_RL,common_sampling_rate,FP_GRF_SF);
                data_Force_FP_LL = resample(data_Force_FP_LL,common_sampling_rate,FP_GRF_SF);
                data_Insole_RL = resample(data_Insole_RL,common_sampling_rate,Insole_GRF_SF);
                data_Insole_LL = resample(data_Insole_LL,common_sampling_rate,Insole_GRF_SF);
                data_EMG = resample(data_EMG,common_sampling_rate,EMG_SF);
                % Remove resampling affect in the extremity
                num_segments_waa = size(data_IMU,1)/common_sampling_rate;
                cropped_length = round((size(data_IMU,1)-(num_segments_waa*common_sampling_rate))/2);
                data_IMU = data_IMU(cropped_length+1:end-cropped_length,:);
                data_CoM = data_CoM(cropped_length+1:end-cropped_length,:);
                data_Force_FP_RL = data_Force_FP_RL(cropped_length+1:end-cropped_length,:);
                data_Force_FP_LL = data_Force_FP_LL(cropped_length+1:end-cropped_length,:);
                data_Insole_RL = data_Insole_RL(cropped_length+1:end-cropped_length,:);
                data_Insole_LL = data_Insole_LL(cropped_length+1:end-cropped_length,:);
                data_EMG = data_EMG(cropped_length+1:end-cropped_length,:);
                num_segments = round(size(data_IMU,1)/(common_sampling_rate*augmentation_factor));
                % Extract statistical information for zscore normalization
                mean_data_IMU = mean(data_IMU);
                mean_data_CoM = mean(data_CoM);
                mean_data_Force_FP_RL = mean(data_Force_FP_RL);
                mean_data_Force_FP_LL = mean(data_Force_FP_LL);
                mean_data_Insole_RL = mean(data_Insole_RL);
                mean_data_Insole_LL = mean(data_Insole_LL);
                mean_data_EMG = mean(data_EMG);
                sd_data_IMU = std(data_IMU);
                sd_data_CoM = std(data_CoM);
                sd_data_Force_FP_RL = std(data_Force_FP_RL);
                sd_data_Force_FP_LL = std(data_Force_FP_LL);
                sd_data_Insole_RL = std(data_Insole_RL);
                sd_data_Insole_LL = std(data_Insole_LL);
                sd_data_EMG = std(data_EMG);
                % Get number of channels for each data stream
                if size(data_IMU,2) < ch_IMU_min
                    ch_IMU_min = size(data_IMU,2);
                end
                if size(data_CoM,2) < ch_CoM_min
                    ch_CoM_min = size(data_CoM,2);
                end
                if size(data_Force_FP_RL,2) < ch_Force_FP_RL_min
                    ch_Force_FP_RL_min = size(data_Force_FP_RL,2);
                end
                if size(data_Force_FP_LL,2) < ch_Force_FP_LL_min
                    ch_Force_FP_LL_min = size(data_Force_FP_LL,2);
                end
                if size(data_Insole_RL,2) < ch_Insole_RL_min
                    ch_Insole_RL_min = size(data_Insole_RL,2);
                end
                if size(data_Insole_LL,2) < ch_Insole_LL_min
                    ch_Insole_LL_min = size(data_Insole_LL,2);
                end
                if size(data_EMG,2) < ch_EMG_min
                    ch_EMG_min = size(data_EMG,2);
                end
                % Create 1 second segments and resample
                data_IMU_segmented = zeros(1000,final_segment_length,size(data_IMU,2));
                data_CoM_segmented = zeros(1000,final_segment_length,size(data_CoM,2));
                data_Force_FP_RL_segmented = zeros(1000,final_segment_length,size(data_Force_FP_RL,2));
                data_Force_FP_LL_segmented = zeros(1000,final_segment_length,size(data_Force_FP_LL,2));
                data_Insole_RL_segmented = zeros(1000,final_segment_length,size(data_Insole_RL,2));
                data_Insole_LL_segmented = zeros(1000,final_segment_length,size(data_Insole_LL,2));
                data_EMG_segmented = zeros(1000,final_segment_length,size(data_EMG,2));
                counter = 0;
                for i=0:num_segments
                    if i*common_sampling_rate*augmentation_factor+final_segment_length > size(data_IMU,1)
                        continue
                    end
                    counter = counter + 1;
                    % uniform_array_resampling = linspace(1,common_sampling_rate*augmentation_factor,common_sampling_rate*augmentation_factor)';
                    data_IMU_temp = data_IMU((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_IMU_segmented(counter,:,:) = normalize(data_IMU_temp,'range');
                    % data_IMU_temp = resample(data_IMU_temp,final_sampling_rate,final_segment_length);
                    % data_IMU_temp = resample(data_IMU_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_IMU_temp = (data_IMU_temp - mean_data_IMU) ./ sd_data_IMU;
                    % data_IMU_segmented(counter,:,:) = normalize(data_IMU_temp(cut_amount+1:end-cut_amount,:),'range');
                    data_CoM_temp = data_CoM((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_CoM_segmented(counter,:,:) = normalize(data_CoM_temp,'range');
                    % data_CoM_temp = resample(data_CoM_temp,final_sampling_rate,final_segment_length);
                    % data_CoM_temp = resample(data_CoM_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_CoM_temp = (data_CoM_temp - mean_data_CoM) ./ sd_data_CoM;
                    % data_CoM_segmented(counter,:,:) = normalize(data_CoM_temp(cut_amount+1:end-cut_amount,:),'range');
                    data_Force_FP_RL_temp = data_Force_FP_RL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_Force_FP_RL_temp(:,1) = normalize(filtButter(data_Force_FP_RL_temp(:,1),final_segment_length,2,[0.01 20],'bandpass'),'range');
                    data_Force_FP_RL_temp(:,4) = normalize(filtButter(data_Force_FP_RL_temp(:,4),final_segment_length,2,[0.01 20],'bandpass'),'range');
                    data_Force_FP_RL_temp(:,5) = normalize(filtButter(data_Force_FP_RL_temp(:,5),final_segment_length,4,[0.01 20],'bandpass'),'range');
                    data_Force_FP_RL_temp(:,6) = normalize(filtButter(data_Force_FP_RL_temp(:,6),final_segment_length,2,[0.01 20],'bandpass'),'range');
                    data_Force_FP_RL_segmented(counter,:,:) = normalize(data_Force_FP_RL_temp,'range');
                    % data_Force_FP_RL_temp = resample(data_Force_FP_RL_temp,final_sampling_rate,final_segment_length);
                    % data_Force_FP_RL_temp = resample(data_Force_FP_RL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_Force_FP_RL_temp = data_Force_FP_RL_temp / weight_temp;
                    % data_Force_FP_RL_temp = (data_Force_FP_RL_temp - mean_data_Force_FP_RL) ./ sd_data_Force_FP_RL;
                    % data_Force_FP_RL_segmented(counter,:,:) = normalize(data_Force_FP_RL_temp(cut_amount+1:end-cut_amount,:),'range');
                    data_Force_FP_LL_temp = data_Force_FP_LL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_Force_FP_LL_temp(:,1) = normalize(filtButter(data_Force_FP_LL_temp(:,1),final_segment_length,2,[0.01 20],'bandpass'),'range');
                    data_Force_FP_LL_temp(:,4) = normalize(filtButter(data_Force_FP_LL_temp(:,4),final_segment_length,2,[0.01 20],'bandpass'),'range');
                    data_Force_FP_LL_temp(:,5) = normalize(filtButter(data_Force_FP_LL_temp(:,5),final_segment_length,4,[0.01 20],'bandpass'),'range');
                    data_Force_FP_LL_temp(:,6) = normalize(filtButter(data_Force_FP_LL_temp(:,6),final_segment_length,2,[0.01 20],'bandpass'),'range');
                    data_Force_FP_LL_segmented(counter,:,:) = normalize(data_Force_FP_LL_temp,'range');
                    % data_Force_FP_LL_temp = resample(data_Force_FP_LL_temp,final_sampling_rate,final_segment_length);
                    % data_Force_FP_LL_temp = resample(data_Force_FP_LL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_Force_FP_LL_temp = data_Force_FP_LL_temp / weight_temp;
                    % data_Force_FP_LL_temp = (data_Force_FP_LL_temp - mean_data_Force_FP_LL) ./ sd_data_Force_FP_LL;
                    % data_Force_FP_LL_segmented(counter,:,:) = normalize(data_Force_FP_LL_temp(cut_amount+1:end-cut_amount,:),'range');
                    data_Insole_RL_temp = data_Insole_RL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_Insole_RL_segmented(counter,:,:) = normalize(data_Insole_RL_temp,'range');
                    % data_Insole_RL_temp = resample(data_Insole_RL_temp,final_sampling_rate,final_segment_length);
                    % data_Insole_RL_temp = resample(data_Insole_RL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_Insole_RL_temp = data_Insole_RL_temp / weight_temp;
                    % data_Insole_RL_temp = (data_Insole_RL_temp - mean_data_Insole_RL) ./ sd_data_Insole_RL;
                    % data_Insole_RL_segmented(counter,:,:) = normalize(data_Insole_RL_temp(cut_amount+1:end-cut_amount,:),'range');
                    data_Insole_LL_temp = data_Insole_LL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_Insole_LL_segmented(counter,:,:) = normalize(data_Insole_LL_temp,'range');
                    % data_Insole_LL_temp = resample(data_Insole_LL_temp,final_sampling_rate,final_segment_length);
                    % data_Insole_LL_temp = resample(data_Insole_LL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_Insole_LL_temp = data_Insole_LL_temp / weight_temp;
                    % data_Insole_LL_temp = (data_Insole_LL_temp - mean_data_Insole_LL) ./ sd_data_Insole_LL;
                    % data_Insole_LL_segmented(counter,:,:) = normalize(data_Insole_LL_temp(cut_amount+1:end-cut_amount,:),'range');
                    data_EMG_temp = data_EMG((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                    data_EMG_segmented(counter,:,:) = normalize(data_EMG_temp,'range');
                    % data_EMG_temp = resample(data_EMG_temp,final_sampling_rate,final_segment_length);
                    % data_EMG_temp = resample(data_EMG_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                    % data_EMG_temp = (data_EMG_temp - mean_data_EMG) ./ sd_data_EMG;
                    % data_EMG_segmented(counter,:,:) = normalize(data_EMG_temp(cut_amount+1:end-cut_amount,:),'range');
                end
                data_IMU_segmented = data_IMU_segmented(1:counter,:,:);
                data_CoM_segmented = data_CoM_segmented(1:counter,:,:);
                data_Force_FP_RL_segmented = data_Force_FP_RL_segmented(1:counter,:,:);
                data_Force_FP_LL_segmented = data_Force_FP_LL_segmented(1:counter,:,:);
                data_Insole_RL_segmented = data_Insole_RL_segmented(1:counter,:,:);
                data_Insole_LL_segmented = data_Insole_LL_segmented(1:counter,:,:);
                data_EMG_segmented = data_EMG_segmented(1:counter,:,:);
                subinfo_current = zeros(size(data_IMU_segmented,1),1) + subinfo;
                speedinfo_current = zeros(size(data_IMU_segmented,1),1) + speed_info;
                %
                global_counter = global_counter + counter;
                data_IMU_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_IMU_segmented,3)) = data_IMU_segmented;
                data_CoM_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_CoM_segmented,3)) = data_CoM_segmented;
                data_Force_FP_RL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Force_FP_RL_segmented,3)) = data_Force_FP_RL_segmented;
                data_Force_FP_LL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Force_FP_LL_segmented,3)) = data_Force_FP_LL_segmented;
                data_Insole_RL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Insole_RL_segmented,3)) = data_Insole_RL_segmented;
                data_Insole_LL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Insole_LL_segmented,3)) = data_Insole_LL_segmented;
                data_EMG_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_EMG_segmented,3)) = data_EMG_segmented;
                subinfo_all(global_counter-counter+1:global_counter,1) = subinfo_current;
                speedinfo_all(global_counter-counter+1:global_counter,1) = speedinfo_current;
            end
        end
    end
end
data_IMU = data_IMU_segmented_all(1:global_counter,:,1:ch_IMU_min);
data_CoM = data_CoM_segmented_all(1:global_counter,:,1:ch_CoM_min);
data_Force_FP_RL = data_Force_FP_RL_segmented_all(1:global_counter,:,1:ch_Force_FP_RL_min);
data_Force_FP_LL = data_Force_FP_LL_segmented_all(1:global_counter,:,1:ch_Force_FP_LL_min);
data_Insole_RL = data_Insole_RL_segmented_all(1:global_counter,:,1:ch_Insole_RL_min);
data_Insole_LL = data_Insole_LL_segmented_all(1:global_counter,:,1:ch_Insole_LL_min);
data_EMG = data_EMG_segmented_all(1:global_counter,:,1:ch_EMG_min);
subinfo = subinfo_all(1:global_counter,1);
speedinfo = speedinfo_all(1:global_counter,1);
% Save
save_path_mat = sprintf('Data_Preprocessed_%s.mat',trial_mode);
save(save_path_mat,'data_IMU','data_CoM','data_Force_FP_RL','data_Force_FP_LL','data_Insole_RL',...
    'data_Insole_LL','data_EMG','subinfo','speedinfo','label_IMU','label_FP','label_Insole','label_EMG','-v7.3');
%% Prepare IMU2VGRF Data for Locomotion - [Jump, Land, Lunge or Squat]
clear;
clc;
% Configurations
common_sampling_rate = 256;  % Resample each stream to this rate
overlapping_percentage = 0;  % in percent, vary between 50 to 90
augmentation_factor = (100-overlapping_percentage)/100;
final_segment_length = 1024;
% Global Constants
weight_array = [75.2, 65.1, 56.0, 80.8, 61.5, 81.0, 61.6, 62.4, 69.0];
ch_IMU_min = Inf;
ch_CoM_min = Inf;
ch_Force_FP_RL_min = Inf;
ch_Force_FP_LL_min = Inf;
ch_Insole_RL_min = Inf;
ch_Insole_LL_min = Inf;
ch_EMG_min = Inf;
trial_mode = "squat";
% Directory Information
rootdir = 'Data';
rootdir_info = dir(rootdir);
foldernames = {rootdir_info.name};
foldernames = foldernames(3:end);
% Global Arrays
data_IMU_segmented_all = zeros(10000,final_segment_length,50);
data_CoM_segmented_all = zeros(10000,final_segment_length,50);
data_Force_FP_RL_segmented_all = zeros(10000,final_segment_length,50);
data_Force_FP_LL_segmented_all = zeros(10000,final_segment_length,50);
data_Insole_RL_segmented_all = zeros(10000,final_segment_length,50);
data_Insole_LL_segmented_all = zeros(10000,final_segment_length,50);
data_EMG_segmented_all = zeros(10000,final_segment_length,50);
subinfo_all = zeros(10000,1);
% Loop through all Folders
% j=1:length(foldernames) % Assuming same folder names are present in both image and mask directories, as it is here
global_counter = 0;
for j=1:length(foldernames)
    foldername = cell2mat(foldernames(j));
    subjdir = fullfile(rootdir,foldername);
    subjdir_info = dir(subjdir);
    subjdir_foldernames = {subjdir_info.name};
    subjdir_foldernames = subjdir_foldernames(3:end);
    subinfo = str2double(extractAfter(foldername,4));
    weight_temp = weight_array(j);
    for jj=1:length(subjdir_foldernames)
        foldername = cell2mat(subjdir_foldernames(jj));
        if foldername == "Locomotion"
            trial_subjdir = fullfile(subjdir,foldername);
            trial_subjdir_info = dir(trial_subjdir);
            trial_subjdir_filenames = {trial_subjdir_info.name};
            trial_subjdir_filenames = trial_subjdir_filenames(3:end);
            for jjj=1:length(trial_subjdir_filenames)
                filename = cell2mat(trial_subjdir_filenames(jjj));
                fullfile_dir = fullfile(trial_subjdir,filename);
                locomotion_type = string(extractBefore(extractAfter(filename,7),'.mat'));
                if strcmp(locomotion_type,trial_mode) == 1
                    disp(fullfile_dir);
                    % Import Data
                    data = load(fullfile_dir).Datastr;
                    % Extract IMU information
                    IMU_SF = data.IMU.IMUFrameRate; % IMU sampling frequency
                    data_IMU = data.IMU.IMUData; % IMU data
                    label_IMU = data.IMU.IMUDataLabel; % IMU labels
                    data_CoM = data.IMU.CoM; % Center of Mass data
                    timesamples_IMU = data.IMU.IMUmvnData.time';
                    % Extract Force Plate Force information
                    FP_GRF_SF = data.Force.FrameRate; % Force Plate (FP) sampling frequency
                    data_Force_FP_RL = data.Force.RightForceData; % FP Right Leg Force Data
                    data_Force_FP_LL = data.Force.LeftForceData; % FP Left Leg Force Data
                    label_FP = data.Force.DataLabel; % FP labels
                    % Extract Insole information
                    Insole_GRF_SF = data.Insole.ResampleRate; % Insole sampling frequency
                    data_Insole_RL = data.Insole.FilteredRight; % Insole Right Leg Force Data
                    data_Insole_LL = data.Insole.FilteredLeft; % Insole Left Leg Force Data
                    label_Insole = data.Insole.DataLabel; % Insole labels
                    % Extract EMG information (Not required in current study)
                    EMG_SF = data.EMG.FrameRate; % EMG sampling frequency
                    data_EMG = data.EMG.Channels; % EMG Data
                    label_EMG = data.EMG.DataLabel; % EMG labels
                    % Get duration for each channel and determine the least
                    data_IMU_duration = floor(length(data_IMU)/IMU_SF);
                    data_Force_FP_RL_duration = floor(length(data_Force_FP_RL)/FP_GRF_SF);
                    data_Force_FP_LL_duration = floor(length(data_Force_FP_LL)/FP_GRF_SF);
                    data_Insole_RL_duration = floor(length(data_Insole_RL)/Insole_GRF_SF);
                    data_Insole_LL_duration = floor(length(data_Insole_LL)/Insole_GRF_SF);
                    data_EMG_duration = floor(length(data_EMG)/EMG_SF);
                    minimum_duration = min([data_IMU_duration,data_Force_FP_RL_duration,data_Force_FP_LL_duration,data_Insole_RL_duration,data_Insole_LL_duration,data_EMG_duration]);
                    % Cut extra points in the end before synchronization
                    data_IMU = data_IMU(end-minimum_duration*IMU_SF+1:end,:);
                    data_CoM = data_CoM(end-minimum_duration*IMU_SF+1:end,:);
                    data_Force_FP_RL = data_Force_FP_RL(end-minimum_duration*FP_GRF_SF+1:end,:);
                    data_Force_FP_LL = data_Force_FP_LL(end-minimum_duration*FP_GRF_SF+1:end,:);
                    data_Insole_RL = data_Insole_RL(end-minimum_duration*Insole_GRF_SF+1:end,:);
                    data_Insole_LL = data_Insole_LL(end-minimum_duration*Insole_GRF_SF+1:end,:);
                    data_EMG = data_EMG(end-minimum_duration*EMG_SF+1:end,:);
                    % Resample data to a uniform rate
                    data_IMU = resample(data_IMU,common_sampling_rate,IMU_SF);
                    data_CoM = resample(data_CoM,common_sampling_rate,IMU_SF);
                    data_Force_FP_RL = resample(data_Force_FP_RL,common_sampling_rate,FP_GRF_SF);
                    data_Force_FP_LL = resample(data_Force_FP_LL,common_sampling_rate,FP_GRF_SF);
                    data_Insole_RL = resample(data_Insole_RL,common_sampling_rate,Insole_GRF_SF);
                    data_Insole_LL = resample(data_Insole_LL,common_sampling_rate,Insole_GRF_SF);
                    data_EMG = resample(data_EMG,common_sampling_rate,EMG_SF);
                    % Remove resampling affect in the extremity
                    num_segments_waa = size(data_IMU,1)/common_sampling_rate;
                    cropped_length = round((size(data_IMU,1)-(num_segments_waa*common_sampling_rate))/2);
                    data_IMU = data_IMU(cropped_length+1:end-cropped_length,:);
                    data_CoM = data_CoM(cropped_length+1:end-cropped_length,:);
                    data_Force_FP_RL = data_Force_FP_RL(cropped_length+1:end-cropped_length,:);
                    data_Force_FP_LL = data_Force_FP_LL(cropped_length+1:end-cropped_length,:);
                    data_Insole_RL = data_Insole_RL(cropped_length+1:end-cropped_length,:);
                    data_Insole_LL = data_Insole_LL(cropped_length+1:end-cropped_length,:);
                    data_EMG = data_EMG(cropped_length+1:end-cropped_length,:);
                    num_segments = round(size(data_IMU,1)/(common_sampling_rate*augmentation_factor));
                    % Extract statistical information for zscore normalization
                    mean_data_IMU = mean(data_IMU);
                    mean_data_CoM = mean(data_CoM);
                    mean_data_Force_FP_RL = mean(data_Force_FP_RL);
                    mean_data_Force_FP_LL = mean(data_Force_FP_LL);
                    mean_data_Insole_RL = mean(data_Insole_RL);
                    mean_data_Insole_LL = mean(data_Insole_LL);
                    mean_data_EMG = mean(data_EMG);
                    sd_data_IMU = std(data_IMU);
                    sd_data_CoM = std(data_CoM);
                    sd_data_Force_FP_RL = std(data_Force_FP_RL);
                    sd_data_Force_FP_LL = std(data_Force_FP_LL);
                    sd_data_Insole_RL = std(data_Insole_RL);
                    sd_data_Insole_LL = std(data_Insole_LL);
                    sd_data_EMG = std(data_EMG);
                    % Get number of channels for each data stream
                    if size(data_IMU,2) < ch_IMU_min
                        ch_IMU_min = size(data_IMU,2);
                    end
                    if size(data_CoM,2) < ch_CoM_min
                        ch_CoM_min = size(data_CoM,2);
                    end
                    if size(data_Force_FP_RL,2) < ch_Force_FP_RL_min
                        ch_Force_FP_RL_min = size(data_Force_FP_RL,2);
                    end
                    if size(data_Force_FP_LL,2) < ch_Force_FP_LL_min
                        ch_Force_FP_LL_min = size(data_Force_FP_LL,2);
                    end
                    if size(data_Insole_RL,2) < ch_Insole_RL_min
                        ch_Insole_RL_min = size(data_Insole_RL,2);
                    end
                    if size(data_Insole_LL,2) < ch_Insole_LL_min
                        ch_Insole_LL_min = size(data_Insole_LL,2);
                    end
                    if size(data_EMG,2) < ch_EMG_min
                        ch_EMG_min = size(data_EMG,2);
                    end
                    % Create 1 second segments and resample
                    data_IMU_segmented = zeros(1000,final_segment_length,size(data_IMU,2));
                    data_CoM_segmented = zeros(1000,final_segment_length,size(data_CoM,2));
                    data_Force_FP_RL_segmented = zeros(1000,final_segment_length,size(data_Force_FP_RL,2));
                    data_Force_FP_LL_segmented = zeros(1000,final_segment_length,size(data_Force_FP_LL,2));
                    data_Insole_RL_segmented = zeros(1000,final_segment_length,size(data_Insole_RL,2));
                    data_Insole_LL_segmented = zeros(1000,final_segment_length,size(data_Insole_LL,2));
                    data_EMG_segmented = zeros(1000,final_segment_length,size(data_EMG,2));
                    counter = 0;
                    for i=0:num_segments
                        if i*common_sampling_rate*augmentation_factor+final_segment_length > size(data_IMU,1)
                            continue
                        end
                        counter = counter + 1;
                        % uniform_array_resampling = linspace(1,common_sampling_rate*augmentation_factor,common_sampling_rate*augmentation_factor)';
                        data_IMU_temp = data_IMU((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_IMU_segmented(counter,:,:) = normalize(data_IMU_temp,'range');
                        % data_IMU_temp = resample(data_IMU_temp,final_sampling_rate,final_segment_length);
                        % data_IMU_temp = resample(data_IMU_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_IMU_temp = (data_IMU_temp - mean_data_IMU) ./ sd_data_IMU;
                        % data_IMU_segmented(counter,:,:) = normalize(data_IMU_temp(cut_amount+1:end-cut_amount,:),'range');
                        data_CoM_temp = data_CoM((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_CoM_segmented(counter,:,:) = normalize(data_CoM_temp,'range');
                        % data_CoM_temp = resample(data_CoM_temp,final_sampling_rate,final_segment_length);
                        % data_CoM_temp = resample(data_CoM_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_CoM_temp = (data_CoM_temp - mean_data_CoM) ./ sd_data_CoM;
                        % data_CoM_segmented(counter,:,:) = normalize(data_CoM_temp(cut_amount+1:end-cut_amount,:),'range');
                        data_Force_FP_RL_temp = data_Force_FP_RL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_Force_FP_RL_segmented(counter,:,:) = normalize(data_Force_FP_RL_temp,'range');
                        % data_Force_FP_RL_temp = resample(data_Force_FP_RL_temp,final_sampling_rate,final_segment_length);
                        % data_Force_FP_RL_temp = resample(data_Force_FP_RL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_Force_FP_RL_temp = data_Force_FP_RL_temp / weight_temp;
                        % data_Force_FP_RL_temp = (data_Force_FP_RL_temp - mean_data_Force_FP_RL) ./ sd_data_Force_FP_RL;
                        % data_Force_FP_RL_segmented(counter,:,:) = normalize(data_Force_FP_RL_temp(cut_amount+1:end-cut_amount,:),'range');
                        data_Force_FP_LL_temp = data_Force_FP_LL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_Force_FP_LL_segmented(counter,:,:) = normalize(data_Force_FP_LL_temp,'range');
                        % data_Force_FP_LL_temp = resample(data_Force_FP_LL_temp,final_sampling_rate,final_segment_length);
                        % data_Force_FP_LL_temp = resample(data_Force_FP_LL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_Force_FP_LL_temp = data_Force_FP_LL_temp / weight_temp;
                        % data_Force_FP_LL_temp = (data_Force_FP_LL_temp - mean_data_Force_FP_LL) ./ sd_data_Force_FP_LL;
                        % data_Force_FP_LL_segmented(counter,:,:) = normalize(data_Force_FP_LL_temp(cut_amount+1:end-cut_amount,:),'range');
                        data_Insole_RL_temp = data_Insole_RL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_Insole_RL_segmented(counter,:,:) = normalize(data_Insole_RL_temp,'range');
                        % data_Insole_RL_temp = resample(data_Insole_RL_temp,final_sampling_rate,final_segment_length);
                        % data_Insole_RL_temp = resample(data_Insole_RL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_Insole_RL_temp = data_Insole_RL_temp / weight_temp;
                        % data_Insole_RL_temp = (data_Insole_RL_temp - mean_data_Insole_RL) ./ sd_data_Insole_RL;
                        % data_Insole_RL_segmented(counter,:,:) = normalize(data_Insole_RL_temp(cut_amount+1:end-cut_amount,:),'range');
                        data_Insole_LL_temp = data_Insole_LL((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_Insole_LL_segmented(counter,:,:) = normalize(data_Insole_LL_temp,'range');
                        % data_Insole_LL_temp = resample(data_Insole_LL_temp,final_sampling_rate,final_segment_length);
                        % data_Insole_LL_temp = resample(data_Insole_LL_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_Insole_LL_temp = data_Insole_LL_temp / weight_temp;
                        % data_Insole_LL_temp = (data_Insole_LL_temp - mean_data_Insole_LL) ./ sd_data_Insole_LL;
                        % data_Insole_LL_segmented(counter,:,:) = normalize(data_Insole_LL_temp(cut_amount+1:end-cut_amount,:),'range');
                        data_EMG_temp = data_EMG((i*common_sampling_rate*augmentation_factor)+1:i*common_sampling_rate*augmentation_factor+final_segment_length,:);
                        data_EMG_segmented(counter,:,:) = normalize(data_EMG_temp,'range');
                        % data_EMG_temp = resample(data_EMG_temp,final_sampling_rate,final_segment_length);
                        % data_EMG_temp = resample(data_EMG_temp,uniform_array_resampling,round(final_sampling_rate/(common_sampling_rate*augmentation_factor)),'pchip');
                        % data_EMG_temp = (data_EMG_temp - mean_data_EMG) ./ sd_data_EMG;
                        % data_EMG_segmented(counter,:,:) = normalize(data_EMG_temp(cut_amount+1:end-cut_amount,:),'range');
                    end
                    data_IMU_segmented = data_IMU_segmented(1:counter,:,:);
                    data_CoM_segmented = data_CoM_segmented(1:counter,:,:);
                    data_Force_FP_RL_segmented = data_Force_FP_RL_segmented(1:counter,:,:);
                    data_Force_FP_LL_segmented = data_Force_FP_LL_segmented(1:counter,:,:);
                    data_Insole_RL_segmented = data_Insole_RL_segmented(1:counter,:,:);
                    data_Insole_LL_segmented = data_Insole_LL_segmented(1:counter,:,:);
                    data_EMG_segmented = data_EMG_segmented(1:counter,:,:);
                    subinfo_current = zeros(size(data_IMU_segmented,1),1) + subinfo;
                    %
                    global_counter = global_counter + counter;
                    data_IMU_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_IMU_segmented,3)) = data_IMU_segmented;
                    data_CoM_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_CoM_segmented,3)) = data_CoM_segmented;
                    data_Force_FP_RL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Force_FP_RL_segmented,3)) = data_Force_FP_RL_segmented;
                    data_Force_FP_LL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Force_FP_LL_segmented,3)) = data_Force_FP_LL_segmented;
                    data_Insole_RL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Insole_RL_segmented,3)) = data_Insole_RL_segmented;
                    data_Insole_LL_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_Insole_LL_segmented,3)) = data_Insole_LL_segmented;
                    data_EMG_segmented_all(global_counter-counter+1:global_counter,:,1:size(data_EMG_segmented,3)) = data_EMG_segmented;
                    subinfo_all(global_counter-counter+1:global_counter,1) = subinfo_current;
                end
            end
        end
    end
end
data_IMU = data_IMU_segmented_all(1:global_counter,:,1:ch_IMU_min);
data_CoM = data_CoM_segmented_all(1:global_counter,:,1:ch_CoM_min);
data_Force_FP_RL = data_Force_FP_RL_segmented_all(1:global_counter,:,1:ch_Force_FP_RL_min);
data_Force_FP_LL = data_Force_FP_LL_segmented_all(1:global_counter,:,1:ch_Force_FP_LL_min);
data_Insole_RL = data_Insole_RL_segmented_all(1:global_counter,:,1:ch_Insole_RL_min);
data_Insole_LL = data_Insole_LL_segmented_all(1:global_counter,:,1:ch_Insole_LL_min);
data_EMG = data_EMG_segmented_all(1:global_counter,:,1:ch_EMG_min);
subinfo = subinfo_all(1:global_counter,1);
% Save
save_path_mat = sprintf('Data_Preprocessed_%s.mat',trial_mode);
save(save_path_mat,'data_IMU','data_CoM','data_Force_FP_RL','data_Force_FP_LL','data_Insole_RL',...
    'data_Insole_LL','data_EMG','subinfo','label_IMU','label_FP','label_Insole','label_EMG','-v7.3');
%% Load Preprocessed Data
clear;
clc;
data_structure = load('Data_Preprocessed_Run.mat');
data_IMU = data_structure.data_IMU;
data_CoM = data_structure.data_CoM;
data_Force_FP_RL = data_structure.data_Force_FP_RL;
data_Force_FP_LL = data_structure.data_Force_FP_LL;
data_Insole_RL = data_structure.data_Insole_RL;
data_Insole_LL = data_structure.data_Insole_LL;
data_EMG = data_structure.data_EMG;
label_IMU = data_structure.label_IMU;
label_FP = data_structure.label_FP;
label_Insole = data_structure.label_Insole;
label_EMG = data_structure.label_EMG;
subinfo = data_structure.subinfo;
speedinfo = data_structure.speedinfo;
%% Plot
sample_num = randi(size(data_IMU,1),1);
% sample_num = sample_num + 10;
segment_length = size(data_IMU,2);
data_IMU_temp = squeeze(data_IMU(sample_num,:,:));
data_Force_FP_RL_temp = squeeze(data_Force_FP_RL(sample_num,:,:));
data_Force_FP_LL_temp = squeeze(data_Force_FP_LL(sample_num,:,:));
subinfo_temp = subinfo(sample_num,1);
speedinfo_temp = speedinfo(sample_num,1);
%
figure;
tiledlayout(4,2,'Padding','compact','TileSpacing','compact');
nexttile;
plot(data_IMU_temp(:,1),'LineWidth',2);
xlim([0 length(data_IMU_temp(:,1))])
title(sprintf('IMU Channel: %s', string(label_IMU(1,1))),'FontSize',12);
nexttile;
x = data_Force_FP_RL_temp(:,1);
y = normalize(filtButter(x,segment_length,2,[0.01 20],'bandpass'),'range');
hold on
plot(x,'LineWidth',2);
plot(y,'LineWidth',2);
hold off
xlim([0 length(data_Force_FP_RL_temp(:,1))])
title(sprintf('Force Plate Channel: %s', string(label_FP(1,1))),'FontSize',12);
nexttile;
plot(data_Force_FP_RL_temp(:,2),'LineWidth',2);
xlim([0 length(data_Force_FP_RL_temp(:,2))])
title(sprintf('Force Plate Channel: %s', string(label_FP(1,2))),'FontSize',12);
nexttile;
plot(data_Force_FP_RL_temp(:,3),'LineWidth',2);
xlim([0 length(data_Force_FP_RL_temp(:,3))])
title(sprintf('Force Plate Channel: %s', string(label_FP(1,3))),'FontSize',12);
nexttile;
plot(data_Force_FP_RL_temp(:,4),'LineWidth',2);
xlim([0 length(data_Force_FP_RL_temp(:,4))])
title(sprintf('Force Plate Channel: %s', string(label_FP(1,4))),'FontSize',12);
nexttile;
x = data_Force_FP_RL_temp(:,4);
y = normalize(filtButter(x,segment_length,2,[0.01 20],'bandpass'),'range');
hold on
plot(x,'LineWidth',2);
plot(y,'LineWidth',2);
hold off
xlim([0 length(data_Force_FP_RL_temp(:,5))])
title(sprintf('Force Plate Channel: %s', string(label_FP(1,5))),'FontSize',12);
nexttile;
plot(data_Force_FP_RL_temp(:,6),'LineWidth',2);
xlim([0 length(data_Force_FP_RL_temp(:,6))])
title(sprintf('Force Plate Channel: %s', string(label_FP(1,6))),'FontSize',12);
sgtitle(sprintf('Subject: %d - Speed: %d', subinfo_temp, speedinfo_temp),'FontSize',16);
%%

