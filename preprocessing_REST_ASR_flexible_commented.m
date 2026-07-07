%% EEG preprocessing with optional REST, ASR, ICA, MARA, ICLabel, and FASTER steps
%
% PURPOSE
% This script preprocesses resting-state EEG recordings stored as EDF file and includes the following
% possible processing steps:
%   1. Load EDF files with EEGLAB/BIOSIG.
%   2. Optionally select a user-defined subset of EEG channels.
%   3. Rename channels by removing a user-defined prefix, for example "EEG ".
%   4. Add standard channel locations.
%   5. Detect and interpolate globally bad channels using kurtosis.
%   6. Remove line noise using CleanLine.
%   7. Optionally re-reference the signal to a point in infinity (REST).
%   8. Apply high-pass filtering.
%   9. Optionally import event files.
%   10. Optionally remove the worst continuous artifacts using the ASR window
%       rejection procedure from clean_rawdata.
%   11. Optionally run ICA and reject components identified by MARA and/or
%       ICLabel.
%   12. Optionally create regular 1-second epochs, detect bad channels within
%       each epoch using FASTER functions, interpolate them, and convert the
%       data back to continuous format.
%   13. Save the final preprocessed EEGLAB .set files and a preprocessing log.
%
% INPUTS
%   - Raw EEG files in EDF format, selected using cfg.rawFilePattern.
%   - Optional event text files with columns: latency, duration, type.
%   - A channel-location file compatible with EEGLAB pop_chanedit.
%
% OUTPUTS
%   - Final preprocessed EEGLAB .set files saved to cfg.outputFinalDir.
%   - Optional intermediate filtered .set files saved to cfg.outputFilteredDir.
%   - A MATLAB .mat file containing:
%       preprocessingLog          main preprocessing information per file
%       badFiles                  files that failed preprocessing checks
%       filesRequiringEventCheck  files with suspicious event timing
%       cfg                       parameter settings used for this run
%
% REQUIRED TOOLBOXES / PLUGINS
%   - MATLAB.
%   - EEGLAB.
%   - BIOSIG plugin for EDF import: pop_biosig.
%   - CleanLine plugin: pop_cleanline.
%   - clean_rawdata plugin if cfg.runASRWindowRejection = true.
%   - REST functions if cfg.runRESTReference = true:
%       dong_getleadfield, rest_refer.
%   - FASTER functions if cfg.runEpochChannelInterpolation = true:
%       single_epoch_channel_properties, min_z_changed, h_epoch_interp_spl.
%   - MARA plugin if cfg.runMARA = true.
%   - ICLabel plugin if cfg.runICLabel = true.
%   - DIPFIT or another EEGLAB-compatible channel-location file.
%
% NOTES
%   - To process all channels present in each file, set cfg.rawChannelLabelsToKeep
%     to an empty cell array: cfg.rawChannelLabelsToKeep = {};
%   - To enforce a common channel set and order across files, define
%     cfg.rawChannelLabelsToKeep with the raw channel labels as they appear in
%     the EDF files.
%   - If REST is enabled, the leadfield is calculated for all currently selected
%     channels. Make sure the channel labels and locations are valid.

clearvars;
clc;

%% ======================= USER PARAMETERS ===============================

% ----------------------- Paths ------------------------------------------
cfg.eeglabPath = 'C:\Program Files\eeglab2021.1';
cfg.extraPluginPaths = { ...
    'C:\Program Files\eeglab2021.1\plugins\FASTER1.2.4\FASTER1.2.4' ...
    };

cfg.rawDataDir = 'E:\DANE IPiN\exported_edf';
cfg.outputFinalDir = 'E:\dane backup\IPiN_preprocessed_REST_ASR';
cfg.outputFilteredDir = 'E:\dane backup\IPiN_preprocessed_REST_ASR\only_REST_filter';
cfg.eventDir = 'E:\DANE IPiN\eventy\events_txt';
cfg.channelLocationFile = 'C:\Program Files\eeglab2021.1\plugins\dipfit\standard_BEM\elec\standard_1005.elc';

% ----------------------- File naming ------------------------------------
% Pattern used to find raw EDF files.
cfg.rawFilePattern = '*base.edf';

% Suffix removed from the EDF file name to create the base name used for
% output files and event-file matching. For example, patient_base.edf becomes
% patient.set and patient.txt when cfg.rawFileSuffixToRemove = '_base.edf'.
cfg.rawFileSuffixToRemove = '_base.edf';

% ----------------------- Channel handling -------------------------------
% Raw channel labels to keep, in the desired order. Leave empty to keep all
% channels available in each file.
cfg.rawChannelLabelsToKeep = { ...
    'EEG Fp1', 'EEG Fp2', 'EEG Fz', 'EEG F3', 'EEG F7', ...
    'EEG F4', 'EEG F8', 'EEG T3', 'EEG T4', 'EEG T5', ...
    'EEG T6', 'EEG C3', 'EEG Cz', 'EEG C4', 'EEG P3', ...
    'EEG Pz', 'EEG P4', 'EEG O1', 'EEG O2' ...
    };

% Prefix removed from channel labels after channel selection. For example,
% "EEG Fp1" becomes "Fp1" when this is set to 'EEG '. Set to '' if labels
% should not be changed.
cfg.channelLabelPrefixToRemove = 'EEG ';

cfg.addChannelLocations = true;

% ----------------------- Main optional processing steps ------------------
cfg.runInitialBadChannelInterpolation = true;
cfg.runCleanLine = true;
cfg.runRESTReference = true;
cfg.runHighPassFilter = true;
cfg.importEvents = true;
cfg.runASRWindowRejection = true;
cfg.runICA = true;
cfg.runMARA = true;
cfg.runICLabel = true;
cfg.runEpochChannelInterpolation = true;

% ----------------------- Saving options ---------------------------------
cfg.saveFilteredIntermediate = true;
cfg.saveFinalFiles = true;
cfg.saveProgressEveryNFiles = 100;

% Set to true only for manual inspection after the pipeline finishes.
cfg.previewExampleFileAfterProcessing = false;
cfg.previewFileName = '';

% ----------------------- General quality checks --------------------------
% Minimum number of samples required for processing. Files shorter than this
% are marked as bad and skipped.
cfg.minSamples = 5000;

% ----------------------- Bad-channel rejection ---------------------------
cfg.initialBadChannelMeasure = 'kurt';
cfg.initialBadChannelThreshold = 6;

% ----------------------- CleanLine parameters ----------------------------
cfg.lineNoiseFrequencyHz = 50;
cfg.cleanLineBandwidth = 2;
cfg.cleanLinePValue = 0.01;
cfg.cleanLinePaddingFactor = 2;
cfg.cleanLineTau = 100;
cfg.cleanLineWindowSizeSec = 4;
cfg.cleanLineWindowStepSec = 1;

% ----------------------- High-pass filter parameters ---------------------
% The original script used EEGLAB's reverse-filter syntax:
% pop_eegfiltnew(EEG, [], highPassHz, filterOrder, 1, [], 0).
cfg.highPassHz = 1.5;
cfg.highPassFilterOrder = 826;

% ----------------------- Event handling ---------------------------------
cfg.eventFileExtension = '.txt';
cfg.eventFields = {'latency' 'duration' 'type'};
cfg.eventSkipLines = 1;
cfg.eventTimeUnit = 1;
cfg.endEventType = 'Koniec';

% If the imported event file does not end with cfg.endEventType, an end event
% is added this many seconds before the end of the recording.
cfg.addMissingEndEventSecondsBeforeEnd = 15;

% If an existing end event occurs earlier than this number of samples before
% the data end, the file is marked for event checking.
cfg.endEventMustBeWithinLastSamples = 5000;

% Events of these types are removed from the final continuous file.
cfg.eventTypesToDelete = {'X', 'boundary'};

% ----------------------- ASR / clean_rawdata parameters ------------------
cfg.asrWindowCriterion = 0.55;
cfg.asrWindowCriterionTolerances = [-Inf 20];

% ----------------------- ICA / MARA / ICLabel parameters -----------------
cfg.icaType = 'runica';
cfg.icaExtended = 1;
cfg.minComponentsToKeep = 3;

% ICLabel class thresholds used by pop_icflag.
% Rows correspond to: Brain, Muscle, Eye, Heart, Line Noise, Channel Noise,
% Other. With default settings, the script rejects Eye components with probability >= 0.5.
cfg.iclabelFlagThresholds = [ ...
    NaN NaN; ...  % Brain
    NaN NaN; ...  % Muscle
    0.5 1;   ...  % Eye
    NaN NaN; ...  % Heart
    NaN NaN; ...  % Line noise
    NaN NaN; ...  % Channel noise
    NaN NaN  ...  % Other
    ];

% ----------------------- FASTER epoch-level interpolation ----------------
cfg.epochLengthSec = 1;
cfg.fasterRejectionOptions.measure = [1 1 1 1];
cfg.fasterRejectionOptions.z = [4 4 3 3];

% A channel is interpolated within an epoch only if at least this many FASTER
% measures exceed threshold.
cfg.minExceededFASTERMeasures = 2;

% Epoch-level interpolation is performed only when at least this many channels
% remain good in the epoch.
cfg.minGoodChannelsPerEpoch = 5;

%% ======================= INITIALIZATION =================================

% Add EEGLAB and plugin folders to the MATLAB path.
if ~isempty(cfg.eeglabPath)
    addpath(cfg.eeglabPath);
end

for p = 1:numel(cfg.extraPluginPaths)
    if exist(cfg.extraPluginPaths{p}, 'dir') == 7
        addpath(cfg.extraPluginPaths{p});
    end
end

% Create output folders if they do not exist yet.
if exist(cfg.outputFinalDir, 'dir') ~= 7
    mkdir(cfg.outputFinalDir);
end

if cfg.saveFilteredIntermediate && exist(cfg.outputFilteredDir, 'dir') ~= 7
    mkdir(cfg.outputFilteredDir);
end

% List raw EDF files to process.
rawFileList = dir(fullfile(cfg.rawDataDir, cfg.rawFilePattern));
nFiles = numel(rawFileList);

% Log columns:
%   1. File base name.
%   2. Initially interpolated bad channels.
%   3. Samples removed by ASR window rejection.
%   4. Removed ICA components.
%   5. Epoch-level interpolated channels.
%   6. Error message, if processing failed.
preprocessingLog = cell(nFiles, 6);
badFiles = {};
filesRequiringEventCheck = {};

%% ======================= MAIN PREPROCESSING LOOP ========================

for fileIdx = 1:nFiles

    isBadFile = false;
    errorMessage = '';

    rawEdfFileName = rawFileList(fileIdx).name;

    % Create a base name for output .set files and event-file matching.
    if ~isempty(cfg.rawFileSuffixToRemove) && ...
            numel(rawEdfFileName) >= numel(cfg.rawFileSuffixToRemove) && ...
            strcmp(rawEdfFileName(end-numel(cfg.rawFileSuffixToRemove)+1:end), cfg.rawFileSuffixToRemove)
        fileBaseName = rawEdfFileName(1:end-numel(cfg.rawFileSuffixToRemove));
    else
        [~, fileBaseName, ~] = fileparts(rawEdfFileName);
    end

    preprocessingLog{fileIdx, 1} = fileBaseName;

    fprintf('Processing file %d/%d: %s\n', fileIdx, nFiles, rawEdfFileName);

    try
        %% ----------------------- Load EDF file --------------------------
        rawEdfPath = fullfile(cfg.rawDataDir, rawEdfFileName);

        % Import the EDF signal only. Events and annotations are imported
        % later from external text files if cfg.importEvents = true.
        EEG = pop_biosig(rawEdfPath, 'importevent', 'off', 'importannot', 'off');

        % Skip very short recordings.
        nSamples = size(EEG.data, 2);
        if nSamples <= cfg.minSamples
            isBadFile = true;
            errorMessage = 'Recording shorter than cfg.minSamples.';
        end

        if ~isBadFile
            %% ----------------------- Channel selection ------------------
            currentChannelLabels = cell(1, EEG.nbchan);
            for chanIdx = 1:EEG.nbchan
                currentChannelLabels{chanIdx} = EEG.chanlocs(chanIdx).labels;
            end

            % If cfg.rawChannelLabelsToKeep is empty, all available channels
            % are retained. Otherwise, the script checks whether all required
            % channels are present and selects them in the user-defined order.
            if ~isempty(cfg.rawChannelLabelsToKeep)
                missingChannels = setdiff(cfg.rawChannelLabelsToKeep, currentChannelLabels);

                if isempty(missingChannels)
                    EEG = pop_select(EEG, 'channel', cfg.rawChannelLabelsToKeep);
                else
                    isBadFile = true;
                    errorMessage = ['Missing required channels: ' strjoin(missingChannels, ', ')];
                end
            end
        end

        if ~isBadFile
            %% ----------------------- Rename channels --------------------
            % Remove the user-defined prefix from channel labels, for example
            % "EEG Fp1" -> "Fp1". This makes labels compatible with standard
            % 10-20 / 10-10 channel-location files.
            if ~isempty(cfg.channelLabelPrefixToRemove)
                for chanIdx = 1:EEG.nbchan
                    oldChannelLabel = EEG.chanlocs(chanIdx).labels;
                    newChannelLabel = strrep(oldChannelLabel, cfg.channelLabelPrefixToRemove, '');
                    EEG = pop_chanedit(EEG, 'changefield', {chanIdx, 'labels', newChannelLabel}, 'rplurchanloc', 1);
                end
            end

            %% ----------------------- Add channel locations --------------
            if cfg.addChannelLocations
                EEG = pop_chanedit(EEG, 'lookup', cfg.channelLocationFile, 'rplurchanloc', 1);
            end

            %% ----------------------- Global bad-channel interpolation ----
            % Detect channels with abnormal kurtosis and interpolate them
            % using spherical spline interpolation.
            initiallyBadChannelIdx = [];
            if cfg.runInitialBadChannelInterpolation
                [~, initiallyBadChannelIdx] = pop_rejchan( ...
                    EEG, ...
                    'threshold', cfg.initialBadChannelThreshold, ...
                    'norm', 'on', ...
                    'measure', cfg.initialBadChannelMeasure);

                if ~isempty(initiallyBadChannelIdx)
                    EEG = eeg_interp(EEG, initiallyBadChannelIdx);
                end
            end
            preprocessingLog{fileIdx, 2} = initiallyBadChannelIdx;

            %% ----------------------- Clean line noise -------------------
            % Remove narrow-band electrical line noise from all currently
            % selected channels.
            if cfg.runCleanLine
                EEG = pop_cleanline( ...
                    EEG, ...
                    'bandwidth', cfg.cleanLineBandwidth, ...
                    'chanlist', 1:EEG.nbchan, ...
                    'computepower', 0, ...
                    'linefreqs', cfg.lineNoiseFrequencyHz, ...
                    'normSpectrum', 0, ...
                    'p', cfg.cleanLinePValue, ...
                    'pad', cfg.cleanLinePaddingFactor, ...
                    'plotfigures', 0, ...
                    'scanforlines', 0, ...
                    'sigtype', 'Channels', ...
                    'tau', cfg.cleanLineTau, ...
                    'verb', 0, ...
                    'winsize', cfg.cleanLineWindowSizeSec, ...
                    'winstep', cfg.cleanLineWindowStepSec);
            end

            %% ----------------------- REST re-reference ------------------
            % REST requires an average-referenced signal before applying the
            % reference electrode standardization transform.
            if cfg.runRESTReference
                EEG = pop_reref(EEG, []);

                restChannelIdx = 1:EEG.nbchan;
                leadfieldMatrix = dong_getleadfield(EEG, restChannelIdx, []);
                restReferencedData = rest_refer(EEG.data, leadfieldMatrix);
                EEG.data = restReferencedData;
            end

            %% ----------------------- High-pass filtering ----------------
            if cfg.runHighPassFilter
                EEG = pop_eegfiltnew(EEG, [], cfg.highPassHz, cfg.highPassFilterOrder, 1, [], 0);
            end

            %% ----------------------- Import and validate events ---------
            if cfg.importEvents
                eventFilePath = fullfile(cfg.eventDir, [fileBaseName cfg.eventFileExtension]);

                if isfile(eventFilePath)
                    EEG = pop_importevent( ...
                        EEG, ...
                        'event', eventFilePath, ...
                        'fields', cfg.eventFields, ...
                        'skipline', cfg.eventSkipLines, ...
                        'timeunit', cfg.eventTimeUnit);

                    if ~isempty(EEG.event) && strcmp(EEG.event(end).type, cfg.endEventType)
                        % The last event should indicate the end of the task
                        % and should occur close to the end of the signal.
                        if EEG.event(end).latency > size(EEG.data, 2) - cfg.endEventMustBeWithinLastSamples
                            if cfg.saveFilteredIntermediate
                                EEG = pop_saveset( ...
                                    EEG, ...
                                    'filename', [fileBaseName '.set'], ...
                                    'filepath', cfg.outputFilteredDir);
                            end
                        else
                            isBadFile = true;
                            filesRequiringEventCheck{end + 1} = fileBaseName;
                            errorMessage = 'End event is too far from the end of the recording.';
                        end
                    else
                        % If the event file does not contain the expected end
                        % marker, add it close to the end of the recording.
                        newEndLatencySec = size(EEG.data, 2) / EEG.srate - cfg.addMissingEndEventSecondsBeforeEnd;

                        EEG = pop_editeventvals( ...
                            EEG, ...
                            'append', {1 [] [] [] [] [] [] []}, ...
                            'changefield', {2, 'latency', newEndLatencySec}, ...
                            'changefield', {2, 'duration', 0}, ...
                            'changefield', {2, 'type', cfg.endEventType});
                        EEG = eeg_checkset(EEG);

                        if cfg.saveFilteredIntermediate
                            EEG = pop_saveset( ...
                                EEG, ...
                                'filename', [fileBaseName '.set'], ...
                                'filepath', cfg.outputFilteredDir);
                        end
                    end
                else
                    isBadFile = true;
                    errorMessage = 'Event file not found.';
                end
            end
        end

        if ~isBadFile
            %% ----------------------- ASR window rejection ----------------
            % This step removes the worst contaminated time windows using the
            % clean_rawdata window-rejection criterion. Other clean_rawdata
            % criteria are switched off to preserve the behavior of the
            % original script.
            removedSampleIdx = [];

            if cfg.runASRWindowRejection
                eventsBeforeASR = EEG.event;
                EEGBeforeASR = EEG; %#ok<NASGU>  % Kept for optional visual inspection with vis_artifacts.

                EEG = pop_clean_rawdata( ...
                    EEG, ...
                    'FlatlineCriterion', 'off', ...
                    'ChannelCriterion', 'off', ...
                    'LineNoiseCriterion', 'off', ...
                    'Highpass', 'off', ...
                    'BurstCriterion', 'off', ...
                    'WindowCriterion', cfg.asrWindowCriterion, ...
                    'BurstRejection', 'off', ...
                    'Distance', 'Euclidian', ...
                    'WindowCriterionTolerances', cfg.asrWindowCriterionTolerances);

                removedSampleIdx = find(EEG.etc.clean_sample_mask == 0);

                % Reinsert events that were lost because their original
                % latency fell inside a removed segment. This preserves the
                % logic of the source script.
                for eventIdx = 1:numel(eventsBeforeASR)
                    removedEventPosition = find(removedSampleIdx == eventsBeforeASR(eventIdx).latency);

                    if ~isempty(removedEventPosition)
                        newEventLatencySec = (eventsBeforeASR(eventIdx).latency - removedEventPosition) / EEG.srate;

                        if newEventLatencySec < size(EEG.data, 2) / EEG.srate && numel(eventsBeforeASR) > 3
                            try
                                EEG = pop_editeventvals( ...
                                    EEG, ...
                                    'append', {1 [] [] [] [] [] [] []}, ...
                                    'changefield', {2, 'latency', newEventLatencySec}, ...
                                    'changefield', {2, 'duration', 0}, ...
                                    'changefield', {2, 'type', eventsBeforeASR(eventIdx).type});
                                EEG = eeg_checkset(EEG);
                            catch
                                isBadFile = true;
                                errorMessage = 'Could not restore event after ASR rejection.';
                            end
                        else
                            isBadFile = true;
                            errorMessage = 'Restored event latency after ASR is outside the recording.';
                        end
                    end
                end
            end

            preprocessingLog{fileIdx, 3} = removedSampleIdx;
        end

        if ~isBadFile
            %% ----------------------- Optional ICA component rejection -----
            % ICA itself, MARA, and ICLabel can be enabled independently.
            % MARA and ICLabel require ICA activations, so they are used only
            % inside this block when cfg.runICA = true.
            artifactComponentIdx = [];

            if cfg.runICA
                eventsBeforeICA = EEG.event;

                EEG = pop_runica( ...
                    EEG, ...
                    'icatype', cfg.icaType, ...
                    'extended', cfg.icaExtended, ...
                    'interrupt', 'on');

                maraArtifactComponentIdx = [];
                eyeComponentIdx = [];

                if cfg.runMARA
                    [maraArtifactComponentIdx, MARAinfo] = MARA(EEG); %#ok<NASGU>
                end

                if cfg.runICLabel
                    % ICLabel is used here only to flag eye-related components,
                    % following the logic of the source script.
                    EEGForICLabel = EEG;
                    EEGForICLabel = iclabel(EEGForICLabel);
                    EEGForICLabel = pop_icflag(EEGForICLabel, cfg.iclabelFlagThresholds);
                    eyeComponentIdx = find(EEGForICLabel.reject.gcompreject > 0);
                end

                artifactComponentIdx = unique([maraArtifactComponentIdx(:); eyeComponentIdx(:)])';
                preprocessingLog{fileIdx, 4} = artifactComponentIdx;

                % Require at least cfg.minComponentsToKeep components to remain
                % after artifact-component removal.
                if numel(artifactComponentIdx) > EEG.nbchan - cfg.minComponentsToKeep
                    isBadFile = true;
                    errorMessage = 'Too many ICA components marked as artifacts.';
                else
                    if ~isempty(artifactComponentIdx)
                        EEG = pop_subcomp(EEG, artifactComponentIdx, 0);
                    end

                    % Restore the event structure saved before ICA.
                    EEG.event = eventsBeforeICA;
                end
            else
                preprocessingLog{fileIdx, 4} = artifactComponentIdx;

                if cfg.runMARA || cfg.runICLabel
                    warning(['cfg.runMARA or cfg.runICLabel is true, but cfg.runICA is false. ' ...
                        'MARA and ICLabel were skipped for this file.']);
                end
            end
        end

        if ~isBadFile
            %% ----------------------- Optional FASTER epoch interpolation --
            % Create regular 1-second epochs, detect bad channels separately
            % within each epoch using FASTER properties, interpolate those
            % channels, and then return to continuous data.
            epochInterpolatedChannels = {};

            if cfg.runEpochChannelInterpolation
                EEG = eeg_regepochs(EEG, 'recurrence', cfg.epochLengthSec, 'rmbase', NaN);

                nChannels = numel(EEG.chanlocs);
                epochInterpolatedChannels = cell(1, EEG.trials);
                epochWasInterpolated = zeros(1, EEG.trials); %#ok<NASGU>

                for epochIdx = 1:EEG.trials
                    epochChannelProperties = single_epoch_channel_properties(EEG, epochIdx, 1:nChannels);
                    exceededThreshold = min_z_changed(epochChannelProperties, cfg.fasterRejectionOptions);

                    badChannelsInEpoch = find(exceededThreshold >= cfg.minExceededFASTERMeasures);

                    if ~isempty(badChannelsInEpoch)
                        epochWasInterpolated(epochIdx) = 1;
                    end

                    % Interpolate only if enough good channels remain in this
                    % epoch. Otherwise, the epoch is left unchanged, matching
                    % the conservative behavior of the source script.
                    if numel(badChannelsInEpoch) < nChannels - cfg.minGoodChannelsPerEpoch
                        epochInterpolatedChannels{epochIdx} = badChannelsInEpoch;
                    else
                        epochInterpolatedChannels{epochIdx} = [];
                    end
                end

                EEG = h_epoch_interp_spl(EEG, epochInterpolatedChannels);

                % Convert the temporary regular epochs back to continuous data.
                EEG = epoch2continuous(EEG);
            end

            preprocessingLog{fileIdx, 5} = epochInterpolatedChannels;

            %% ----------------------- Remove unnecessary events -----------
            eventsToDelete = [];

            for eventIdx = 1:numel(EEG.event)
                if ismember(EEG.event(eventIdx).type, cfg.eventTypesToDelete)
                    eventsToDelete(end + 1) = eventIdx; %#ok<SAGROW>
                end
            end

            if ~isempty(eventsToDelete)
                EEG = pop_editeventvals(EEG, 'delete', eventsToDelete);
            end
        end

        %% ----------------------- Save final file -------------------------
        if ~isBadFile && cfg.saveFinalFiles
            EEG = pop_saveset( ...
                EEG, ...
                'filename', [fileBaseName '.set'], ...
                'filepath', cfg.outputFinalDir);
        end

    catch ME
        isBadFile = true;
        errorMessage = ME.message;
    end

    if isBadFile
        badFiles{end + 1} = rawEdfFileName; %#ok<SAGROW>
    end

    preprocessingLog{fileIdx, 6} = errorMessage;

    % Save progress periodically so partial results are not lost if the script
    % stops during a long batch.
    if cfg.saveProgressEveryNFiles > 0 && rem(fileIdx, cfg.saveProgressEveryNFiles) == 0
        save(fullfile(cfg.outputFinalDir, 'preprocessing_log_partial.mat'), ...
            'preprocessingLog', 'badFiles', 'filesRequiringEventCheck', 'cfg');
    end
end

%% ======================= SAVE FINAL LOG =================================

% Backward-compatible aliases for variables used in the original script.
epo_removed = preprocessingLog;
bad_files = badFiles;
files_check_events = filesRequiringEventCheck;

save(fullfile(cfg.outputFinalDir, 'preprocessing_log.mat'), ...
    'preprocessingLog', 'badFiles', 'filesRequiringEventCheck', ...
    'epo_removed', 'bad_files', 'files_check_events', 'cfg');

%% ======================= OPTIONAL VISUAL CHECK ==========================

if cfg.previewExampleFileAfterProcessing
    if isempty(cfg.previewFileName)
        finalSetFiles = dir(fullfile(cfg.outputFinalDir, '*.set'));

        if ~isempty(finalSetFiles)
            cfg.previewFileName = finalSetFiles(1).name;
        else
            error('No .set files were found in cfg.outputFinalDir for preview.');
        end
    end

    EEG = pop_loadset(cfg.previewFileName, cfg.outputFinalDir);
    pop_eegplot(EEG, 1, 1, 1);
end
