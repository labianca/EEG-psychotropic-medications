function [matched_indices, chi_stat_all] = balance_groups_meds(metadata, data_meds, nmatch, nmatch2)
%==========================================================================
% BALANCE_GROUPS_MEDS FUNCTION
%
% Creates matched participant indices for medication-group comparisons.
%
% For each medication class, the function compares:
%   - participants taking the medication class, and
%   - participants not taking the medication class.
%
% The smaller group is used as the reference group. For each participant in
% the reference group, the function searches for matching participants in the
% larger group based on:
%   1. sex,
%   2. age rounded to 5-year bins,
%   3. diagnosis.
%
% If several candidate matches are available, the function repeatedly samples
% possible matched sets and selects the set with the best balance in the
% distribution of other medication classes, estimated using a chi-square
% goodness-of-fit statistic.
%
% INPUTS
%   metadata : matrix of participant metadata, with rows corresponding to
%              participants. Expected columns:
%              column 1 = sex
%              column 2 = age
%              column 3 = diagnosis
%
%   data_meds : binary matrix of medication labels, with rows corresponding
%               to participants and columns corresponding to medication
%               classes. A value of 1 indicates that the participant takes
%               the medication class; 0 indicates that they do not.
%
%   nmatch  : number of independent matched samples to generate for each
%             medication class.
%
%   nmatch2 : number of random candidate matched sets tested within each
%             matching iteration. Larger values may produce better-balanced
%             matched groups, but increase computation time.
%
% OUTPUT
%   matched_indices : cell array such that matched_indices{l}{i} contains
%                     matched participant indices for medication class l and
%                     matching repetition i.
%
%                     The first column contains indices from the medication
%                     group, and the second column contains indices from
%                     the matched comparison group.
%
%   chi_stat_all :    goodness of fit statistics for diagnostic purposes 

%==========================================================================

% Initialise output cell array.
matched_indices = {};

% Store chi-square statistics for candidate matches.
chi_stat_all = {};

% Loop over medication classes.
for med_idx = 1:size(data_meds, 2)

    % Continue only if at least one participant takes this medication class.
    if sum(data_meds(:, med_idx)) > 0

        % Identify participants who do not take and who do take this
        % medication class.
        idx_nonmed = find(data_meds(:, med_idx) == 0);
        idx_med    = find(data_meds(:, med_idx) == 1);

        % The matching procedure selects controls for the smaller group.
        % If the medication group is larger than the non-medication group,
        % reverse the groups so that the smaller group is always used as
        % the reference group.
        groups_reversed = false;

        if numel(idx_med) > numel(idx_nonmed)
            idx_nonmed_backup = idx_nonmed;
            idx_nonmed = idx_med;
            idx_med = idx_nonmed_backup;
            groups_reversed = true;
        end

        % Create a matrix used for matching and for balancing the
        % distribution of other medication classes.
        %
        % Rows correspond to participants; columns include:
        %   1. sex,
        %   2. age,
        %   3:end medication-class indicators.
        matching_predictors = [metadata(:, 1), metadata(:, 2), data_meds];

        % Extract predictors for the two comparison groups.
        % The transposition is used because the later indexing expects variables in rows and participants in columns.
        predictors_nonmed = matching_predictors(idx_nonmed, :)';
        predictors_med    = matching_predictors(idx_med, :)';

        % Extract diagnosis labels for both groups.
        diagnosis_nonmed = metadata(idx_nonmed, 3);
        diagnosis_med    = metadata(idx_med, 3);

        % Build simplified matching descriptors:
        %   column 1 = sex,
        %   column 2 = age rounded into 5-year bins,
        %   column 3 = diagnosis.
        match_descriptors_nonmed = [ ...
            predictors_nonmed(1, :)', ...
            round(predictors_nonmed(2, :)' / 5), ...
            diagnosis_nonmed];

        match_descriptors_med = [ ...
            predictors_med(1, :)', ...
            round(predictors_med(2, :)' / 5), ...
            diagnosis_med];

        % Preallocate chi-square statistics for this medication class.
        chi_stat_all{med_idx} = zeros(nmatch, 6);

        % Generate several independent matched samples for this medication
        % class. These repetitions can later be used to estimate the
        % stability of model coefficients or effect sizes across matchings.
        for match_rep = 1:nmatch

            % Initialise the output for this medication class and repetition.
            matched_indices{med_idx}{match_rep} = [];

            % Candidate matches for each participant in the smaller group.
            candidate_matches = {};

            % Number of candidate matches available for each participant.
            n_candidate_matches = zeros(numel(idx_med), 1);

            % Count participants for whom no diagnosis-matched control was
            % found. This value is not used later, but is useful for
            % inspecting whether matching failed often.
            n_unmatched = 0;

            % Search for possible matches for each participant in the
            % smaller/reference group.
            for participant_idx = 1:numel(idx_med)

                % Descriptor of the current participant:
                % sex, 5-year age bin, and diagnosis.
                participant_descriptor = [ ...
                    predictors_med(1, participant_idx), ...
                    round(predictors_med(2, participant_idx) / 5), ...
                    diagnosis_med(participant_idx)];

                % First, try to find exact matches on sex, age bin, and
                % diagnosis in the larger/comparison group.
                exact_matches = find(sum(match_descriptors_nonmed == participant_descriptor, 2) == 3);

                if numel(exact_matches) >= 1

                    % Store all exact candidate matches. The final match is
                    % selected later while also trying to balance other
                    % medication classes and avoid repeated participants.
                    candidate_matches{participant_idx} = exact_matches;
                    n_candidate_matches(participant_idx) = numel(exact_matches);

                else
                    % If no exact match is available, try to find
                    % participants with the same diagnosis and then choose
                    % the one closest in sex and age-bin values.
                    diagnosis_matches = find(participant_descriptor(3) == match_descriptors_nonmed(:, 3));

                    if numel(diagnosis_matches) >= 1

                        % Compute distance based on sex and age-bin values.
                        similarity_distance = sum( ...
                            abs(participant_descriptor(:, 1:2) - match_descriptors_nonmed(diagnosis_matches, 1:2)), ...
                            2);

                        % Retain candidates with the minimum distance.
                        nearest_matches = diagnosis_matches(similarity_distance == min(similarity_distance));

                        % Randomly select one of the nearest candidates.
                        random_idx = randperm(numel(nearest_matches), 1);
                        selected_match = nearest_matches(random_idx);

                        candidate_matches{participant_idx} = selected_match;
                        n_candidate_matches(participant_idx) = 1;

                    else
                        % No diagnosis-matched participant was available.
                        % This case should be inspected if it occurs often.
                        n_unmatched = n_unmatched + 1;
                        candidate_matches{participant_idx} = [];
                        n_candidate_matches(participant_idx) = 0;
                    end
                end
            end

            % Participants with fewer candidate matches are processed first
            % so they have a better chance of receiving a unique match.
            [~, sorted_participant_idx] = sort(n_candidate_matches);

            % Store candidate matched sets.
            candidate_matched_sets = zeros(nmatch2, numel(sorted_participant_idx));

            % Store the chi-square statistic for each candidate matched set.
            chi_square_values = zeros(nmatch2, 1);

            % Generate nmatch2 random matched sets and evaluate each one.
            for random_rep = 1:nmatch2

                % Select one candidate match for each reference participant.
                for participant_idx = sorted_participant_idx'

                    if ~isempty(candidate_matches{participant_idx})

                        random_idx = randperm(numel(candidate_matches{participant_idx}), 1);
                        n_attempts = 0;

                        % Try to avoid selecting the same participant more
                        % than once within the same matched set.
                        while ismember(candidate_matches{participant_idx}(random_idx), candidate_matched_sets(random_rep, :)) && ...
                                n_attempts < numel(candidate_matches{participant_idx})
                            random_idx = randperm(numel(candidate_matches{participant_idx}), 1);
                            n_attempts = n_attempts + 1;
                        end

                        candidate_matched_sets(random_rep, participant_idx) = ...
                            candidate_matches{participant_idx}(random_idx);

                    else
                        % If no candidate was found, assign a placeholder.
                        % The placeholder will be corrected below by adding
                        % random additional participants if necessary.
                        candidate_matched_sets(random_rep, participant_idx) = 1;
                    end
                end

                % Remove duplicate selected participants.
                selected_controls = unique(candidate_matched_sets(random_rep, :));

                % If duplicates caused the number of selected controls to be
                % too small, randomly add unused participants from the larger
                % group until the selected set has the required size.
                while numel(unique(selected_controls)) < numel(idx_med)

                    available_indices = setdiff(1:numel(idx_nonmed), selected_controls);
                    n_to_add = numel(idx_med) - numel(unique(selected_controls));

                    selected_controls = [ ...
                        unique(selected_controls), ...
                        available_indices(randperm(numel(available_indices), n_to_add))];
                end

                candidate_matched_sets(random_rep, :) = selected_controls;

                % Extract medication indicators for the selected controls.
                % Rows correspond to medication classes; columns correspond
                % to selected participants.
                selected_medication_profiles = predictors_nonmed( ...
                    3:size(predictors_nonmed, 1), ...
                    candidate_matched_sets(random_rep, :));

                % Build a vector listing medication classes present in the
                % selected control set. This is used to evaluate whether the
                % selected controls are balanced with respect to other
                % medications.
                medication_distribution = [];

                for selected_idx = 1:numel(candidate_matched_sets(random_rep, :))
                    active_medications = find(selected_medication_profiles(:, selected_idx) == 1);
                    medication_distribution = [medication_distribution; active_medications];
                end

                % Add all medication-class labels once so that expected bins
                % are represented in the chi-square test.
                medication_distribution = [ ...
                    medication_distribution; ...
                    (1:(size(matching_predictors, 2) - 2))'];

                % Define a simple approximately uniform expected
                % distribution across medication classes.
                expected_counts = ...
                    (int32(numel(candidate_matched_sets(random_rep, :)) / (size(matching_predictors, 2) - 2)) + 1) * ...
                    int32(ones(1, size(matching_predictors, 2) - 2));

                % Evaluate how close the distribution of medication classes
                % in this candidate matched set is to the expected
                % distribution. Lower chi-square values are preferred.
                [~, ~, chi2stat] = chi2gof( ...
                    medication_distribution, ...
                    'Expected', double(expected_counts));

                chi_square_values(random_rep) = chi2stat.chi2stat;
            end

            % Replace NaN chi-square values with 0
            chi_square_values(isnan(chi_square_values)) = 0;

            % Sort candidate matched sets by chi-square statistic and choose
            % the most balanced one.
            [~, sorted_candidate_idx] = sort(chi_square_values);
            best_match_idx = sorted_candidate_idx(1);

            % Try to avoid selecting a candidate matched set that overlaps
            % too strongly with the match already stored for the same
            % medication class and repetition.

            n_attempts_best = 1;
            previous_matches = unique(matched_indices{med_idx}{match_rep});

            while numel(intersect(unique(idx_nonmed(candidate_matched_sets(best_match_idx, :))), previous_matches)) > ...
                    size(match_descriptors_med, 1) - 2 && ...
                    n_attempts_best < min(12, numel(sorted_candidate_idx))

                best_match_idx = sorted_candidate_idx(1 + n_attempts_best);
                n_attempts_best = n_attempts_best + 1;
            end

            % Select the best candidate matched set.
            selected_controls = candidate_matched_sets(best_match_idx, :);

            % Diagnostic check: the final number of selected controls should
            % match the number of participants in the reference group.
            if numel(selected_controls) ~= size(match_descriptors_med, 1)
                warning(['The number of selected controls does not match ', ...
                         'the number of reference participants for medication ', ...
                         'class %d, repetition %d.'], med_idx, match_rep);
            end

            % Store matched indices using original indices within the
            % input subset. If the groups were reversed above, the output is
            % reversed accordingly so exposed patients are always in the first column.
            if groups_reversed
                matched_indices{med_idx}{match_rep} = [idx_nonmed(selected_controls), idx_med];
            else
                matched_indices{med_idx}{match_rep} = [idx_med, idx_nonmed(selected_controls)];
            end
        end

    else
        % If no participant takes this medication class, return empty
        % matches for all requested repetitions.
        for match_rep = 1:nmatch
            matched_indices{med_idx}{match_rep} = [];
        end
    end
end

end
