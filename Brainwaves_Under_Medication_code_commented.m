%==========================================================================
% MEDICATION EFFECT ANALYSIS IN RESTING-STATE EEG
%
% Purpose:
%   This script analyses associations between psychotropic medication groups
%   and resting-state EEG features. It prepares the data, performs PCA-based
%   dimensionality reduction, creates matched medication comparison groups 
%  (for medication versus other drugs [OvR] and medication versus drug-naive patients [DN] comparison schemes),
%   fits linear mixed-effects models for each retained principal component,
%   computes effect sizes, applies multiple-comparison correction, and saves
%   the resulting statistics.
%
% Expected input:
%   An HDF5 file containing:
%       /data       - EEG features, size: N_patients x N_features
%       /data_meds  - binary medication labels, size: N_patients x N_med_groups
%       /diagnosis  - numeric diagnostic labels, size: N_patients x 1
%       /age        - patient age, size: N_patients x 1
%       /sex        - numeric sex code, size: N_patients x 1
%       /year       - date (year) of EEG recording, size: N_patients x 1
%       /site       - recording/hospital site, numeric, size: N_patients x 1
%       /drug_naive - OPTIONAL. binary labels indicating which patients do not take any psychotropic drugs, even beyond classes in data_meds, size: N_patients x 1
%
% Expected output:
%   - Cleaned EEG and medication data saved as HDF5
%   - PCA decomposition saved as MAT file
%   - Matched participant indices for both comparison schemes
%   - P-values, corrected p-values, t-statistics, and effect sizes
%   - Linear mixed-effects model coefficients and model-level statistics
%
% Notes for reuse:
%   - All input variables should be numeric.
%   - Missing values should be encoded as NaN.
%   - The custom functions balance_groups_meds, balance_groups_meds_DN, and bonf_holm 
%    (available in the folder helpers) must be available on the MATLAB path.
%==========================================================================

%% ---------- Set analysis parameters ----------
% These parameters define the main analysis settings. The current values
% correspond to the settings used in the publication analysis.
exclude_multidrug = 0;   % Exclude patients taking drugs from more than one medication group?
                         % 0 = keep all patients; 1 = remove multi-drug patients.

nmatch = 10;             % Number of independent matched comparisons to generate.
                         % Multiple matchings can be used to assess stability of results.

nmatch2 = 500;           % Number of random candidate matchings evaluated inside the
                         % balancing function. Higher values can improve matching quality
                         % but increase computation time.

min_var_expl = 90;       % Minimum cumulative variance explained by retained PCA components.
                         % The script retains the smallest number of PCs explaining at least
                         % this percentage of variance.

SD_threshold = 20;       % Outlier threshold in standard deviation units.
                         % Participants with any z-scored feature exceeding this threshold
                         % are removed from further analysis.

%% ---------- Define input and output paths ----------
% Folder containing the HDF5 input file.
datapath = 'C:\Users\mszponar\Documents\kody itp\classification2\medicines\data_meds_REST_ASR_IPiN_i_Nowo\DANE_popr';

% Folder where outputs will be saved.
out_path = 'C:\Users\mszponar\Documents\kody itp\classification2\medicines\data_meds_REST_ASR_IPiN_i_Nowo\DANE_popr\out_folder_test';

% Name of the HDF5 file containing EEG features, medication labels, and metadata.
datafname = 'datafile_name.h5';


%% ---------- Load data from HDF5 file ----------
% Construct the full path to the HDF5 data file.
datafile = fullfile(datapath, datafname);

% Load medication-group labels.
% Each column corresponds to one medication group; values should be binary
% indicators, e.g. 1 = patient takes medication from this group, 0 = no.
data_meds = h5read(datafile, '/data_meds');

% Load EEG feature matrix.
% Rows correspond to patients; columns correspond to EEG features.
data = h5read(datafile, '/data');

% Load patient metadata.
diagnosis = h5read(datafile, '/diagnosis');  % Numeric diagnostic labels.
age       = h5read(datafile, '/age');        % Patient age.
sex       = h5read(datafile, '/sex');        % Numeric sex code.
year      = h5read(datafile, '/year');       % Year of EEG recording.
site      = h5read(datafile, '/site');       % Recording/hospital site.

% Load drug-naive labels if available. Otherwise, derive them from the
% medication matrix. Drug-naive participants are then defined as those who do
% not take any medication group listed in data_meds.

info = h5info(datafile);
dataset_names = {info.Datasets.Name};

if ismember('drug_naive', dataset_names)
  % Read pre-computed drug-naive labels from the HDF5 file.
  drug_naive = h5read(datafile, '/drug_naive'); 
else
  % Create drug-naive labels from medication data. Participants with no medication group assigned are considered drug-naive.
  drug_naive = sum(data_meds, 2) == 0
end

%% ---------- Prepare metadata ----------
% Centre the recording year around its mean. This improves interpretability
% of the regression intercept and reduces unnecessary scaling differences
% between predictors.
year_c = year - nanmean(year);

% Combine metadata variables into one matrix for easier indexing later.
% Column order:
%   1 = sex
%   2 = age
%   3 = diagnosis
%   4 = centred recording year
%   5 = site
metadata = horzcat(sex, age, diagnosis, year_c, site);

%% ---------- Remove participants with missing data ----------
% Identify participants with missing metadata.
missing1 = find(sum(isnan(metadata), 2) > 0);

% Identify participants with missing medication labels.
missing2 = find(sum(isnan(data_meds), 2) > 0);

% Identify participants with missing EEG features.
missing3 = find(any(isnan(data), 2));

% Combine all missing-data indices and remove duplicates.
missing = unique([missing1; missing2; missing3]);

% Remove participants with missing data from all aligned matrices.
data(missing, :) = [];
data_meds(missing, :) = [];
metadata(missing, :) = [];

%% ---------- Z-score EEG features and remove feature outliers ----------
% Standardise EEG features across participants.
% zdata has the same size as data, but each feature is centred and scaled.
[zdata, mu, sigma] = zscore(data);

% Identify and remove participants with at least one feature above the outlier threshold.
rowsToDelete = find(any(abs(zdata) > SD_threshold, 2));
data(rowsToDelete, :) = [];
data_meds(rowsToDelete, :) = [];
metadata(rowsToDelete, :) = [];

% Recompute z-scored features after outlier removal.
% This ensures PCA is performed on the final cleaned dataset.
[zdata, mu, sigma] = zscore(data);

%% ---------- Perform principal component analysis ----------
% PCA is performed on standardised EEG features to reduce dimensionality.
% PCA_score contains the principal component scores used later as dependent
% variables in the mixed-effects models.
[PCA_coeff, PCA_score, latent, tsquared, explained, mu2] = pca(zdata);

%% ---------- Save cleaned data and PCA results ----------
% Create a subfolder for cleaned data and PCA output.
mkdir(fullfile(out_path, 'cleaned_data'))

% Define the output HDF5 file for cleaned data.
new_data = fullfile(out_path, 'cleaned_data', 'data_new.h5');

% Save cleaned raw EEG features.
h5create(new_data, "/data", size(data))
h5write(new_data, "/data", data)

% Save cleaned and z-scored EEG features.
h5create(new_data, "/zdata", size(zdata))
h5write(new_data, "/zdata", zdata)

% Save cleaned medication labels.
h5create(new_data, "/data_meds", size(data_meds))
h5write(new_data, "/data_meds", data_meds)

% Save PCA decomposition and summary information.
save(fullfile(out_path, 'cleaned_data', 'PCA_data.mat'), ...
    'PCA_coeff', 'PCA_score', 'latent', 'tsquared', 'explained', 'mu2', '-v7.3')

%% ---------- Optional exclusion of multi-drug patients ----------
% For sensitivity analysis, remove patients taking medications from more
% than one medication class. This step is controlled by exclude_multidrug.
if exclude_multidrug == 1
    % Identify patients with more than one active medication-group label.
    multidrug = find(sum(data_meds, 2) > 1);

    % Display the number of excluded patients.
    length(multidrug)

    % Remove multi-drug patients from all aligned matrices.
    data(multidrug, :) = [];
    data_meds(multidrug, :) = [];
    metadata(multidrug, :) = [];
end

%% ---------- Create matched comparison groups ----------
% This section creates matched participant indices for two comparison schemes:
%
%   1. OvR: one medication group versus the rest
%   2. DN:  medication group versus participants not receiving medication
%
% Matching is performed separately within each site and then combined across
% sites. This helps preserve site structure and reduce confounding by site.

matched_indices_OvR = {};
matched_indices_DN = {};
matched_indices_by_site = {};
matched_indices_by_site_DN = {};

idx_set = {};
unique_sites = unique(metadata(:, 5));
for s = 1:length(unique_sites)
    % Identify participants from each recording/hospital site.
    idx_set{s} = find(metadata(:, 5) == unique_sites(s));
   
    % Match participants within each site for the one-vs-rest comparison.
    matched_indices_by_site{s}, ~ = balance_groups_meds( ...
        metadata(idx_set{s}, :), ...
        data_meds(idx_set{s}, :), ...
        nmatch, ...
        nmatch2);
    
    % Match participants within each site for the one-vs-drug-naive comparison.
    matched_indices_by_site_DN{s}, ~ = balance_groups_meds_DN( ...
        metadata(idx_set{s}, :), ...
        data_meds(idx_set{s}, :), ...
        drug_naive, ...
        nmatch, ...
        nmatch2);
end

% Convert site-specific indices back to indices in the full dataset and
% concatenate matched participants across all sites.
for l = 1:size(data_meds, 2)
   for i = 1:nmatch
       
       matched_indices_OvR{l}{i} = [];
       matched_indices_DN{l}{i} = [];

       for s = 1:length(unique_sites)
            matched_indices_OvR{l}{i} = vertcat( ...
                matched_indices_OvR{l}{i}, ...
                idx_set{s}(matched_indices_by_site{s}{l}{i}));
        
            matched_indices_DN{l}{i} = vertcat( ...
                matched_indices_DN{l}{i}, ...
                idx_set{s}(matched_indices_by_site_DN{s}{l}{i}));
       end
   end
end


% Save both sets of matched indices.
save(fullfile(out_path, 'matched_indices.mat'), ...
    'matched_indices_OvR', 'matched_indices_DN')

%% ---------- Calculate effect sizes and model statistics ----------
% Retain the minimum number of principal components required to explain at
% least min_var_expl percent of the total variance.
nPCA = find(cumsum(explained) > min_var_expl, 1);

% Run the downstream statistical analysis separately for each comparison
% scheme: OvR and DN.
for comparison = ["OvR", "DN"]

    % Select the matched index structure corresponding to the current
    % comparison scheme.
    if comparison == "OvR"
        matched_indices = matched_indices_OvR;
    elseif comparison == "DN"
        matched_indices = matched_indices_DN;
    end

    % Initialise containers for statistical results.
    P_all = {};          % Uncorrected p-values for each medication group and matching.
    T_all = {};          % t-statistics for group effects.
    Hedges_g_all = {};   % Hedges' g effect sizes.
    Cohens_d_all = {};   % Cohen's d effect sizes.
    P_corr = {};         % Holm-Bonferroni corrected p-values.
    LME_all = {};        % Linear mixed-effects model coefficients and diagnostics.

    %% ---------- Loop over medication groups ----------
    for l = 1:size(data_meds, 2)

        % Preallocate matrices for one medication group.
        % Rows = retained principal components; columns = matching repetitions.
        p_li = zeros(nPCA, nmatch);
        t_li = zeros(nPCA, nmatch);
        h_li = zeros(nPCA, nmatch);
        d_li = zeros(nPCA, nmatch);

        %% ---------- Loop over repeated matched samples ----------
        for i = 1:nmatch

            % Create binary group labels for the current matched sample.
            % The first half of the vector is coded as medication group = 1,
            % and the second half as comparison group = 0.

            group = [ones(length(matched_indices{l}{i}), 1); ...
                     zeros(length(matched_indices{l}{i}), 1)];

            % Extract centred recording year for both matched groups.
            year_c = [metadata(matched_indices{l}{i}(:, 1), 4); ...
                      metadata(matched_indices{l}{i}(:, 2), 4)];

            % Extract site labels for both matched groups.
            site = [metadata(matched_indices{l}{i}(:, 1), 5); ...
                    metadata(matched_indices{l}{i}(:, 2), 5)];

            %% ---------- Loop over retained PCA components ----------
            for c = 1:nPCA

                % Extract PCA component scores for the current matched groups.
                % These scores are the dependent variable in the LME model.
                feature = [PCA_score(matched_indices{l}{i}(:, 1), c); ...
                           PCA_score(matched_indices{l}{i}(:, 2), c)];

                % Build the model table.
                tbl = table(feature, group, site, year_c);

                % Treat group and site as categorical predictors.
                tbl.group = categorical(tbl.group);
                tbl.site  = categorical(tbl.site);

                % Fit a linear mixed-effects model:
                %   fixed effects: medication group, recording year
                %   random effect: site-level intercept
                %
                % The group coefficient estimates the adjusted difference
                % between the medication and comparison groups for this PCA component.
                lme = fitlme(tbl, 'feature ~ group + year_c + (1|site)');

                % Locate the coefficient corresponding to the medication group effect.
                idx_g = strcmp(lme.Coefficients.Name, 'group_1');

                % Extract t-statistic and p-value for the group effect.
                t_group = lme.Coefficients.tStat(idx_g);
                p_group = lme.Coefficients.pValue(idx_g);

                %% ---------- Calculate effect sizes ----------
                % Extract the group coefficient and model residual variance.
                beta = lme.Coefficients.Estimate(idx_g);
                sigma = sqrt(lme.MSE);
                df = lme.Coefficients.DF(idx_g);

                % Cohen's d: standardised group difference using residual SD.
                d = beta / sigma;

                % Hedges' g: Cohen's d corrected for small-sample bias.
                J = 1 - 3 ./ (4 * df - 1);
                g = J .* d;

                % Store effect sizes and test statistics.
                d_li(c, i) = d;
                h_li(c, i) = g;
                p_li(c, i) = p_group;
                t_li(c, i) = t_group;

                % Save model coefficients and selected model-level statistics.
                %This variable contains coefficients also for recording
                %year effect and confidence intervals for beta coefficients
                % Storing MSE and DF allows recalculation of effect sizes later.
                LME_all{l}{i}{c}.coeff = lme.Coefficients;
                LME_all{l}{i}{c}.MSE = lme.MSE;
                LME_all{l}{i}{c}.DF = df;

            end

            %% ---------- Correct for multiple comparisons ----------
            % Correct p-values across retained PCA components using the
            % Holm-Bonferroni procedure to control the family-wise error rate.
            [corrected_p, h] = bonf_holm(p_li(:, i), 0.05);

            % Store significant component indices and summary counts.
            ist_p_bh{l, i} = find(h == 1);
            sum_ist_bh(l, i) = sum(h == 1, 'all');
            P_corr{l, i} = corrected_p;
            sum_p_uncor(l, i) = sum(p_li(:, i) < 0.05, 'all');

        end

        % Store results for the current medication group.
        P_all{l} = p_li;
        Cohens_d_all{l} = d_li;
        Hedges_g_all{l} = h_li;
        T_all{l} = t_li;

    end

    %% ---------- Save statistical results ----------
    % Save p-values, corrected p-values, t-statistics, and effect sizes.
    savename = strcat('Effects_medicines_PCA_', comparison, '.mat');
    save(fullfile(out_path, savename), ...
        'P_all', 'P_corr', 'T_all', 'Hedges_g_all', 'Cohens_d_all')

    % Save linear mixed-effects model coefficients and model diagnostics.
    savename2 = strcat('Regression_model_medicines_PCA_', comparison, '.mat');
    save(fullfile(out_path, savename2), 'LME_all');

end
