%%%
%  MLS-line with 3DFIN


traj = {'G'};

for i = 1:length(traj)
    walk = traj{i};




    %% plots in analysis (44)
    include_ids = [600032, 600031, 600016, 600015, 599969, 599110, 599109, ...
        599107, 599106, 599101, 599099, 599096, 599093, 599090, ...
        599088, 599086, 599084, 599082, 599081, 599079, 599074, ...
        599070, 599059, 599057, 599052, 599042, 599038, 599034, ...
        599030, 599026, 599024, 598978, 573871, 573870, 573863, ...
        573861, 573860, 573858, 572638, 549171, 549168, 547671, ...
        547667, 547666];
    % Manaully inspected through Riscan. Fallen trees have been removed from
    % field tree before. This is just to make thigs transparent...
    %% Trees that werent segmented
    lukupuu = readtable('C:/Users/laurliik/OneDrive - University of Eastern Finland/Faro_orbis_results/Lukupuu_maastolaserkohteet.csv');

    fallen_and_clipped_trees = double([8398022 8398039 8349829 8358929 8358928 8349921 8358564 8356483 7500205 8355038 8356953 6519974 8354944 6451494 8357385 6451498 8398061 8398019 8398019 8393844 8360359 8360358 8360347 8357384 8350847 8349886 8343053 6520039 8351484 8398040 8450544 8358603 8355665 8355664 8355660 8355659 8355654 8355653 8352732 7500127 8398062 8393798 8393787 8358603 8357358 8357356 8357142 8356254 8355656 8353547 8352732 8352804 8352811 8351002 7500127 7498823 6520020 6520019 6520018 6520017 6451474 6450759 6450519 6450506 6450550 6519965 8354891 8358532 8358533 8358525 8356224 8355655 8353547 8352812 8351035 8351028 8350926 8350404 6520032 6519027 6450504 6450544]);

    fallen_trees = double([8398022 8398039 8349829 8358929 8358928 8349921 8358564 8356483 7500205 8355038 ...
        8356953 6519974 8354944 6451494 8357385 6451498 8398061 8398019 8398019 8393844 ...
        8360359 8360358 8360347 8357384 8350847 8349886 8343053 6520039 8351484 8398040 ...
        8354891 8356224 6519027
        ]);

    clipped_trees = double([
        8358532 8358533 8358525 8355655 8353547 8352812 8351035 8351028 8350926 ...
        8350404 6520032 6450504 6450544 8450544 8358603 8355665 8355664 8355660 ...
        8355659 8355654 8355653 8352732 7500127 8398062 8393798 8393787 8358603 8357358 ...
        8357356 8357142 8356254 8355656 8353547 8352732 8352804 8352811 8351002 7500127 ...
        7498823 6520020 6520019 6520018 6520017 6451474 6450759 6450519 6450506 6450550 ...
        6519965
        ]);


    clipped_trees_table =lukupuu(ismember(lukupuu.PuuId, clipped_trees), :);
    fallen_trees_table =lukupuu(ismember(lukupuu.PuuId, fallen_trees), :);


    %% Modify reference data (notFallenAndClippedTable = used in this study)

    % Remove trees that were fallen or clipped away from the shape.
    fallenAndClippedTable = lukupuu(ismember(lukupuu.PuuId, fallen_and_clipped_trees), :);
    notFallenAndClippedTable = lukupuu(~ismember(lukupuu.PuuId, fallen_and_clipped_trees), :);

    % include only trees of used plots.
    fallenAndClippedTable = fallenAndClippedTable(ismember(fallenAndClippedTable.KoealaId, include_ids), :);
    notFallenAndClippedTable = notFallenAndClippedTable(ismember(notFallenAndClippedTable.KoealaId, include_ids), :);

    %delete trees < 5cm
    fallenAndClippedTable = fallenAndClippedTable(fallenAndClippedTable.D > 5, :);
    notFallenAndClippedTable = notFallenAndClippedTable(notFallenAndClippedTable.D > 5, :);

    % Delete PuuLk 8
    fallenAndClippedTable(fallenAndClippedTable.PuuLk == 8, :) = [];
    notFallenAndClippedTable(notFallenAndClippedTable.PuuLk == 8, :) = [];

    % delete trees ulkona
    notFallenAndClippedTable(notFallenAndClippedTable.Ulkona == 1, :) = [];
    fallenAndClippedTable(fallenAndClippedTable.Ulkona == 1, :) = [];

    %%

    %%
    Fieldtree_all = notFallenAndClippedTable;
    Fieldtree_all_unprocessed = Fieldtree_all;





    %% Bring plots (areas are needed)

    plots = readtable('C:\Users\laurliik\OneDrive - University of Eastern Finland\Faro_orbis_results\kuviot_kessu.csv');

    %% Filter plots and Ulkona = 1 class! (this is done already)


    total_trees_amount = height(Fieldtree_all);%for detection rate
    plots = plots(ismember(plots.KoealaId, include_ids),:);



   matlab_treelist_combined_all_total = readtable('C:\Users\laurliik\Riegl Scans\3Dfin_A\matlab_treelist_combined_with_threedfin_A.csv');

    matlab_treelist_combined_all_total.threedfin_dbh = matlab_treelist_combined_all_total.threedfin_dbh * 100;

    matlab_treelist_combined_all_total = matlab_treelist_combined_all_total(matlab_treelist_combined_all_total.distance_to_threedfin_cm < 5,:);


    % matlab_treelist_combined_all_total = combined with 3dfin
    matlab_treelist_combined_all_total_unfilter = matlab_treelist_combined_all_total;

    % Filter based on distance

    valid_idx = ~isnan(matlab_treelist_combined_all_total.threedfin_dbh) & ...
            matlab_treelist_combined_all_total.threedfin_dbh ~= 0;

    matlab_treelist_combined_all_total = matlab_treelist_combined_all_total(valid_idx, :);
    matlab_treelist_combined_all_total_measured = matlab_treelist_combined_all_total;

    %% Logaritmic model to filter bad h/dbh pairs.

    treeTable = matlab_treelist_combined_all_total;

    zeroDBH = treeTable.threedfin_dbh == 0;
    treeTable.threedfin_dbh(zeroDBH) = NaN;

    threedfin_dbh_cm = treeTable.threedfin_dbh ;
    threedfin_hmax_cm = treeTable.threedfin_hmax; % Convert height from meters to cm
    threedfin_hmax_cm(threedfin_hmax_cm == 0) = NaN;  % Avoid division by zero
    %
    % Fit a logarithmic model
    validRows = (threedfin_dbh_cm > 0) & (threedfin_hmax_cm > 0);
    dbh_valid = threedfin_dbh_cm(validRows);
    hmax_valid = threedfin_hmax_cm(validRows);

    % Transform to logarithms
    log_DBH = log(dbh_valid);
    log_Hmax = log(hmax_valid);

    % Fit linear regression on log-transformed data
    mdl_log = fitlm(log_DBH, log_Hmax);

    % Predict heights
    hmax_predicted = exp(predict(mdl_log, log(threedfin_dbh_cm)));

    % Compute residuals and filter out extreme cases
    residuals = threedfin_hmax_cm - hmax_predicted;
    isOutlier = abs(residuals) > 2 * std(residuals, 'omitnan'); % 2-sigma rule



    outliertrees_matlab = matlab_treelist_combined_all_total(isOutlier, :); %outliertrees as own table

    outliertrees_matlab_table.(sprintf('walk_%s', walk)) = outliertrees_matlab;


    treeTableClean = treeTable(~isOutlier, :);
    matlab_treelist_combined_all_total = treeTableClean;



    %%
    

  matlab_treelist_combined_all_total.threedfin_dbh = matlab_treelist_combined_all_total.dbh;
matlab_treelist_combined_all_total.threedfin_dbh = matlab_treelist_combined_all_total.TLS_dbh*100;

    %% calculate tree finding ratios // Dr // Detection rate against reference

        
threedfin_found_trees = sum(~isnan(double(string(matlab_treelist_combined_all_total_unfilter.threedfin_match)))); % trees which have matchID

tree_segmented_ratio = (threedfin_found_trees/ total_trees_amount )*100; % Ratio of matched trees

threedfin_measured_trees = sum(~isnan(matlab_treelist_combined_all_total_unfilter.threedfin_dbh) & ...
                               matlab_treelist_combined_all_total_unfilter.threedfin_dbh ~= 0);

threedfin_measured_ratio =  (threedfin_measured_trees/ total_trees_amount )*100;; % Trees found in MATLAB


threedfin_total_trees_amount = height(matlab_treelist_combined_all_total);
    tree_finding_ratio = (threedfin_total_trees_amount / total_trees_amount) * 100; %(DR)

    % Display the result
    fprintf('%s Tree Finding Ratio: %.2f%%\n', walk, tree_finding_ratio);
    fprintf('%s Tree Segmentation Ratio: %.2f%%\n', walk, tree_segmented_ratio);
    fprintf('%s Tree Measured Ratio: %.2f%%\n', walk, threedfin_measured_ratio);




    %% This RMSE not needed
    

    % %% Calculate RMSE and Bias for Height for whole population
    % % Filter out rows with NaNs for Height. EVERY THIRD FIELD TREE HAS HEIGHT
    % % MEASUREMENT!
    % validH = ~isnan(matlab_treelist_combined_all_total.field_H);
    % 
    % validfieldH  = matlab_treelist_combined_all_total.field_H(validH);
    % validmatlabH = matlab_treelist_combined_all_total.threedfin_hmax(validH);
    % 
    % RMSE_h  = sqrt(mean((validfieldH - validmatlabH).^2));
    % Bias_h  = mean(validmatlabH - validfieldH);
    % mean_field_H = mean(validfieldH);
    % RMSE_pros = (RMSE_h / mean_field_H) * 100;
    % 
    % 
    % fprintf('Height RMSE = %.3f, Height Bias = %.3f\n', RMSE_h, Bias_h);
    % 


%%

    %% DBH: RMSE and Bias for DBH at thresholds >= 5, 7, and 10 cm for whole population

    thresholds = [5, 7, 10];

    nThresh = numel(thresholds);
    DBH_results = struct('Threshold', num2cell(thresholds), 'RMSE', [], 'Bias', [], 'SampleSize', []);
    H_results   = struct('Threshold', num2cell(thresholds), 'RMSE', [], 'Bias', [], 'SampleSize', []);

    % Loop over thresholds for DBH
    for i = 1:nThresh
        t = thresholds(i);
        mask = matlab_treelist_combined_all_total.field_D >= t & ...
            matlab_treelist_combined_all_total.threedfin_dbh > 0 & ...
            matlab_treelist_combined_all_total.field_D > 0 & ...
            ~isnan(matlab_treelist_combined_all_total.field_D) & ...
            ~isnan(matlab_treelist_combined_all_total.threedfin_dbh);
        fieldDBH  = matlab_treelist_combined_all_total.field_D(mask);
        threedfin_DBH = matlab_treelist_combined_all_total.threedfin_dbh(mask);

        RMSE_dbh = sqrt(mean((fieldDBH - threedfin_DBH).^2));
        Bias_dbh = mean(threedfin_DBH - fieldDBH);
        mean_field_D = mean(fieldDBH);
        RMSE_dbh_percent_matlab = (RMSE_dbh / mean_field_D) * 100;


        DBH_results(i).RMSE_dbh_percent_matlab = RMSE_dbh_percent_matlab;
        DBH_results(i).RMSE = RMSE_dbh;
        DBH_results(i).Bias = Bias_dbh;
        DBH_results(i).SampleSize = numel(fieldDBH);


        fprintf('DBH stats: Path: %s: DBH >= %d cm: RMSE = %.3f, Bias = %.3f, SampleSize = %d\n', ...
            walk, t, RMSE_dbh, Bias_dbh, numel(fieldDBH));
    end

    % Loop over thresholds for Height
    for i = 1:nThresh
        t = thresholds(i);
        mask = matlab_treelist_combined_all_total.field_D >= t & ...
            matlab_treelist_combined_all_total.threedfin_dbh > 0 & ...
            matlab_treelist_combined_all_total.field_H > 0 & ...
            ~isnan(matlab_treelist_combined_all_total.field_H) & ...
            ~isnan(matlab_treelist_combined_all_total.threedfin_hmax);
        fieldH  = matlab_treelist_combined_all_total.field_H(mask);
        threedfinH = matlab_treelist_combined_all_total.threedfin_hmax(mask);

        RMSE_h = sqrt(mean((fieldH - threedfinH).^2));
        Bias_h = mean(threedfinH - fieldH);
        mean_field_H = mean(fieldH);
        RMSE_h_percent_matlab = (RMSE_h / mean_field_H) * 100;


        H_results(i).RMSE_h_percent_matlab = RMSE_h_percent_matlab;
        H_results(i).RMSE = RMSE_h;
        H_results(i).Bias = Bias_h;
        H_results(i).SampleSize = numel(fieldH);

        fprintf('Height stats: Path: %s: DBH >= %d cm: RMSE = %.3f, Bias = %.3f, SampleSize = %d\n', ...
            walk, t, RMSE_h, Bias_h, numel(fieldH));
    end

    % Store the results in your whole_population structure using the variable 'walk'
    whole_population.(walk).DBH = DBH_results;
    whole_population.(walk).Height = H_results;

    %% calculate tree finding ratios // Dr // Detection rate against reference

        
threedfin_found_trees = sum(~isnan(double(string(matlab_treelist_combined_all_total_unfilter.threedfin_match)))); % trees which have matchID

    tree_segmented_ratio = (threedfin_found_trees/ total_trees_amount )*100; % Ratio of matched trees

threedfin_measured_trees = sum(~isnan(matlab_treelist_combined_all_total_unfilter.threedfin_dbh) & ...
                               matlab_treelist_combined_all_total_unfilter.threedfin_dbh ~= 0);

threedfin_measured_ratio =  (threedfin_measured_trees/ total_trees_amount )*100;; % Trees found in MATLAB


threedfin_total_trees_amount = height(matlab_treelist_combined_all_total);
    tree_finding_ratio = (threedfin_total_trees_amount / total_trees_amount) * 100; %(DR)

    % Display the result
    fprintf('%s Tree Finding Ratio: %.2f%%\n', walk, tree_finding_ratio);
    fprintf('%s Tree Segmentation Ratio: %.2f%%\n', walk, tree_segmented_ratio);
    fprintf('%s Tree Measured Ratio: %.2f%%\n', walk, threedfin_measured_ratio);    %% calculate tree finding ratios // Dr // Detection rate against reference

   
    %% Correlation and R^2 for whole population
    % Compute Correlation for DBH

    for i = 1:nThresh
        t = thresholds(i);
        mask = matlab_treelist_combined_all_total.field_D >= t & ...
            matlab_treelist_combined_all_total.threedfin_dbh > 0 & ...
            matlab_treelist_combined_all_total.field_D > 0 & ...
            ~isnan(matlab_treelist_combined_all_total.field_D) & ...
            ~isnan(matlab_treelist_combined_all_total.threedfin_dbh);
        fieldDBH  = matlab_treelist_combined_all_total.field_D(mask);
        threedfin_DBH = matlab_treelist_combined_all_total.threedfin_dbh(mask);

        r_dbh =  corr(threedfin_DBH, fieldDBH);
        r_dbh_2= r_dbh^2;

        SS_res = sum((fieldDBH - threedfin_DBH).^2);
        SS_tot = sum((fieldDBH - mean(fieldDBH)).^2);
        r2_dbh = 1 - (SS_res / SS_tot);




        DBH_results(i).r_dbh = r_dbh;
        DBH_results(i).r_dbh_2 = r_dbh_2;
        DBH_results(i).r2_dbh = r2_dbh;


    end

    % Loop over thresholds for Height
    for i = 1:nThresh
        t = thresholds(i);
        mask = matlab_treelist_combined_all_total.field_D >= t & ...
            matlab_treelist_combined_all_total.threedfin_dbh > 0 & ...
            matlab_treelist_combined_all_total.field_H > 0 & ...
            ~isnan(matlab_treelist_combined_all_total.field_H) & ...
            ~isnan(matlab_treelist_combined_all_total.threedfin_hmax);
        fieldH  = matlab_treelist_combined_all_total.field_H(mask);
        threedfinH = matlab_treelist_combined_all_total.threedfin_hmax(mask);

        r_h =  corr(threedfinH, fieldH);
        r_h_2= r_h^2;

        SS_res = sum((fieldH - threedfinH).^2);
        SS_tot = sum((fieldH - mean(fieldH)).^2);
        r2_h = 1 - (SS_res / SS_tot);



        H_results(i).r_h = r_h;
        H_results(i).r_h_2 = r_h_2;
        H_results(i).r2_h = r2_h;
        H_results(i).SampleSize_h = numel(fieldH);

    end

    % Store the results in your whole_population structure using the variable 'walk'
    whole_population.(walk).DBH = DBH_results;
    whole_population.(walk).Height = H_results;


    %% For plot level comparision, choose PuuLk 1-5

    tablePuuLk6to22 = matlab_treelist_combined_all_total(ismember(matlab_treelist_combined_all_total.field_PuuLk, 6:22),:);
    matlab_treelist_combined_all_total_1to5 = matlab_treelist_combined_all_total(ismember(matlab_treelist_combined_all_total.field_PuuLk, 1:5),:);
    Fieldtree_all_1to5 = Fieldtree_all(ismember(Fieldtree_all.PuuLk, 1:5),:);
    Fieldtree_all_6to22 = Fieldtree_all(ismember(Fieldtree_all.PuuLk, 6:22),:); % other tree classes - not used for now

    %% Calculate plot level summaries. Groupsummaries function works well in these cases.
    % CHECK !! PLOT LEVEL REFERENCE HAS ESTIMATES INCLUDED
    % hdom with all field_h

    plot_level_threedfin = groupsummary(...
        matlab_treelist_combined_all_total_1to5, ...
        "field_KoealaId", ...
        ["mean","std","min","max"], ...
        ["threedfin_dbh","threedfin_hmax"] ...
        );

    % Summarize the Field data
    plot_level_field = groupsummary(...
        Fieldtree_all_1to5, ...
        "KoealaId", ...
        ["mean","std","min","max"], ...
        ["D","H"] ...
        );

    % Now join them based on the group ID
    plot_level = outerjoin( ...
        plot_level_threedfin, plot_level_field, ...
        'LeftKeys', 'field_KoealaId', ...  % Key in plot_level_matlab
        'RightKeys', 'KoealaId', ...       % Key in plot_level_field
        'MergeKeys', true, ...             % Ensures matching values are merged
        'Type', 'full' ...                 % Keeps all groups, even if missing in one table
        );

    plot_level = renamevars(plot_level, 'field_KoealaId_KoealaId', 'field_KoealaId');
    plot_level = renamevars(plot_level, 'GroupCount_plot_level_threedfin', 'threedfin_found_trees');
    plot_level = renamevars(plot_level, 'GroupCount_plot_level_field', 'field_trees');


    %% RMSE IS CALCULATED FOR ONLY FOUND TREES... BIAS TOO
    % Initialize RMSE column
    plot_level.threedfin_RMSE_dbh = NaN(height(plot_level), 1);
    plot_level.threedfin_RMSE_dbh_percentage = NaN(height(plot_level), 1);
    % Get unique plots
    unique_plots = unique(matlab_treelist_combined_all_total_1to5.field_KoealaId);

    for i = 1:length(unique_plots)
        plot_id = unique_plots(i);

        % Extract values for the current plot
        idx = matlab_treelist_combined_all_total_1to5.field_KoealaId == plot_id;
        field_D_values = matlab_treelist_combined_all_total_1to5.field_D(idx);
        threedfin_dbh_values = matlab_treelist_combined_all_total_1to5.threedfin_dbh(idx);

        % Remove NaN values to avoid calculation issues
        valid_idx = ~isnan(field_D_values) & ~isnan(threedfin_dbh_values);
        field_D_values = field_D_values(valid_idx);
        threedfin_dbh_values = threedfin_dbh_values(valid_idx);

        % Compute RMSE if there are valid values
        if ~isempty(field_D_values) && ~isempty(threedfin_dbh_values)
            RMSE_value = sqrt(mean((field_D_values - threedfin_dbh_values).^2));
            plot_level.threedfin_RMSE_dbh(i) = RMSE_value;
        end
    end

    plot_level.threedfin_RMSE_dbh_percentage = (plot_level.threedfin_RMSE_dbh ./ plot_level.mean_D) * 100;

    plot_level.threedfin_Bias_dbh = NaN(height(plot_level), 1);

    % Get unique plots
    unique_plots = unique(matlab_treelist_combined_all_total_1to5.field_KoealaId);

    % Loop through each plot and calculate bias
    for i = 1:length(unique_plots)
        plot_id = unique_plots(i);

        % Extract values for the current plot
        idx = matlab_treelist_combined_all_total_1to5.field_KoealaId == plot_id;
        field_D_values = matlab_treelist_combined_all_total_1to5.field_D(idx);
        threedfin_dbh_values = matlab_treelist_combined_all_total_1to5.threedfin_dbh(idx);

        % Remove NaN values to avoid calculation issues
        valid_idx = ~isnan(field_D_values) & ~isnan(threedfin_dbh_values);
        field_D_values = field_D_values(valid_idx);
        threedfin_dbh_values = threedfin_dbh_values(valid_idx);

        % Compute bias if there are valid values
        if ~isempty(field_D_values) && ~isempty(threedfin_dbh_values)
            Bias_value = mean(threedfin_dbh_values - field_D_values);
            plot_level.threedfin_Bias_dbh(i) = Bias_value;
        end
    end

    %% Calculate RMSE of height per plot

    % Initialize RMSE column
    plot_level.threedfin_RMSE_h = NaN(height(plot_level), 1);

    % Get unique plots
    unique_plots = unique(matlab_treelist_combined_all_total_1to5.field_KoealaId);

    % Loop through each plot and calculate RMSE
    for i = 1:length(unique_plots)
        plot_id = unique_plots(i);

        % Extract values for the current plot
        idx = matlab_treelist_combined_all_total_1to5.field_KoealaId == plot_id;
        field_H_values = matlab_treelist_combined_all_total_1to5.field_H(idx); %field measured
        threedfin_h_values = matlab_treelist_combined_all_total_1to5.threedfin_hmax(idx); %matlab measured

        % Remove NaN values to avoid calculation issues
        valid_idx = ~isnan(field_H_values) & ~isnan(threedfin_h_values);
        field_H_values = field_H_values(valid_idx);
        threedfin_h_values = threedfin_h_values(valid_idx);

        % Compute RMSE if there are valid values
        if ~isempty(field_H_values) && ~isempty(threedfin_h_values)
            RMSE_value = sqrt(mean((field_H_values - threedfin_h_values).^2));
            plot_level.threedfin_RMSE_h(i) = RMSE_value;
        end
    end


    plot_level.threedfin_RMSE_h_percentage = (plot_level.threedfin_RMSE_h ./ plot_level.mean_H) * 100;

    % Initialize H Bias column
    plot_level.threedfin_Bias_h= NaN(height(plot_level), 1);

    % Get unique plots
    unique_plots = unique(matlab_treelist_combined_all_total_1to5.field_KoealaId);

    % Loop through each plot and calculate H Bias
    for i = 1:length(unique_plots)
        plot_id = unique_plots(i);

        % Extract values for the current plot
        idx = matlab_treelist_combined_all_total_1to5.field_KoealaId == plot_id;
        field_H_values = matlab_treelist_combined_all_total_1to5.field_H(idx);  % field measured height
        threedfin_h_values = matlab_treelist_combined_all_total_1to5.threedfin_hmax(idx);            % MATLAB measured height

        % Remove NaN values to avoid calculation issues
        valid_idx = ~isnan(field_H_values) & ~isnan(threedfin_h_values);
        field_H_values = field_H_values(valid_idx);
        threedfin_h_values = threedfin_h_values(valid_idx);

        % Compute bias if there are valid values
        if ~isempty(field_H_values) && ~isempty(threedfin_h_values)
            % Compute bias as the mean difference: field_H - hmax
            Bias_value = mean(threedfin_h_values - field_H_values);
            plot_level.threedfin_Bias_h(i) = Bias_value;
        end
    end



    %% Create 'areas' table with proper variable names
    areas = table(plots.KoealaId, plots.Pala, 'VariableNames', {'PlotId', 'area'});

    % Rename 'PlotId' in areas to 'field_KoealaId'
    areas = renamevars(areas, 'PlotId', 'field_KoealaId');

    % Now both tables have the same key name: 'field_KoealaId'
    plot_level = outerjoin(plot_level, areas, ...
        'Keys', 'field_KoealaId', ...
        'Type', 'left', ...
        'MergeKeys', true);

    %% stems per hectare matlab

    plot_level.threedfin_stems_per_ha = (plot_level.threedfin_found_trees./plot_level.area)*10000 ;




    %% Calculate Basal Area for matlab
    % Identify rows where ba == 0, because some trees were processed through
    % fallback
   
    matlab_treelist_combined_all_total_1to5.threedfin_ba = pi * (matlab_treelist_combined_all_total_1to5.threedfin_dbh / 200).^2;



    threedfin_ba_sum = groupsummary(matlab_treelist_combined_all_total_1to5, ...
        "field_KoealaId", "sum", "threedfin_ba");

    threedfin_ba_sum = renamevars(threedfin_ba_sum, "sum_threedfin_ba", "threedfin_Total_BA");

    plot_level = outerjoin(plot_level, threedfin_ba_sum(:, ["field_KoealaId","threedfin_Total_BA"]), ...
        "Keys",       "field_KoealaId", ...
        "Type",       "left", ...
        "MergeKeys",  true);

    plot_level.threedfin_ba_per_ha = (plot_level.threedfin_Total_BA./plot_level.area)*10000;

    %% Calculate Basal Area for ALL field AND STEMS PER HECTARE

    %Stems per hectare
    plot_level.field_stems_per_ha = (plot_level.field_trees./plot_level.area)*10000 ;


    Fieldtree_all_1to5.field_ba = pi * (Fieldtree_all_1to5.D/200).^2;

    ba_fieldsum = groupsummary(Fieldtree_all_1to5, ...
        "KoealaId", "sum", "field_ba");

    ba_fieldsum = renamevars(ba_fieldsum, "sum_field_ba", "field_Total_BA");
    ba_fieldsum = renamevars(ba_fieldsum, "KoealaId", "field_KoealaId");
    ba_fieldsum = renamevars(ba_fieldsum, "GroupCount", "fieldstems");



    plot_level = outerjoin(plot_level, ba_fieldsum(:, ["field_KoealaId","field_Total_BA"]), ...
        "Keys",       "field_KoealaId", ...
        "Type",       "left", ...
        "MergeKeys",  true);

    plot_level.field_ba_per_ha = (plot_level.field_Total_BA./plot_level.area)*10000;

    %calculate stems per hectare from group count that has been used in ba per
    %hectare calculations


    %% Volume per hectare (this is not good, model gives better results
    % 
    % vol_sum = groupsummary(matlab_treelist_combined_all_total_1to5, ...
    %     "field_KoealaId", "sum", "vol");
    % 
    % vol_sum = renamevars(vol_sum, "sum_vol", "Total_vol");
    % 
    % plot_level = outerjoin(plot_level, vol_sum(:, ["field_KoealaId","Total_vol"]), ...
    %     "Keys",       "field_KoealaId", ...
    %     "Type",       "left", ...
    %     "MergeKeys",  true);
    % 
    % plot_level.vol_per_ha = (plot_level.Total_vol./plot_level.area)*10000;

    %% volume for field trees is just from plot level list


    field_vol = table(plots.KoealaId, plots.V, 'VariableNames', {'field_KoealaId', 'V'});
    field_vol = renamevars(field_vol, "V", "field_vol");

    % Now both tables have the same key name: 'field_KoealaId'
    plot_level = outerjoin(plot_level, field_vol, ...
        'Keys', 'field_KoealaId', ...
        'Type', 'left', ...
        'MergeKeys', true);

    %% bias_corrected RMSE = sqrt( RMSE^2 - BIAS^2 )



    %% PPA:lla painotettu lpm ja pituus
matlab_treelist_combined_all_total_1to5.threedfin_ba = pi * (matlab_treelist_combined_all_total_1to5.threedfin_dbh / 200).^2;


    matlab_treelist_combined_all_total_1to5.dbhTimesBA = ...
        matlab_treelist_combined_all_total_1to5.threedfin_dbh .* ...
        matlab_treelist_combined_all_total_1to5.threedfin_ba;

    plotSums_dbh = groupsummary(matlab_treelist_combined_all_total_1to5, ...
        "field_KoealaId", ...
        "sum", ...
        ["dbhTimesBA", "threedfin_ba"]);

    plotSums_dbh.threedfin_Dg = plotSums_dbh.sum_dbhTimesBA ./ plotSums_dbh.sum_threedfin_ba;


    plot_level = outerjoin(plot_level, plotSums_dbh(:, {'field_KoealaId', 'threedfin_Dg'}), ...
        'Keys', 'field_KoealaId', ...
        'Type', 'left', ...
        'MergeKeys', true);


    %% for Hg

    matlab_treelist_combined_all_total_1to5.hTimesBA = matlab_treelist_combined_all_total_1to5.threedfin_hmax .* matlab_treelist_combined_all_total_1to5.threedfin_ba;

    plotSums_H = groupsummary(matlab_treelist_combined_all_total_1to5, ...
        "field_KoealaId", ...
        "sum", ...
        ["hTimesBA", "threedfin_ba"]);

    plotSums_H.threedfin_Hg = plotSums_H.sum_hTimesBA ./ plotSums_H.sum_threedfin_ba;


    plot_level = outerjoin(plot_level, plotSums_H(:, {'field_KoealaId', 'threedfin_Hg'}), ...
        'Keys', 'field_KoealaId', ...
        'Type', 'left', ...
        'MergeKeys', true);



    %% HG and Dg for field
    Fieldtree_all_1to5.field_dbhTimesBA = Fieldtree_all_1to5.D .* Fieldtree_all_1to5.field_ba;

    plotSums_field_dbh = groupsummary(Fieldtree_all_1to5, "KoealaId", "sum", ...
        ["field_dbhTimesBA", "field_ba"]);

    plotSums_field_dbh.field_Dg = plotSums_field_dbh.sum_field_dbhTimesBA ./ plotSums_field_dbh.sum_field_ba;

    plotSums_field_dbh = renamevars(plotSums_field_dbh, 'KoealaId', 'field_KoealaId');

    plot_level = outerjoin(plot_level, plotSums_field_dbh(:, {'field_KoealaId', 'field_Dg'}), ...
        'Keys', 'field_KoealaId', ...
        'Type', 'left', ...
        'MergeKeys', true);

    % Field_Hg (remove NAN or 0 heights). OR USE HEst?

    valid_rows = ~isnan(Fieldtree_all_1to5.H) & Fieldtree_all_1to5.H > 0;
    Fieldtree_all_1to5_h = Fieldtree_all_1to5(valid_rows, :);

    Fieldtree_all_1to5_h.field_hTimesBA = Fieldtree_all_1to5_h.H .* Fieldtree_all_1to5_h.field_ba;

    plotSums_field_h = groupsummary(Fieldtree_all_1to5_h, ...
        "KoealaId", ...
        "sum", ...
        ["field_hTimesBA", "field_ba"]);

    plotSums_field_h.field_Hg = plotSums_field_h.sum_field_hTimesBA ./ plotSums_field_h.sum_field_ba;

    plotSums_field_h = renamevars(plotSums_field_h, 'KoealaId', 'field_KoealaId');

    plot_level = outerjoin(plot_level, plotSums_field_h(:, {'field_KoealaId', 'field_Hg'}), ...
        'Keys', 'field_KoealaId', ...
        'Type', 'left', ...
        'MergeKeys', true);

   % clear plotSums_field_h n:tresh plot_id mean_field_H 
    %% Recognition percentage  // Detection rate(Dr)

    plot_level.finding_ratio = plot_level.threedfin_found_trees ./ plot_level.field_trees;
    %% Volume with Laasasenaho

    % Puulajikoodit:
    % 01 Mänty, 12 Kontortamänty, 23 Serbiankuusi
    % 02 Kuusi, 13 Kynäjalava, 24 Tammi
    % 03 Rauduskoivu, 14 Lehtikuusi, 25 Tuomi
    % 04 Hieskoivu ,15 Metsälehmus, 26 Vaahtera
    % 05 Haapa, 16 Mustakuusi 27, Visakoivu
    % 06 Harmaaleppä, 17 Paju, 28 Vuorijalava
    % 07 Tervaleppä, 18 Pihlaja, 29 Lehtipuu
    % 08 Muu havupuu, 19 Pihta, 30 Havupuu
    % 09 Muu lehtipuu, 20 Raita, 31 GPS-piste
    % 10 Douglaskuusi, 21 Saarni
    % 11 Kataja, 22 Sembramänty
    % For safety, create or initialize a new column for volume:
    matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen = NaN(height(matlab_treelist_combined_all_total_1to5), 1);

    %% 1) MÄNTY

    % Logical index for rows with species = 'mänty'
    Manty = (matlab_treelist_combined_all_total_1to5.treesp == 1)& (matlab_treelist_combined_all_total_1to5.threedfin_hmax > 1.3);
    % Apply the pine volume formula to those rows
    matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen(Manty) = ...
        0.036089 .* (matlab_treelist_combined_all_total_1to5.threedfin_dbh(Manty)).^2.01395 .* ...
        (0.99676 .^ matlab_treelist_combined_all_total_1to5.threedfin_dbh(Manty)) .* ...
        (matlab_treelist_combined_all_total_1to5.threedfin_hmax(Manty)).^2.07025 .* ...
        (matlab_treelist_combined_all_total_1to5.threedfin_hmax(Manty) - 1.3).^-1.07209;

    %% 2) KUUSI
    Kuusi =  (matlab_treelist_combined_all_total_1to5.treesp == 2)& (matlab_treelist_combined_all_total_1to5.threedfin_hmax > 1.3);;
    matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen(Kuusi) = ...
        0.022927 .* (matlab_treelist_combined_all_total_1to5.threedfin_dbh(Kuusi)).^1.91505 .* ...
        (0.99146 .^ matlab_treelist_combined_all_total_1to5.threedfin_dbh(Kuusi)) .* ...
        (matlab_treelist_combined_all_total_1to5.threedfin_hmax(Kuusi)).^2.82541 .* ...
        (matlab_treelist_combined_all_total_1to5.threedfin_hmax(Kuusi) - 1.3).^-1.53547;

    %% 3) KOIVU ja muut
    Koivu = (matlab_treelist_combined_all_total_1to5.treesp >= 3) & ...
        (matlab_treelist_combined_all_total_1to5.treesp <= 22) & (matlab_treelist_combined_all_total_1to5.threedfin_hmax > 1.3);
    matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen(Koivu) = ...
        0.011197 .* (matlab_treelist_combined_all_total_1to5.threedfin_dbh(Koivu)).^2.10253 .* ...
        (0.98600 .^ matlab_treelist_combined_all_total_1to5.threedfin_dbh(Koivu)) .* ...
        (matlab_treelist_combined_all_total_1to5.threedfin_hmax(Koivu)).^3.98519 .* ...
        (matlab_treelist_combined_all_total_1to5.threedfin_hmax(Koivu) - 1.3).^-2.65900;



    %% IF under h<1,3 then inf -> give 0
    small_trees = matlab_treelist_combined_all_total_1to5.threedfin_hmax <= 1.3;

    % Set volume to 0 for these small trees
    matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen(small_trees) = 0;

    % Convert Laasasen volume from dm³ to m³
    matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen_m3 = ...
        matlab_treelist_combined_all_total_1to5.threedfin_volumeLaasasen / 1000;

    %% Total volume Laasasenaho
    vol_sum_Laasasen = groupsummary(matlab_treelist_combined_all_total_1to5, ...
        "field_KoealaId", "sum", "threedfin_volumeLaasasen_m3");

    vol_sum_Laasasenaho = renamevars(vol_sum_Laasasen, "sum_threedfin_volumeLaasasen_m3", "threedfin_total_volumeLaasasen");

    plot_level = outerjoin(plot_level, vol_sum_Laasasenaho(:, ["field_KoealaId","threedfin_total_volumeLaasasen"]), ...
        "Keys",       "field_KoealaId", ...
        "Type",       "left", ...
        "MergeKeys",  true);


    plot_level.threedfin_Laasasen_vol_per_ha = (plot_level.threedfin_total_volumeLaasasen./plot_level.area)*10000;




    %% hdom field all


    plot_level.hdom_field_all = NaN(height(plot_level),1);     % empty column

    for k = 1:height(plot_level)

        pid      = plot_level.field_KoealaId(k);      % current plot ID
        area_m2  = plot_level.area(k);                % area of that plot

        % How many trees correspond to 100 stems · ha⁻¹ on THIS plot?
        n_dom = max(1, round( area_m2 * 100 / 10000 ));  % always at least 1

        % Extract tree records of this plot
        rows = Fieldtree_all_1to5.KoealaId == pid;
        tr   = Fieldtree_all_1to5(rows, :);

        % Protect against plots with fewer trees than n_dom
        if isempty(tr);  continue;  end
        n_dom = min(n_dom, height(tr));

        % Sort by DBH (descending) and take the first n_dom
        tr    = sortrows(tr, 'D', 'descend');
        sel   = tr.HEst(1:n_dom);                         % heights of dominant trees

        % Dominant height for the plot
        plot_level.hdom_field_all(k) = mean(sel, 'omitnan');
    end
    %% hdom MLS


    plot_level.hdom_threedfin = NaN(height(plot_level),1);     % empty column

    for k = 1:height(plot_level)

        pid      = plot_level.field_KoealaId(k);      % current plot ID
        area_m2  = plot_level.area(k);                % area of that plot

        % How many trees correspond to 100 stems · ha⁻¹ on THIS plot?
        n_dom = max(1, round( area_m2 * 100 / 10000 ));  % always at least 1

        % Extract tree records of this plot
        rows = matlab_treelist_combined_all_total_1to5.field_KoealaId == pid;
        tr   = matlab_treelist_combined_all_total_1to5(rows, :);

        % Protect against plots with fewer trees than n_dom
        if isempty(tr);  continue;  end
        n_dom = min(n_dom, height(tr));

        % Sort by DBH (descending) and take the first n_dom
        tr_2    = sortrows(tr, 'threedfin_dbh', 'descend');
        sel_2   = tr_2.threedfin_hmax(1:n_dom);                         % heights of dominant trees

        % Dominant height for the plot
        plot_level.hdom_threedfin(k) = mean(sel_2, 'omitnan');
    end

    %% Plot level Residuals

    plot_level.Dg_residual_1 = plot_level.threedfin_Dg-plot_level.field_Dg;
    plot_level.Hg_residual_1 = plot_level.threedfin_Hg- plot_level.field_Hg ;
    plot_level.BA_residual_1 = plot_level.threedfin_ba_per_ha-plot_level.field_ba_per_ha;
    plot_level.TPH_residual_1 = plot_level.threedfin_stems_per_ha-plot_level.field_stems_per_ha;
    plot_level.Hdom_residual_1 = plot_level.hdom_threedfin-plot_level.hdom_field_all;
    plot_level.vol_residual_1 = plot_level.threedfin_Laasasen_vol_per_ha-plot_level.field_vol;

    %% Plot level Residuals (this is wrong way around)
    % 
    % plot_level.Dg_residual = plot_level.field_Dg-plot_level.Dg;
    % plot_level.Hg_residual = plot_level.field_Hg-plot_level.Hg;
    % plot_level.BA_residual = plot_level.field_ba_per_ha-plot_level.threedfin_ba_per_ha;
    % plot_level.TPH_residual = plot_level.field_stems_per_ha-plot_level.threedfin_stems_per_ha;
    % 






    %% Correlations for all trees (Tree classes 1-5 NOT filtered)

    r_dbh = corr(matlab_treelist_combined_all_total.threedfin_dbh, matlab_treelist_combined_all_total.field_D);

    valid_idx = ~isnan(matlab_treelist_combined_all_total.threedfin_hmax) & ~isnan(matlab_treelist_combined_all_total.field_H);
    r_h = corr(matlab_treelist_combined_all_total.threedfin_hmax(valid_idx), matlab_treelist_combined_all_total.field_H(valid_idx));

    %correlation for plot level
    r_Hg = corr(plot_level.threedfin_Hg, plot_level.field_Hg);

    r_Dg = corr(plot_level.threedfin_Dg, plot_level.field_Dg);

    r_TPH = corr(plot_level.threedfin_stems_per_ha, plot_level.field_stems_per_ha);

    r_BA = corr(plot_level.threedfin_ba_per_ha, plot_level.field_ba_per_ha);

    r_Vol = corr(plot_level.threedfin_Laasasen_vol_per_ha, plot_level.field_vol);
    % Create a table with variable names and correlation values
    correlation_table = table(...
        ["DBH vs Field D";
        "Hmax vs Field H";
        "Hg vs Field Hg";
        "Dg vs Field Dg";
        "Stems per ha vs Field Stems per ha";
        "BA per ha vs Field BA per ha";
        "Volume per ha vs Field Volume"], ...
        [r_dbh; r_h; r_Hg; r_Dg; r_TPH; r_BA; r_Vol], ...
        'VariableNames', {'Variable_Pair', 'Pearson_Correlation'});
    correlation_tables.(walk) = correlation_table;% Display the table in a readable format



    %% Save and structurize plot_level table
    % 
    % outputFilePath = sprintf('C:/Users/laurliik/OneDrive - University of Eastern Finland/Faro_orbis_results/results_3/%s_plot_level_results.csv', walk);
    % writetable(plot_level, outputFilePath);
    plot_levels.(walk) = plot_level;
    % 

    %% ALL PLOTS RMSE plot level
    % this is needed for plot level comparision


    stats_plot_level.RMSE_Dg = sqrt(mean((plot_level.threedfin_Dg - plot_level.field_Dg).^2));
    stats_plot_level.RMSE_Hg = sqrt(mean((plot_level.threedfin_Hg - plot_level.field_Hg).^2));
    stats_plot_level.RMSE_BA = sqrt(mean((plot_level.threedfin_ba_per_ha - plot_level.field_ba_per_ha).^2));
    stats_plot_level.RMSE_TPH = sqrt(mean((plot_level.threedfin_stems_per_ha - plot_level.field_stems_per_ha).^2));
    stats_plot_level.RMSE_vol = sqrt(mean((plot_level.threedfin_Laasasen_vol_per_ha - plot_level.field_vol).^2));
    stats_plot_level.RMSE_hdom = sqrt(mean((plot_level.hdom_field_all - plot_level.hdom_threedfin).^2));


    stats_plot_level.bias_Dg = mean(plot_level.threedfin_Dg - plot_level.field_Dg);
    stats_plot_level.bias_Hg = mean(plot_level.threedfin_Hg - plot_level.field_Hg);
    stats_plot_level.bias_BA = mean(plot_level.threedfin_ba_per_ha - plot_level.field_ba_per_ha);
    stats_plot_level.bias_TPH = mean(plot_level.threedfin_stems_per_ha - plot_level.field_stems_per_ha);
    stats_plot_level.bias_vol = mean(plot_level.threedfin_Laasasen_vol_per_ha - plot_level.field_vol);
    stats_plot_level.bias_hdom = mean(plot_level.hdom_threedfin - plot_level.hdom_field_all);


    stats_plot_level.r_Hg = corr(plot_level.threedfin_Hg, plot_level.field_Hg);
    stats_plot_level.r_Dg = corr(plot_level.threedfin_Dg, plot_level.field_Dg);
    stats_plot_level.r_TPH = corr(plot_level.threedfin_stems_per_ha, plot_level.field_stems_per_ha);
    stats_plot_level.r_BA = corr(plot_level.threedfin_ba_per_ha, plot_level.field_ba_per_ha);
    stats_plot_level.r_vol = corr(plot_level.threedfin_Laasasen_vol_per_ha, plot_level.field_vol);
    stats_plot_level.r_hdom = corr(plot_level.hdom_threedfin, plot_level.hdom_field_all);

    stats_plot_level.RMSE_Dg_percentage = stats_plot_level.RMSE_Dg /mean(plot_level.field_Dg);
    stats_plot_level.RMSE_Hg_percentage = stats_plot_level.RMSE_Hg /mean(plot_level.field_Hg);
    stats_plot_level.RMSE_BA_percentage = stats_plot_level.RMSE_BA /mean(plot_level.field_ba_per_ha);
    stats_plot_level.RMSE_TPH_percentage = stats_plot_level.RMSE_TPH /mean(plot_level.field_stems_per_ha);
    stats_plot_level.RMSE_vol_percentage = stats_plot_level.RMSE_vol /mean(plot_level.field_vol);
    stats_plot_level.RMSE_hdom_percentage = stats_plot_level.RMSE_hdom /mean(plot_level.hdom_field_all);

    stats_plot_levels.(walk) = stats_plot_level;



    %% structurize all measured trees (NOT filetered 1-5 classes) for tree level analysis)
    All_measured_trees.(walk) = matlab_treelist_combined_all_total;


    end

clearvars -except All_measured_trees plot_levels whole_population stats_plot_levels Fieldtree_all_unprocessed



