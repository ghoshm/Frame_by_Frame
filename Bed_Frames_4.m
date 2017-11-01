% Bed_Frames_4 

% V2 - aims to correct errors as it loads each Excel sheet 
    % This is far easier than trying to handle the errors once the data is 
    % Combined as the combined set is so large that logically indexing is 
    % difficult/impossible 
    
% V3 - Ease of combining experiments 
    % Aims to apply logical indicies to each bout (rather than relying on
        % time boundaries) 
    % Aims to seperate analysis from figures and stats 

% V4 - incorporates Jason's corrections
    % Log axis for distribution figures 
    
%% Hard Coded Variables 
lines_per_sheet = 50000; % Specify the number of data points per Excel Sheet 
%% Required Scripts 

% dir2 - marcus.ghosh.11@ucl.ac.uk 

% Nat Sort Files -
    %http://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort 
    
%ProgressBar 
   %http://uk.mathworks.com/matlabcentral/fileexchange/6922-progressbar
   
%% Load Data from Excel Sheets 
tic
folder_path = uigetdir; % Choose your experiment
folder_open = dir2(folder_path); % Open this folder
disp(horzcat('Running File ',folder_path)); % Report file choice  

% Pre-allocation
time = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % time {1}
pause(30); % Wait for memory 
data_type = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Data type {2} 
pause(30); % Wait for memory 
fish_id = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Fish id {3}
pause(30); % Wait for memory 
delta_px = nan(size(folder_open,1)*lines_per_sheet,1,'single'); % Delta px {4}
pause(30); % Wait for memory 
sheet_names = cell(size(folder_open,1),1); % Excel sheet names  

% Ordering by File Name 
for f = 1:size(folder_open,1) % For each excel file 
    sheet_names{f} = folder_open(f).name; % Take it's name 
end % Note that these will be in "computer" order 
    % Ie. 1-10, 100, 1000 etc 

[~,O] = natsortfiles(sheet_names); % Use natsortfiles to sort by file name
    clear sheet_names; % Clear Sheet names 
    
a = 1; % Start a counter 
progress = 0; % Start a timer
data_type_errors = 0; % Start a counter  
order_errors = 0; % Start a counter 
order_errors_size = []; % pre-allocate an empty vector
progressbar('Files') %Initialise progress bars 
for f = O' % For each Excel file
    
    % Load data 
        % Raw data structure 
            % 1 - Time 
            % 2 - Data type 
            % 3 - Fish ID 
            % 4 - Delta px
    fid = fopen(strcat(folder_path,'\',folder_open(f).name)); % Open it
    if f == O(1) % For the first file Skip the first blank lines
        raw_text = textscan(fid,... % Read in the data
            '%*f32 %f32 %*f32 %f32 %*1s %s %f32','headerlines',192 + 1);
    else % For the rest of the files 
        raw_text = textscan(fid,... % Read in the data
            '%*f32 %f32 %*f32 %f32 %*1s %s %f32','headerlines',1);
    end
    
    % Error Handling 
    
    % Dealing with data_type errors (e)
        % Note that this includes cropping the end of the experiment off 
    found = find(raw_text{2} ~= 101); % Find errors
    if isempty(found) ~= 1 % If there are errors 
        for e = 1:size(raw_text,2) % For each data type being imported
            raw_text{e}(found) = []; % Remove errors
        end
        data_type_errors = data_type_errors + 1; % Add to error counter  
    end
    clear found;
    
    % Dealing with ordering errors 
    fish_id_order = str2num(char(raw_text{3})); % Convert str to num
    fish_id_order_check = diff(fish_id_order); % Diff these values
    found = find(fish_id_order_check ~= 1 & fish_id_order_check ~= -191 &...
        isnan(fish_id_order_check) == 0,1,'first'); % Find the first frame drops 
    
    while isempty(found) == 0 % While there are dropped frames 
        
        % Use a sliding window to find where the pattern normalises
        % Cut backwards
        found_clip_b = [191 1];
        while sum(fish_id_order_check(found - found_clip_b(1):found - found_clip_b(2))) ~= 191
            found_clip_b = found_clip_b + 1; % Slide window
            
            % Catch Running past the start of the file exception
            if found - found_clip_b(1) <= 1
                found_clip_b = [found_clip_b(1) + 2 found_clip_b(1) + 2];
                break
            end 
        end
        
        % Cut forwards
        found_clip_f = [1 191];
        while sum(fish_id_order_check(found + found_clip_f(1):found + found_clip_f(2))) ~= 191
            found_clip_f = found_clip_f + 1; % Slide window 
            
            % Catch Running past the end of the file exception 
             if found + found_clip_f(2) >= size(fish_id_order_check,1)
                 found_clip_f = [found_clip_f(2)+2 found_clip_f(2)+2]; 
                 break 
             end 
        end
        
        % Now set values between these sections to NaN
        fish_id_order(found - found_clip_b(2)+2:found + found_clip_f(1)-1) = NaN;
        order_errors_size(order_errors+1,1) = ...
            size(found - found_clip_b(2)+2:found + found_clip_f(1)-1,2); 
            % Store the size of the removed data 
        clear found_clip_b found_clip_f
        
        fish_id_order_check = diff(fish_id_order); % Diff the new fish order
        found = find(fish_id_order_check ~= 1 & fish_id_order_check ~= -191 &...
            isnan(fish_id_order_check) == 0,1,'first'); % Find other frame drops
        order_errors = order_errors + 1; % Add to the order errors counter 
    end
        
    % Data Storage 
    if a == 1 % For the first file 
        time(a:size(raw_text{1},1),1) = raw_text{1}; % Time 
        data_type(a:size(raw_text{1},1),1) = raw_text{2}; % Data Type 
        fish_id(a:size(raw_text{1},1),1) = fish_id_order; % Fish Id order 
        delta_px(a:size(raw_text{1},1),1) = raw_text{4}; % Delta px 
            a = a + size(raw_text{1},1); % Add to counter 
    else % For all other files 
        time(a:a+size(raw_text{1},1)-1,1) = raw_text{1}; % Time 
        data_type(a:a+size(raw_text{1},1)-1,1) = raw_text{2}; % Data Type 
        fish_id(a:a+size(raw_text{1},1)-1,1) = fish_id_order; % Fish Id order 
        delta_px(a:a+size(raw_text{1},1)-1,1) = raw_text{4}; % Delta px 
            a = a + size(raw_text{1},1); % Add to counter
    end
        
    % Clean up 
    fclose(fid); clear fish_id_order fish_id_order_check found raw_text; % Clear variables
    
    progress = progress + 1; % Add to timer 
    progressbar(progress/size(folder_open,1)); %Update the Fish progressbar
    
end
toc 

% Report back on errors 
disp(horzcat(num2str(data_type_errors),' Files with Data Type Errors')); 
disp(horzcat(num2str(order_errors),' Order Errors')); 
for e = 1:size(order_errors_size,1) % For each order error 
    disp(horzcat('Order Error ',num2str(e),' Size = ',num2str(order_errors_size(e)))); 
end 
disp(horzcat('Ran File ',folder_path)); % Report file choice  

clear a e ans data_type data_type_errors f fid folder_open folder_path ...
    lines_per_sheet O order_errors order_errors_size progress

%% Reshape The Data 

% Check the number of frames per fish is correct 
frames_per_fish = zeros(1,max(fish_id)); 
for f = 1:max(fish_id) % For each fish 
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    frames_per_fish(f) = size(found,1); % Store number of frames 
    disp(horzcat('Found frames for fish Number ',num2str(f),' of ',...
        num2str(max(fish_id)))); % Report on progress 
end 

if min(frames_per_fish) ~= max(frames_per_fish) % If the number of frames is not equal 
   error('Data formatted Incorrectly'); % Call an error  
end 

% Delta px sq 
for f = 1:max(fish_id) % For each fish  
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    
    if f == 1 % For the first fish 
        delta_px_sq = nan(frames_per_fish(1),max(fish_id),'single'); 
        % Pre-allocate 
    end
    
    delta_px_sq(:,f) = delta_px(found); % Take delta_px values
    disp(horzcat('Collected frames for fish ',num2str(f),...
        ' of ',num2str(max(fish_id)))); % Report on progress 
end 
clear f found delta_px

% Time sq 
    % Note this is separated from delta px sq for ease of memory handling 
for f = 1:max(fish_id) % For each fish  
    clear found; 
    found = find(fish_id == f); % Find it's data points 
    
    if f == 1 % For the first fish 
        time_sq = nan(frames_per_fish(1),max(fish_id),'single'); 
        % Pre-allocate 
    end
    
    time_sq(:,f) = time(found); % Take time values 
    disp(horzcat('Collected time for fish ',num2str(f),...
        ' of ',num2str(max(fish_id)))); % Report on progress 
end
clear f found time fish_id frames_per_fish

%% Add in time 

%Load a XLS (Perl Batch Data File)
[filename, pathname] = uigetfile('*.txt',...
    'Select a Perl Batch Output (DATA) File'); % Select a geno file
if isequal(filename,0) % If no file is selected
    error('No File Selected') % Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) % Show selected filename
end

fid = fopen(strcat(pathname,filename)); % Open it
formatSpec = '%*884s%f%[^\n\r]';
dataArray = textscan(fid, formatSpec, 3-3+1, 'Delimiter',...
    '', 'WhiteSpace', '', 'HeaderLines', 2, 'ReturnOnError', false, 'EndOfLine', '\r\n');
start_time = dataArray{1}; % Extract start time (Matlab generated code) 
fclose(fid); % Close it 

clear ans dataArray fid filename formatSpec pathname

% Light Boundaries Calculation
time_sq_max = ((max(time_sq'))/(60*60))'; % Max time @ each frame (hours from start) 

a = 1; % start a counter
time_counter = 0; % Start a counter 
boundary = 14 - start_time; % Assumes the experiment starts during the day
while time_counter < time_sq_max(end) - 10 % (10 allows for at least a night)  
    lb(a,1) = knnsearch(time_sq_max,boundary); % Find the best match 
    if mod(a,2) == 1 % If odd
        boundary = boundary + 10; % Add night hours
    else
        boundary = boundary + 14; % Add day hours
    end
    time_counter = time_sq_max(lb(a,1)); % Set time counter 
    disp(horzcat('Found Light Boundary = ',num2str(a))); % Report progress
    a = a + 1; % Add to counter
end

% Day vs Night
dn = ones(size(time_sq,1),1,'single'); % Pre-allocate 
for t = 1:2:size(lb,1) % For each night boundary
    dn(lb(t):lb(t+1)-1) = 0;   
end 

clear a boundary start_time time_counter time time_sq t 

%% Normalise the data 

% First work out pixel noise 
    % This is possible only if one of the boxes is empty 
%     scrap = delta_px_sq(:,1:96); 
%     prctile(scrap(:)',[0 10 20 30 40 50 60 70 80 90])
%     prctile(scrap(:)',[0 10 20 30 40 50 60 70 80 90])

delta_px_sq(:,1:96) = []; % Remove the unused box 
delta_px_sq = delta_px_sq - 1; % Remove pixel noise 

% Now normalise for each fish 
% scrap = delta_px_sq; scrap(scrap == 0) = NaN; 
% delta_px_sq = delta_px_sq./nanmedian(scrap);
% clear scrap 

%% Group the data by condition 

 %Load a geno_list 
 [filename, pathname] = uigetfile('*.txt', 'Select a Genotype List'); % Select a geno file
 if isequal(filename,0) % If no file is selected
     error('No File Selected') % Show Error
 else %If selected
     disp(['User selected ', fullfile(pathname, filename)]) % Show selected filename
 end
 
geno_list = importdata(strcat(pathname,filename),'\t',2); % Load genolist 

% Generate group tags 
group_tags = nan(size(delta_px_sq,2),1); % Pre-allocate
for g = 1:size(geno_list.data,2) % For each group 
    group_tags(geno_list.data(1:find(isnan(geno_list.data(:,g))==0,1,'last')...
        ,g)) = g; % Assign group membership  
end 
delta_px_sq(:,isnan(group_tags)) = []; % Remove data   
group_tags(isnan(group_tags)) = []; % Remove blank values 

clear filename g pathname 

%% Extract Parameters from Data (< Re-shaped data)

% Variables 
tic
% Calculate an approximate frame rate 
fps = round(1/((nanmean(diff(time_sq_max)))*(60*60)),1); 
    % Diff time_sq_max to find time between frames 
    % Take a mean, Convert to mins then hours 
    % Divide one by this value and round -> frames per second
    
parameters = {'Active Bout Length','Active Bout Mean',...
    'Active Bout Variance','Active Bout Total',...
    'Active Bout Minimum','Active Bout Maximum','Number of Active Bouts',...
    'Total Time Active','Total Activity','Inactive Bout Length','Number of Inactive Bouts',...
    'Total Time Inactive'}; % Specify parameters  

% Specify how to convert frame values to parameter units (for figures)
    % Note that these are denominators (will later be used for division) 
unit_conversion(1,:) = [fps 1 1 1 1 1 1 (fps*3600) 1 fps 1 (fps*3600)];
unit_conversion(2,:) = [fps 1 1 1 1 1 1 1 1 fps 1 1];  

% Specify units (for figures) 
units = {'Seconds','Delta Px','Delta Px','Delta Px','Delta Px','Delta Px',...
    'No.','Hours','Delta Px','Seconds','No.','Hours'};
units_2 = {'Seconds','Delta Px','Delta Px','Delta Px','Delta Px','Delta Px',...
    'No.',horzcat('Seconds/',num2str(fps),'s'),'Delta Px','Seconds',...
    'No.',horzcat('Seconds/',num2str(fps),'s')};

% Specify Smoothing operation (0 = mean, 1 = total, 2 = max) - for figures 
parameter_smooth(1:size(parameters,2)) = 0; 
parameter_smooth(7:8) = 1; parameter_smooth(9) = 2; parameter_smooth(11:12) = 1; 

% Pre-allocate 
wake_cells = cell(1,size(delta_px_sq,2)); % Wake Cells (bout parameters)  
sleep_cells = cell(1,size(delta_px_sq,2)); % Sleep Cells (bout parameters) 

for p = 1:size(parameters,2) % For each parameter
    parameter_time{p} = nan(size(delta_px_sq),'single'); % Parameters across time 
end 

% Finding transitions
delta_px_sq_scrap = delta_px_sq; 
delta_px_sq_scrap(delta_px_sq_scrap > 0) = 1; % Find active frames 
delta_px_sq_scrap = diff(delta_px_sq_scrap); % Diff to find transitions  
    % 1 = inactive to active 
    % -1 = active to inactive 
 
parfor f = 1:size(delta_px_sq,2) % For each fish 
        % Note this this runs apporximately twice as fast as just using a
        % For loop 
    
    % Starts - ensures no bouts are lost at the start  
    if  delta_px_sq(1,f) > 0 % If active in first bin  
        wake_cells{1,f}(:,1) = [1 ; find(delta_px_sq_scrap(:,f) == 1)+1]; % Find active bout starts
        sleep_cells{1,f}(:,1) = (find(delta_px_sq_scrap(:,f) == -1)+1); % Find sleep bout starts 
    else % Ie. if inactive in first bin 
        wake_cells{1,f}(:,1) = find(delta_px_sq_scrap(:,f) == 1)+1; % Find active bout starts 
        sleep_cells{1,f}(:,1) = [1 ; (find(delta_px_sq_scrap(:,f) == -1)+1)]; % Find sleep bout starts 
    end 
    
    % Ends - ensures no bouts are lost at the end 
    if delta_px_sq(end,f) > 0 % If active in last bin 
        wake_cells{1,f}(:,2) = [find(delta_px_sq_scrap(:,f) == - 1);...
            size(delta_px_sq,1)]; % Find active bout ends
        sleep_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == 1); % Find sleep bout ends
    else 
        wake_cells{1,f}(:,2) = find(delta_px_sq_scrap(:,f) == - 1); 
        sleep_cells{1,f}(:,2) = [(find(delta_px_sq_scrap(:,f) == 1)) ; size(delta_px_sq,1)]; % Find sleep bout ends
    end
    
    % Parameter extraction 
    wake_cells{1,f}(:,3:8) = NaN; % Pre-allocate 
    wake_cells{1,f}(:,3) = (wake_cells{1,f}(:,2)+1) - wake_cells{1,f}(:,1); % Wake Bout Length 
    sleep_cells{1,f}(:,3) = (sleep_cells{1,f}(:,2)+1) - sleep_cells{1,f}(:,1); % Sleep Bout Length 

    % Active bouts 
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        wake_cells{1,f}(b,4) = nanmean(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Mean
        wake_cells{1,f}(b,5) = nanvar(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Variance
        wake_cells{1,f}(b,6) = nansum(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Total
        wake_cells{1,f}(b,7) = nanmin(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Min
        wake_cells{1,f}(b,8) = nanmax(delta_px_sq(wake_cells{1,f}(b,1):...
            wake_cells{1,f}(b,2),f)); % Max
    end
    
end 
 
% Parameter Time 
for f = 1:size(delta_px_sq,2) % For each fish
    
    % Active Bouts 
    for b = 1:size(wake_cells{1,f},1) % For each active bout
        parameter_time{1}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,3); % Fill in bout length
        parameter_time{2}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,4); % Fill in bout Mean
        parameter_time{3}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,5); % Fill in bout Variance 
        parameter_time{4}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,6); % Fill in bout Total  
        parameter_time{5}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,7); % Fill in bout Minimum 
        parameter_time{6}(wake_cells{1,f}(b,1),f) = ...
            wake_cells{1,f}(b,8); % Fill in bout Maximum 
        parameter_time{7}(wake_cells{1,f}(b,1),f) = 1; % No. of Active bouts 
        parameter_time{8}(wake_cells{1,f}(b,1):wake_cells{1,f}(b,2),f) = 1; % Total Time Active  
        parameter_time{9}(wake_cells{1,f}(b,1):wake_cells{1,f}(b,2),f) = ...
            nansum(wake_cells{1,f}(1:b,6)); % Total Activity
    end
    
    % Inactive Bouts 
    for b = 1:size(sleep_cells{1,f},1) % For each sleep bout
        parameter_time{10}(sleep_cells{1,f}(b,1),f) = ...
            sleep_cells{1,f}(b,3); % Fill in bout length
        parameter_time{11}(sleep_cells{1,f}(b,1),f) = 1; % No. of Inactive Bouts
        parameter_time{12}(sleep_cells{1,f}(b,1):sleep_cells{1,f}(b,2),f) = 1; % Total Time Inactive
    end
    
    disp(horzcat('Calculated parameters across time for fish = ',num2str(f),...
    ' of ',num2str(size(delta_px_sq,2)))); % Report progress 
end

toc
clear b delta_px_sq_scrap f p   

%% Statistics & Plots - Variables 

% Hard Code your periods of interest 
days = [1 2 3 4]; % Hard code 
nights = [1 2 3]; % Hard code 

% Determine day/night order
    % Note that dn currently assumes the experiment starts during the day
lb = [1 ; lb]; % Add 1 to lb
if dn(1) == 1 % If the experiment starts in the day 
    lb_days = lb(1:2:size(lb,1)); % Assign day start values (in frames) 
    lb_nights = lb(2:2:size(lb,1)); % Assign night start values (in frames) 
    days_crop = 1:2:size(lb,1); nights_crop = 2:2:size(lb,1); 
            % Assign logical indicies  
else
    lb_days = lb(2:2:size(lb,1)); % Assign day start values (in frames) 
    lb_nights = lb(1:2:size(lb,1)); % Assign night start values (in frames) 
    days_crop = 2:2:size(lb,1); nights_crop = 1:2:size(lb,1); 
            % Assign logical indicies  
end 

% Determine which windows are of interest 
time_window(1) = min([days_crop(days) nights_crop(nights)]);  
time_window(2) = max([days_crop(days) nights_crop(nights)]); 

% Determine first night  
if min(days_crop(days)) < min(nights_crop(nights))
    first_night = 2; 
else 
    first_night = 1; 
end

% Colors 

cmap_2 = hsv(max(group_tags)*2); % Generate a 2x color map 
cmap = cmap_2(1:2:size(cmap_2,1),:); % Extract main colors 

night_color = [0.9608 0.9608 0.9608]; % For background  

% Choose plot binning 
time_bins = 60*15; 
    
% Group sizes 
for g = 1:max(group_tags) % For each group 
    group_sizes(g) = size(find(group_tags == g),1); 
end 

clear g

%% Parameters - Generating Averages (in frames) 
parameter_matrix = nan(size(wake_cells,2),size(parameters,2),...
    size(lb,1)); % Fish x parameters x time windows
parameter_indicies = cell(2,size(wake_cells,2)); % wake/sleep x fish
lb = [lb ; size(delta_px_sq,1)]; % Add end to lb

for f = 1:size(wake_cells,2) % For each fish
    for t = 1:size(parameter_matrix,3) % For each time window
        % Wake bouts
        clear time_start time_stop;
        % Find the first bout that starts within the window
        time_start = find(wake_cells{1,f}(:,1) >= lb(t),1,'first');
        % Find the last bout that starts within the window
        if t+1 < size(lb,1) % For most windows
            time_stop = find(wake_cells{1,f}(:,1) < lb(t+1),1,'last');
        else % For the last window
            time_stop = find(wake_cells{1,f}(:,1) <= lb(t+1),1,'last');
        end
        
        % Store logical index
        parameter_indicies{1,f} = [parameter_indicies{1,f} ; ...
            ones(size(time_start:time_stop,2),1)*t];
        
        % Extract bout parameters (1-6)
        parameter_matrix(f,1:(size(wake_cells{1,f},2)-2),t)...
            = nanmean(wake_cells{1,f}(time_start:time_stop,3:end)); % Means
        % Number of bouts (7)
        parameter_matrix(f,7,t) = size(wake_cells{1,f}...
            (time_start:time_stop,3:end),1);
        % Total time active (8) - sum of lengths
        parameter_matrix(f,8,t) = nansum(wake_cells{1,f}...
            (time_start:time_stop,3),1);
        % Total activity (9) - sum of activity
        parameter_matrix(f,9,t) = nansum(wake_cells{1,f}...
            (time_start:time_stop,6),1);
        
        % sleep bouts (10-12)
        clear time_start time_stop;
        % Find the first bout that starts within the window
        time_start = find(sleep_cells{1,f}(:,1) >= lb(t),1,'first');
        % Find the last bout that starts within the window
        if t+1 < size(lb,1) % For most windows
            time_stop = find(sleep_cells{1,f}(:,1) < lb(t+1),1,'last');
        else % For the last window
            time_stop = find(sleep_cells{1,f}(:,1) <= lb(t+1),1,'last');
        end
        
        % Store logical index
        parameter_indicies{2,f} = [parameter_indicies{2,f} ; ...
            ones(size(time_start:time_stop,2),1)*t];
        
        % Sleep Bout Length (10)
        parameter_matrix(f,10,t)...
            = nanmean(sleep_cells{1,f}(time_start:time_stop,3));
        % Number of bouts (11)
        parameter_matrix(f,11,t) = size(sleep_cells{1,f}...
            (time_start:time_stop,3:end),1);
        % Total time inactive (12) - sum of lengths
        parameter_matrix(f,12,t) = nansum(sleep_cells{1,f}...
            (time_start:time_stop,3),1);
        
    end
end

% Re-group data for ease of comparisons
parameter_comparisons = cell(1,size(parameters,2)); % Pre-allocate
for p = 1:size(parameter_comparisons,2) % For each parameter
    parameter_comparisons{p}(1:max(group_sizes),1:max(group_tags),...
        1:size(parameter_matrix,3)) = NaN; % {parameters} Most fish per group x
    % x each group x time windows
end

for p = 1:size(parameter_comparisons,2) % For each parameter
    for g = 1:max(group_tags) % For each group
        for t = 1:size(parameter_matrix,3) % For each time window
            parameter_comparisons{p}(1:group_sizes(g),g,t) = ...
                parameter_matrix(group_tags == g,p,t);
        end
    end
end

clear f t time_start time_stop p g

%% Smoothing data into seconds

% Pre-allocate
delta_px_sq_sec = nan(size(1:(fps):size(delta_px_sq,1),2),...
    size(delta_px_sq,2),'single'); % time (seconds) x fish
delta_px_sq_sec_smooth = nan(size(1:(fps):size(delta_px_sq,1),2),...
    size(delta_px_sq,2),'single'); % time (seconds) x fish
dn_sec = nan(size(1:(fps):size(delta_px_sq,1),2),1,'single'); % time x 1
for p = 1:size(parameters,2) % For each parameter
    parameter_time_sec_smooth{p} = nan(size(delta_px_sq_sec_smooth),'single');
    % {parameters} time (seconds) x fish
end

% Smooth each fish's activity & parameters
for f = 1:size(delta_px_sq,2) % For each fish
    
    a = 1; % Start a counter
    for t = 1:(fps):size(delta_px_sq,1) % For each second
        if t + (fps-1) < size(delta_px_sq,1) % Check to prevent running off the end
            delta_px_sq_sec(a,f) = nanmean(delta_px_sq(t:t+(fps-1),f)); % Bin activity
            
            for p = 1:size(parameters,2) % For each parameter
                if parameter_smooth(p) == 0 % For most parameters
                    parameter_time_sec_smooth{p}(a,f) = ...
                        nanmean(parameter_time{p}(t:t+(fps-1),f)); % Mean
                elseif parameter_smooth(p) == 1
                    parameter_time_sec_smooth{p}(a,f) = ...
                        nansum(parameter_time{p}(t:t+(fps-1),f)); % Sum
                elseif parameter_smooth(p) == 2
                    parameter_time_sec_smooth{p}(a,f) = ...
                        nanmax(parameter_time{p}(t:t+(fps-1),f)); % Max
                end
            end
            
            if f == 1 % For the first fish
                dn_sec(a,1) = mode(dn(t:t+(fps-1),1));
                % Take the most common light value within each bin
                % Will need to create indicies for each time window here
                % (04.08.17)
            end
        
        a = a + 1; % Add to counter
        
        else 
            delta_px_sq_sec(a,f) = 0; % This prevents smoothing errors
            dn_sec(a,1) = dn_sec(a-1,1); % Assume last value 
        end
        
    end
    
    % Smooth the activity data
    delta_px_sq_sec_smooth(:,f) = smooth(delta_px_sq_sec(:,f),time_bins);
    disp(horzcat('Smoothed fish ',num2str(f),' of ' ,...
        num2str(size(delta_px_sq,2)))); % Report progress
end

% Determine Day/Night Transitions in seconds
% Note that the binning will make this not 100% accurate
lb_sec = [1 ; find(diff(dn_sec) ~= 0) + 1; size(delta_px_sq_sec_smooth,1)];

clear a f p t 

%% Combining Experiments (< Parameter Extracted Data)

 %Load Data - using multiselect
 [filename, pathname] = uigetfile('*.mat', 'Select files','MultiSelect','on'); %Select files
 if isequal(filename,0) %If no file is selected
     error('No Files Selected') %Show Error
 else %If selected
     disp(['User selected ', fullfile(pathname, filename)]) %Show selected filenames
 end
 
 experiment_tags = [];
 parameter_time_sec_smooth = [];
 lb_sec = [];
 group_tags = [];
 wake_cells = [];
 sleep_cells = [];
 parameter_indicies = [];
 parameter_matrix = [];
 parameter_comparisons = cell(1,12);
 for f = 1:size(filename,2) %For each file
     clear experiment; 
     experiment = load(strcat(pathname,filename{f})); %Load the mat file
     experiment_tags = [experiment_tags ; ones(size(experiment.group_tags,1),1)*f]; 
     group_tags = [group_tags ; experiment.group_tags]; 
             % Allocate experiment tags 
     
     % Merge variables 
     parameter_time_sec_smooth = [parameter_time_sec_smooth ; experiment.parameter_time_sec_smooth];  
     lb_sec = [lb_sec experiment.lb_sec]; 
     wake_cells = [wake_cells experiment.wake_cells]; 
     sleep_cells = [sleep_cells experiment.sleep_cells]; 
     parameter_indicies = [parameter_indicies experiment.parameter_indicies]; 
     parameter_matrix = [parameter_matrix ; experiment.parameter_matrix]; 
     delta_px_sq_sec_smooth{f,1} = experiment.delta_px_sq_sec_smooth; 
     
     % Nab variables 
     if f == 1 % For the first file 
         parameters = experiment.parameters; 
         time_window = experiment.time_window; 
         time_bins = experiment.time_bins; 
         unit_conversion = experiment.unit_conversion; 
         cmap = experiment.cmap; 
         cmap_2 = experiment.cmap_2; 
         night_color = experiment.night_color; 
         nights = experiment.nights; 
         nights_crop = experiment.nights_crop; 
         geno_list = experiment.geno_list; 
         units = experiment.units; 
         units_2 = experiment.units_2; 
         parameter_smooth = experiment.parameter_smooth; 
         first_night = experiment.first_night; 
         days_crop = experiment.days_crop; 
         days = experiment.days;
     end 
        
     for p = 1:size(parameters,2) % For each parameter 
        parameter_comparisons{p} = [parameter_comparisons{p} ; ...
            experiment.parameter_comparisons{p}]; 
     end 
     
 end
 
 clear pathname filename geno f;

% Aligning data
    % Parameter time 
parameter_time_sec_smooth(end+1,:) = cell(1,size(parameters,2)); 
start = max(lb_sec(2,:));
for p = 1:size(parameters,2) % For each parameter 
    parameter_time_sec_smooth{end,p} = nan(max(lb_sec(end,:)),...
        size(group_tags,1),'single'); % Pre-allocate  
    
    for f = 1:max(experiment_tags) % For each experiment 
        parameter_time_sec_smooth{end,p}(start-lb_sec(2,f)+1:...
            start-lb_sec(2,f)+lb_sec(end,f),experiment_tags==f) = ...
            parameter_time_sec_smooth{f,p}; 
    end 
end 
parameter_time_sec_smooth = parameter_time_sec_smooth(end,:); % Keep merged data 

% Activity 
delta_px_sq_sec_smooth{end+1,1} = nan(max(lb_sec(end,:)),...
    size(group_tags,1),'single'); % Pre-allocate
for f = 1:max(experiment_tags) % For each experiment 
    delta_px_sq_sec_smooth{end,1}(start-lb_sec(2,f)+1:...
        start-lb_sec(2,f)+lb_sec(end,f),experiment_tags==f) = ...
        delta_px_sq_sec_smooth{f,1}; 
end 
delta_px_sq_sec_smooth = delta_px_sq_sec_smooth{end,1}; % Keep merged data 

% Store time windows 
[~,p] = find(lb_sec == max(lb_sec(end,:))); % Find new time windows 
lb_sec = lb_sec(:,p); % Keep longest time windows 

% Re-determine Group sizes
for g = 1:max(group_tags) % For each group
    group_sizes(g) = size(find(group_tags == g),1);
end

% Remove experiment tags - Optional  
 experiment_tags(:,1) = 1;

% Possible color shading 
%cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5))

% If just using a single experiment 
%experiment_tags = ones(size(group_tags)); 

% Alternative colors 
cmap(1,:) = [135 206 250]/255; % light sky blue
cmap_2(1,:) = cmap;
cmap_2(2,:) = [25 25 112]/255; % midnight blue 

% Selecting a time window 
days = [2 3]; nights = [2 3]; 
time_window(1) = min([days_crop(days) nights_crop(nights)]);  
time_window(2) = max([days_crop(days) nights_crop(nights)]); 

% Clean up
clear f p start  
%% Parameters Across Time 

% Removing time points with no data
for p = 1:size(parameters,2) % For each parameter
    scrap = []; % Clear scrap
    
    for g = 1:max(group_tags) % For each group
        scrap = [scrap ; find(sum(isnan(parameter_time_sec_smooth{p}(:,group_tags == g)),2) == ...
            size(find(group_tags == g),1))]; % Find time points without data
    end
    
    scrap = unique(scrap); % Clear repeated values
    scrap = sort(scrap); % Sort to ascending order
    
    parameter_time_sec_smooth_y{p}(:,1) = 1:size(parameter_time_sec_smooth{p},1); 
        % Generate a time-line
    parameter_time_sec_smooth_y{p}(scrap,:) = []; % Remove times with no data  
    parameter_time_sec_smooth{p}(scrap,:) = []; % Remove rows with no data 
    
end

clear scrap p g 

%% Distributions

% Calculating the boundaries for each fit 
dist_boundaries = nan(size(wake_cells,2)*2,size(parameters,2),'single');
    % fish*2 (max/min) x parameters 

a = 1; % Start a counter 
for f = 1:size(wake_cells,2) % For each fish  
    dist_boundaries(a:a+1,1:6) = minmax(wake_cells{1,f}(:,3:end)')'; % Wake Bout Parameters 
    dist_boundaries(a:a+1,7:9) = minmax((squeeze(parameter_matrix(f,7:9,:))')')'; % Totals  
    dist_boundaries(a:a+1,10) = minmax(sleep_cells{1,f}(:,3)')'; % Sleep Bout Length  
    dist_boundaries(a:a+1,11:12) = minmax((squeeze(parameter_matrix(f,11:12,:))')')'; % Totals
    a = a + 2; 
end 

% Pre-allocation 
parameter_dists = cell(1,size(parameters,2)); 
for p = find(parameter_smooth == 0) % For wake bout parameters 
    parameter_dists{p} = nan(size(wake_cells,2),size(min(dist_boundaries(:,p)):...
        max(dist_boundaries(:,p)),2),size(parameter_matrix,3),'single'); 
        % {parameters} fish x parameter range x time 
end 

tic
% Distribution Fitting
counter = 0; % Start counter  
for p = find(parameter_smooth == 0) % For most parameters
    
    if ismember(p,1:6) == 1  % For wake parameters
        for f = 1:size(wake_cells,2) % For each fish
            for t = 1:size(parameter_matrix,3) % For each time window                
                
                if isempty(find(parameter_indicies{1,f} == t)) == 0 % If there are bouts
                    pd = fitdist(wake_cells{1,f}(parameter_indicies{1,f}==t,p+2),'kernel','Width',1); % Fit
                    parameter_dists{p}(f,:,t) = pdf(pd,min(dist_boundaries(:,p)):max(dist_boundaries(:,p)));
                else % If there are no bouts
                    parameter_dists{p}(f,:,t) = zeros(1,size(min(dist_boundaries(:,p)):...
                        max(dist_boundaries(:,p)),2)); % Fill with zeros 
                end
            end
        end
              
    elseif p == 10 % For sleep bout length
        for f = 1:size(sleep_cells,2) % For each fish
            for t = 1:size(parameter_matrix,3) % For each time window
                
                if isempty(find(parameter_indicies{2,f} == t)) == 0 % If there are bouts
                    pd = fitdist(sleep_cells{1,f}(parameter_indicies{2,f}==t,3),'kernel','Width',1); % Fit
                    parameter_dists{p}(f,:,t) = pdf(pd,min(dist_boundaries(:,p)):max(dist_boundaries(:,p)));
                else % If there are no bouts 
                    parameter_dists{p}(f,:,t) = zeros(1,size(min(dist_boundaries(:,p)):...
                        max(dist_boundaries(:,p)),2)); % Fill with zeros 
                end
                
            end
        end
        
    end
    
    counter = counter + 1; % Add to counter 
    disp(horzcat('Parameter ',num2str(counter),' of ',...
        num2str(size(find(parameter_smooth == 0),2)),' Fit'));
end
 toc 
 
clear a counter f p pd t time_start time_stop

%% Figure - Parameter Means 

figure; 
for p = 1:size(parameter_comparisons,2) % For each parameter
    subplot(3,4,p); hold on; clear scrap; counter = 1; 
    title(parameters{p}); % Add title
    
    % Plot parameters
    for e = 1:max(experiment_tags) % For each experiment 
        for g = 1:max(group_tags) % For each group
            legend_lines(g) = errorbar((nanmean(squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))),1)/unit_conversion(1,p)),...
                (nanstd(squeeze(parameter_comparisons{p}(experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))))...
                /unit_conversion(1,p))./sqrt(size(find(experiment_tags == e & group_tags == g),1)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);

            if e == 1 % For the first experiment 
            legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
                num2str(size(find(group_tags == g),1)),')');
            % Append the group size to each group name
            end 
            
            % To accurately determine axis scale
            % Add and subtract the std from each mean then determine
            % The highest & lowest value for each group
            scrap(1,counter) = max(nanmean(squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))),1)/unit_conversion(1,p) + ...
                (nanstd(squeeze(parameter_comparisons{p}(experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))))...
                /unit_conversion(1,p))./sqrt(size(find(experiment_tags == e & group_tags == g),1))); 
            scrap(2,counter) = min(nanmean(squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))),1)/unit_conversion(1,p) - ...
                (nanstd(squeeze(parameter_comparisons{p}(experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))))...
                /unit_conversion(1,p))./sqrt(size(find(experiment_tags == e & group_tags == g),1)));
            
            counter = counter + 1; % Add to counter 
        end
    end 
    
    % Add night patches
    y_lims = [(min(scrap(2,:)) - min(scrap(2,:))*0.05) ...
        (max(scrap(1,:)) + max(scrap(1,:))*0.05)]; % Add a bit of space either side 
        
    a = 1; night_start = first_night; % Start counters  
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
            1 (y_lims(2)-y_lims(1))],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1; night_start = night_start + 2; % Add to counters
    end
    
    % Figure Looks
    if p == 4 % For the 4th parameter
        [~,~,~,~] = legend(legend_cell,'Location','northwest'); % Generate axis
        legend('boxoff'); % Turn legend off  
    end
    axis([0.5 size([days_crop(days) nights_crop(nights)],2)+0.5 ...
        y_lims]); % Set axis 
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format 
    xlabel('Time (Days/Nights)','Fontsize',12); % X Labels 
    set(gca, 'XTick', []); % Turn off X-Ticks 
    ylabel(units(p),'Fontsize',12); % Y Labels 

end
  
clear p scrap g legend_lines y_lims a n r night_start count

%% Figure - Parameters Across Time  

figure; hold on; clear legend_lines
for p = 1:size(parameters,2) % For each parameter
    subplot(3,4,p); hold on; title(parameters{p}); clear scrap; clear y;
    
    y = find(parameter_time_sec_smooth_y{p} >= lb_sec(time_window(1))...
        & parameter_time_sec_smooth_y{p} < lb_sec(time_window(2)+1)); 

    for g = 1:max(group_tags) % For each group
        
        legend_lines(g) = plot(parameter_time_sec_smooth_y{p}(y,1),...
            smooth(nanmean(parameter_time_sec_smooth{p}(y,group_tags == g),2)...
            ,time_bins)/unit_conversion(2,p),'color',cmap(g,:)); % Plot 
                
        legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
            num2str(size(find(group_tags == g),1)),')');
        % Append the group size to each group name
            
        % Find the top & bottom 
            % Excluding values that smooth outside of the data
        scrap(1,g) = max(legend_lines(g).YData((time_bins):...
            (size(legend_lines(g).YData,2) - time_bins)));
        scrap(2,g) = min(legend_lines(g).YData((time_bins):...
            (size(legend_lines(g).YData,2) - time_bins)));
        
    end
    
    a = 1;
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[lb_sec(nights_crop(nights(n))) ...
            min(scrap(2,:)) - (min(scrap(2,:))*0.05)...
            (lb_sec(nights_crop(nights(n))+1)-1) - lb_sec(nights_crop(nights(n)))...
            max(scrap(1,:)) + (max(scrap(1,:))*0.05)],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1;
    end
    
    if p == 4 % For the 4th parameter - add a legend 
        [~,~,~,~] = legend(legend_cell,'Location','northwest');
        legend('boxoff'); 
    end
    
    axis([(parameter_time_sec_smooth_y{p}(y(1),1) + time_bins)...
        (parameter_time_sec_smooth_y{p}(y(end),1) - time_bins) ...
        min(scrap(2,:)) - (min(scrap(2,:))*0.05) ...
        max(scrap(1,:)) + (max(scrap(1,:))*0.05)]);
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12);
    xlabel('Time (Days/Nights)','Fontsize',12);
    set(gca, 'XTick', []);
    ylabel(units_2(p),'Fontsize',12);
end

clear a g legend_lines legend_cell n p r scrap y 

%% Figure - Activity   

figure; hold on; clear legend_lines  
for e = 1:max(experiment_tags) % For each experiment
    for g = 1:max(group_tags) % For each group
        
        legend_lines(g) = shadedErrorBar(lb_sec(time_window(1)):...
            (lb_sec(time_window(2)+1)-1),nanmean(delta_px_sq_sec_smooth...
            (lb_sec(time_window(1)):(lb_sec(time_window(2)+1)-1),group_tags == g & experiment_tags == e),2),...
            nanstd(delta_px_sq_sec_smooth(lb_sec(time_window(1)):(lb_sec(time_window(2)+1)-1),...
            group_tags == g & experiment_tags == e)')/sqrt(size(find(group_tags == 1 & experiment_tags == 1),1)),...
            {'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5))});
        
        legend_cols(g) = legend_lines(g).mainLine; % Store color
        
        legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
            num2str(size(find(group_tags == g),1)),')');
        % Append the group size to each group name
        
        if g == max(group_tags) && e == max(experiment_tags) % After the last group
            a = 1; % Start counter
            for n = 1:size(nights,2) % For each night
                y_lims = ylim; % Find the axis limits
                r(a) = rectangle('Position',[lb_sec(nights_crop(nights(n))) 0,...
                    (lb_sec(nights_crop(nights(n))+1)-1) - lb_sec(nights_crop(nights(n))) y_lims(2)],...
                    'FaceColor',night_color,'Edgecolor',[1 1 1]);
                uistack(r(a),'bottom'); % Send to back
                a = a + 1; % Add to counter
            end
        end
    end
end 
box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
[~,icons,plots,~] = legend(legend_cols,legend_cell,'Location','northwest');
legend('boxoff'); 
set(icons(1:max(group_tags)),'Fontsize',32) ; set(plots,'LineWidth',3);
axis([lb_sec(time_window(1)) (lb_sec(time_window(2)+1)-1) 0 y_lims(2)]);
x = lb_sec(time_window(1)):(60*60*12):(lb_sec(time_window(2)+1)-1); 
set(gca,'XTick',x); 
set(gca,'XTickLabel',{(0:size(x,2)-1)*12})
xlabel('Time (Hours)','Fontsize',32);
ylabel('Delta Px (a.u)','Fontsize',32);

clear a g h icons plots str legend_cell legend_cols legend_lines n r x y_lims 

%% Figure - Parameter Distributions V3
    % Uses a log axis for plots 
    
figure; 
for p = 1:size(parameters,2) % For each parameter
    subplot(3,4,p); hold on; col = 1; counter = 1; 
    clear data legend_lines legend_cols legend_cell r sample crop y_lims spread_cols;
    title(parameters{p}); % Add title
    
    if ismember(p,[7:9 11:12]) == 0 % For most parameters 
        
        % Take day & night distributions across time windows  
        data{1,1} = reshape(permute(parameter_dists{p}(:,:,days_crop(days)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,days_crop(days)),2));
        data{1,2} = reshape(permute(parameter_dists{p}(:,:,nights_crop(nights)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,nights_crop(nights)),2));
        
        % Expand group tags to account for multiple days/nights 
        data{2,1} = repmat(group_tags,[size(data{1,1},1)/size(group_tags,1),1]);
        data{2,2} = repmat(group_tags,[size(data{1,2},1)/size(group_tags,1),1]);
        
        % Expand experiment tags to account for multiple days/nights 
        data{3,1} = repmat(experiment_tags,[size(data{1,1},1)/size(experiment_tags,1),1]);
        data{3,2} = repmat(experiment_tags,[size(data{1,2},1)/size(experiment_tags,1),1]);
        
        % V2
        % Find a cut off for each figure 
            % This helps with visualisation & prevents the generation of 
                % large patch objects which often crash Matlab 
            %crop = find(nanmean([data{1,1} ; data{1,2}]) < 0.01,1,'first');
        
       % V3 
       crop = size(data{1,1},2); % crop = all data
       
       for e = 1:max(experiment_tags) % For each experimet 
           col = 1;  
        for g = 1:max(group_tags) % For each group
            for t = 1:2 % For day/night
                if t == 1 % for day 
                    sample = size(find(data{2,t} == g & data{3,t} == e),1)/...
                        size(days,2);
                else % for night
                    sample = size(find(data{2,t} == g & data{3,t} == e),1)/...
                        size(nights,2);
                end
                    % Find the number of fish
                
                legend_lines(col) = shadedErrorBar((1:crop)/unit_conversion(1,p),...
                    nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,1:crop)),...
                    nanstd(data{1,t}(data{2,t} == g & data{3,t} == e,1:crop))/...
                    sqrt(sample),...
                    {'color',cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5))}); 
                % Plot - note that this is now in appropriate units
                    % (eg.seconds) 
                
                legend_cols(col) = legend_lines(col).mainLine; % Save color

                % Append time window to group names 
                if e == 1 % For the first experiment 
                    if first_night == 2 % If starting in the day
                        if t == 1 % For the days
                            legend_cell{col} = horzcat(geno_list.colheaders{g},' - Day'); % Tag
                        else % For the nights
                            legend_cell{col} = horzcat(geno_list.colheaders{g},' - Night'); % Tag
                        end
                    else % If starting in the night
                        if t == 1 % For the nights
                            legend_cell{col} = horzcat(geno_list.colheaders{g},' - Night'); % Tag
                        else % For the days
                            legend_cell{col} = horzcat(geno_list.colheaders{g},' - Day'); % Tag
                        end
                    end
                end
                % Determine Axis boundaries 
                y_lims(1,counter) = max(legend_lines(col).mainLine.YData);
                y_lims(2,counter) = min(legend_lines(col).mainLine.YData);
                
                col = col + 1; % Add to color 
                counter = counter + 1; % Add to counter 
            end
            
        end
       end 
        axis([1/unit_conversion(1,p) crop/unit_conversion(1,p) ...
            min(y_lims(2,:)) max(y_lims(1,:))]); % Set axis limits 
        set(gca,'XTick',...
            [1/unit_conversion(1,p), 1/unit_conversion(1,p)*10,...
                crop/unit_conversion(1,p)]); % set x tick labels 
        % Set decimal places depending on units 
        if unit_conversion(1,p) > 1 
             xtickformat('%.2f');
        else 
            xtickformat('%.0f');
        end 
        set(gca,'XScale','log'); % set log axis 
        xlabel(units(p),'Fontsize',12); % X labels 
        ylabel('Probability Density','Fontsize',12); % Y label  
    
    else % For the other parameters 
        for e = 1:max(experiment_tags)
          spread_cols = plotSpread(squeeze(nanmean(permute(parameter_comparisons{p}(experiment_tags(group_tags == g) == e,:,days_crop(days)),[1 3 2]),2))/...
                unit_conversion(1,p),...
                'distributionColors',cmap_2(1:2:size(cmap_2,1),:)+(1-cmap_2(1:2:size(cmap_2,1),:))*(1-(1/e^.5)),'showMM',2); % Scatter day data  
            spread_cols{2}.LineWidth = 3; % Change marker width 
            spread_cols = plotSpread(squeeze(nanmean(permute(parameter_comparisons{p}(experiment_tags(group_tags == g) == e,:,nights_crop(nights)),[1 3 2]),2))/...
                unit_conversion(1,p),...
                'distributionColors',cmap_2(2:2:size(cmap_2,1),:)+(1-cmap_2(2:2:size(cmap_2,1),:))*(1-(1/e^.5)),'showMM',2); % Scatter night data 
            spread_cols{2}.LineWidth = 3; % Change marker width
        end 
        ylabel(units(p),'Fontsize',12); % Y labels 
        set(gca,'xticklabel',geno_list.colheaders,'Fontsize',12); % Name each group  
    end
    if p == 4 % For the 4th parameter - add a legend
        [~,~,~,~] = legend(legend_cols,legend_cell,'Location','northeast');
        legend('boxoff'); 
    end
    
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12);
    
end

clear col g p spread_cols t xl

%% Stats - Two Way ANOVA 

% Define groups 
anova_group = []; % group  
anova_experiment = []; % experiments 
for g = 1:max(group_tags) % For each group 
    anova_group = [anova_group  ones(max(group_sizes),1)*g]; % Allocate tag 
    anova_experiment = [anova_experiment ; experiment_tags(group_tags == g)]; % Experiment tags
end 
anova_group = repmat(anova_group(:)',[1,size([days_crop(days) nights_crop(nights)],2)]); % Expand for all days/nights  
anova_experiment = repmat(anova_experiment',[1,size([days_crop(days) nights_crop(nights)],2)]); % Experiment 

anova_time = []; % time 
for t = time_window(1):time_window(2) % For each time window 
    anova_time = [anova_time ones(max(group_sizes),max(group_tags))*mod(t,2)]; 
    % Allocate alternating zeros and ones to each time window 
end 
anova_time = anova_time(:)'; % Vectorise 

if size(days_crop(days),2) == size(nights_crop(nights),2) % If there are an equal number of windows 
    anova_development = []; % development
    anova_development = zeros(1,size(anova_group,2)); % Pre-allocate 
    d = 1:size(anova_development,2)/(size(time_window(1):time_window(2),2)/2):...
        size(anova_development,2); % divide into "24h" windows 
    for t = 1:size(d,2)-1 
        anova_development(d(t):d(t+1)-1) = t; 
    end 
end 

% Pre-allocation (Need to test to find x)  
for p = 1:size(parameters,2) % For each parameter 
    twa.p(1:15,p) = NaN; % Comparisons x parameters  
    twa.stats{p} = []; % Stats structure 
end 

% Calculations - account for experiment tags
for p = 1:size(parameters,2) % For each parameter
    clear scrap;
    scrap = parameter_comparisons{p}(:,:,time_window(1):time_window(2));
    scrap = scrap(:)'; % Vectorise  
    
    if p == 1 % For the first parameter remove NaN values 
        anova_group(isnan(scrap)) = []; 
        anova_time(isnan(scrap)) = []; 
        if exist('anova_development','var') == 1 % Check if variable exists 
            anova_development(isnan(scrap)) = []; 
        end 
    end 
    
    scrap(isnan(scrap)) = []; 
    
    if size(days_crop(days),2) == size(nights_crop(nights),2) % If comparing development 
        if max(experiment_tags) > 1 % If comparing experiments 
            [twa.p(:,p),~,twa.stats{p}] = anovan(scrap,...
                {anova_group,anova_time,anova_development,anova_experiment},'display','off','model','full');
        else % Development but no experiments 
            [twa.p(1:7,p),~,twa.stats{p}] = anovan(scrap,...
                {anova_group,anova_time,anova_development},'display','off','model','full');
        end
    else % Without development 
        if max(experiment_tags) > 1 % With experiments 
            [twa.p(1:7,p),~,twa.stats{p}] = anovan(scrap,...
                {anova_group,anova_time,anova_experiment},'display','off','model','full'); % Try without
        else % Without experiments 
            [twa.p(1:3,p),~,twa.stats{p}] = anovan(scrap,...
                {anova_group,anova_time},'display','off','model','full'); % Try without
        end
        
    end
    
end

clear anova_development anova_group anova_time anova_experiment p scrap t

%% Stats - KW 

for p = 1:size(parameters,2) % For each parameter
    clear data; c = 1; d = 1; % Start counters 
    if ismember(p,[7:9 11:12]) == 0 % For most parameters 
        data{1,1} = reshape(permute(parameter_dists{p}(:,:,days_crop(days)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,days_crop(days)),2)); % Day data
        data{1,2} = reshape(permute(parameter_dists{p}(:,:,nights_crop(nights)),[1 3 2]),...
            [],size(parameter_dists{p}(:,:,nights_crop(nights)),2)); % Night data
        
        data{2,1} = repmat(group_tags,[size(data{1,1},1)/size(group_tags,1),1]); % Day tags
        data{2,2} = repmat(group_tags,[size(data{1,2},1)/size(group_tags,1),1]); % Night tags
        
        % Expand experiment tags to account for multiple days/nights
        data{3,1} = repmat(experiment_tags,[size(data{1,1},1)/size(experiment_tags,1),1]);
        data{3,2} = repmat(experiment_tags,[size(data{1,2},1)/size(experiment_tags,1),1]);
        
        % Compare Day vs Night for each group
        for e = 1:max(experiment_tags) % For each experiment 
            for g = 1:max(group_tags) % For each group
                [kw.time_h(d,p),kw.time_p(d,p)] = kstest2(nanmean(data{1,1}(data{2,1} == g & data{3,1} == e,:)),...
                    nanmean(data{1,2}(data{2,2} == g & data{3,2} == e,:))); % Compare day vs night distributions
                d = d + 1; 
                for n = g+1:max(group_tags) % for each group comparison
                    for t = 1:2 % For each time window
                        [kw.group_h{p}(c,t),kw.group_p{p}(c,t)] = kstest2(nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,:)),...
                            nanmean(data{1,t}(data{2,t} == n & data{3,t} == e,:))); % Make group-wise comparisons
                    end
                    c = c + 1;
                end
            end
        end 
        
    else
        data{1,1} = squeeze(nanmean(permute(parameter_comparisons{p}(:,:,days_crop(days)),[1 3 2]),2)); % Day data
        data{1,2} = squeeze(nanmean(permute(parameter_comparisons{p}(:,:,nights_crop(nights)),[1 3 2]),2)); % Night data
             
        for e = 1:max(experiment_tags) % For each experiment
            for g = 1:max(group_tags) % For each group
                [kw.time_h(d,p),kw.time_p(d,p)] = kstest2(data{1,1}(experiment_tags(group_tags == g) == e,g),...
                    data{1,2}(experiment_tags(group_tags == g) == e,g)); % Compare day vs night
                d = d + 1; % add to counter
                for n = g+1:max(group_tags) % for each group comparison
                    for t = 1:2 % For each time window
                        [kw.group_h{p}(c,t),kw.group_p{p}(c,t)] = kstest2(data{1,t}(experiment_tags(group_tags == g) == e,g),...
                            data{1,t}(experiment_tags(group_tags == g) == e,n)); % Make each comparison
                    end
                    c = c + 1; % add to counter
                end
            end
        end
    end
end
