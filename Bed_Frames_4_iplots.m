% Bed_Frames_4_i_plots

%% Figure - Parameter Means (I) 

figure; 
for p = 1:size(parameter_comparisons,2) - 2 % For each parameter
    subplot(2,5,p); hold on; clear scrap; counter = 1; 
    set(gca,'FontName','Calibri'); % Set Font  
    title(parameters{p}); % Add title
    
    % Plot parameters
    for e = 1:max(experiment_tags) % For each experiment 
        for g = 1:max(group_tags) % For each group
      
            plot((squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2)))/unit_conversion(1,p))',...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/(e+4)^.5)),'linewidth',1.5)
            
            legend_lines(g) = plot((nanmean(squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))),1)/unit_conversion(1,p)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);

            if e == 1 % For the first experiment 
            legend_cell{g} = horzcat(geno_list.colheaders{g},', n = ',...
                num2str(size(find(group_tags == g),1)));
            % Append the group size to each group name
            end 
            
            % To accurately determine axis scale
            % Add and subtract the std from each mean then determine
            % The highest & lowest value for each group
            scrap(1,counter) = max(max((squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2)))/unit_conversion(1,p))));
            scrap(2,counter) = min(min((squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2)))/unit_conversion(1,p))));
            
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
    if p == 5 % For the 4th parameter
        [~,~,~,~] = legend(legend_cell,'Location','northwest'); % Generate axis
        legend('boxoff'); % Turn legend off  
    end
    axis([1 size([days_crop(days) nights_crop(nights)],2) ...
        y_lims]); % Set axis 
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format 
    xlabel('Time (Days/Nights)','Fontsize',12); % X Labels 
    set(gca, 'XTick', []); % Turn off X-Ticks 
    ylabel(units(p),'Fontsize',12); % Y Labels 

end
  
clear p scrap g legend_lines y_lims a n r night_start count

%% Figure - Activity (Cropped) 
    % Crops the activity so that only points with data from all experiments
        % are shown 

figure; hold on; clear legend_lines; set(gca,'FontName','Calibri'); % Set Font  
for e = 1:max(experiment_tags) % For each experiment
    for g = 1:max(group_tags) % For each group
        
        legend_lines(g) = shadedErrorBar(lb_sec(time_window(1)):...
            (lb_sec(time_window(2)+1)-1),nanmean(delta_px_sq_sec_smooth...
            (lb_sec(time_window(1)):(lb_sec(time_window(2)+1)-1),group_tags == g & experiment_tags == e),2),...
            nanstd(delta_px_sq_sec_smooth(lb_sec(time_window(1)):(lb_sec(time_window(2)+1)-1),...
            group_tags == g & experiment_tags == e)')/sqrt(size(find(group_tags == g & experiment_tags == e),1)),...
            'lineprops',{'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5))});
        
        legend_cols(g) = legend_lines(g).mainLine; % Store color
        
        legend_cell{g} = horzcat(geno_list.colheaders{g},', n = ',...
            num2str(size(find(group_tags == g),1)));
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
axis([find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'first') + time_bins...
    find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'last')- time_bins 0 y_lims(2)]);

x = find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'first') + time_bins:...
    (60*60*12):find(sum(isnan(delta_px_sq_sec_smooth'))==0,1,'last')- time_bins; 
set(gca,'XTick',x); 
set(gca,'XTickLabel',{(0:size(x,2)-1)*12})
xlabel('Time (Hours)','Fontsize',32);
ylabel('Delta Px','Fontsize',32);

clear a g h icons plots str legend_cell legend_cols legend_lines n r x y_lims 


%% Figure - Parameter Distributions V3 (I)
% Uses a log axis for plots
% Plot's Std Rather than SEM

figure;
for p = 1:size(parameters,2) - 2 % For each parameter
    subplot(2,5,p); hold on; col = 1; counter = 1;
    clear data legend_lines legend_cols legend_cell r sample crop y_lims spread_cols;
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12);
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
                    % SEM
                    %                 if t == 1 % for day
                    %                     sample = size(find(data{2,t} == g & data{3,t} == e),1)/...
                    %                         size(days,2);
                    %                 else % for night
                    %                     sample = size(find(data{2,t} == g & data{3,t} == e),1)/...
                    %                         size(nights,2);
                    %                 end
                    %                     Find the number of fish
                    
                    plot((1:crop)/unit_conversion(1,p),...
                        data{1,t}(data{2,t} == g & data{3,t} == e,1:crop),...
                        'color',cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/(e+4)^.5)),'linewidth',1.5)
                    
                    legend_lines(col) = plot((1:crop)/unit_conversion(1,p),...
                        nanmean(data{1,t}(data{2,t} == g & data{3,t} == e,1:crop)),...
                        'color',cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5)),'linewidth',3);
                    
                    % Plot - note that this is now in appropriate units
                    % (eg.seconds)
                    
                    legend_cols(col,:) = legend_lines(col).Color; % Save color
                    
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
                    y_lims(1,counter) = max(legend_lines(col).YData);
                    y_lims(2,counter) = min(legend_lines(col).YData);
                    
                    col = col + 1; % Add to color
                    counter = counter + 1; % Add to counter
                end
                
            end
        end
        axis([1/unit_conversion(1,p) crop/unit_conversion(1,p) ...
            min(y_lims(2,:)) max(y_lims(1,:))]); % Set axis limits
         try
            set(gca,'XTick',...
                [1/unit_conversion(1,p), 1/unit_conversion(1,p)*10,...
                crop/unit_conversion(1,p)]); % set x tick labels
        catch
            set(gca,'XTick',...
                [1/unit_conversion(1,p),(0.5*crop)/unit_conversion(1,p),...
                crop/unit_conversion(1,p)]); % set x tick labels      
        end
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
        col = 1;
        for g = 1:max(group_tags) % for each group
            clear data;
            % Day
            data{1} = parameter_comparisons{p}(:,g,days_crop(days));
            data{1}(isnan(data{1})) = [];
            data{1} = nanmean(reshape(data{1},[group_sizes(g),size(days_crop(days),2)]),2)/...
                unit_conversion(1,p);
            
            % Night
            data{2} = parameter_comparisons{p}(:,g,nights_crop(nights));
            data{2}(isnan(data{2})) = [];
            data{2} = nanmean(reshape(data{2},[group_sizes(g),size(nights_crop(nights),2)]),2)/...
                unit_conversion(1,p);
            
            % Plot
            for e = 1:max(experiment_tags)
                spread_cols = plotSpread(data{1}(experiment_tags(group_tags == g) == e),...
                    'xValues',g,'distributionColors',...
                    cmap_2(col,:)+(1-cmap_2(col,:))*(1-(1/e^.5)),'showMM',2);
                spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
                spread_cols = plotSpread(data{2}(experiment_tags(group_tags == g) == e),...
                    'xValues',g,'distributionColors',...
                    cmap_2(col+1,:)+(1-cmap_2(col+1,:))*(1-(1/e^.5)),'showMM',2);
                spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = [1 0.5 0]; % Change marker properties
            end
            
            col = col + 2;
            
        end
        ylabel(units(p),'Fontsize',12); % Y labels
        set(gca,'xtick',1:max(group_tags)); % Set x ticks
        set(gca,'xticklabel',geno_list.colheaders,'Fontsize',12); % Name each group    end
        %     if p == 5 % For the 5th parameter - add a legend
        %         [~,icons,~,~] = legend(legend_cols,legend_cell,'Location','northeast');
        %         legend('boxoff');
        %         set(icons(1:max(group_tags)*2),'Fontsize',13);
        %         set(icons((max(group_tags)*2)+1:2:end),'LineWidth',3)
        %     end
        
    end
end 
    clear col g p spread_cols t xl icons data 