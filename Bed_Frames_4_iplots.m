% Bed_Frames_4_i_plots

%% Figure - Parameter Means 

figure; 
for p = 1:size(parameter_comparisons,2) % For each parameter
    subplot(3,4,p); hold on; clear scrap; counter = 1; 
    title(parameters{p}); % Add title
    
    % Plot parameters
    for e = 1:max(experiment_tags) % For each experiment 
        for g = 1:max(group_tags) % For each group
      
            plot((squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2)))/unit_conversion(1,p))',...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/(e+4)^.5)))
            
            legend_lines(g) = plot((nanmean(squeeze(parameter_comparisons{p}...
                (experiment_tags(group_tags == g) == e,g,time_window(1):time_window(2))),1)/unit_conversion(1,p)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);

            if e == 1 % For the first experiment 
            legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
                num2str(size(find(group_tags == g),1)),')');
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
    if p == 4 % For the 4th parameter
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

