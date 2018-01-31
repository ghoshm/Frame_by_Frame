%% Video Bout Features Figure 

% Import Data as matrix from 
% "D:\Behaviour\SleepWake\Videos\171013_18_19\Bonsai\Crop_dp.txt"

scrap = Cropdp(495:508); 
scrap(scrap < 256) = 0; 
scrap = [scrap ; zeros(6,1)]; 

% Specify Colors 
cmap(1,:) = [135 206 250]/255; % light sky blue
cmap_2(1,:) = cmap;
cmap_2(2,:) = [25 25 112]/255; % midnight blue 

% Specify frame rate 
fps = 25; 

% Specify "Example Frames" 
example_frames = 10:14; 

% Plotting - Data 
figure; hold on;  set(gca, 'FontName', 'Calibri'); 
plot(scrap,'color',cmap(1,:),'linewidth',3); 

% Plotting - Marks 
scatter(example_frames,scrap(example_frames),108,...
    'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:));
plot([2 10],[(0.5*10^4),(0.5*10^4)],'k','linewidth',3); 
plot([2 2],[(0.5*10^4)-1250,(0.5*10^4)+1250],'k','linewidth',3); 
plot([10 10],[(0.5*10^4)-1250,(0.5*10^4)+1250],'k','linewidth',3); 
plot([11 13],[nanmean(scrap(11:13)),nanmean(scrap(11:13))],...
    'k','linewidth',3); 
plot([11 11],[nanmean(scrap(11:13))-1250,nanmean(scrap(11:13))+1250],...
        'k','linewidth',3); 
plot([13 13],[nanmean(scrap(11:13))-1250,nanmean(scrap(11:13))+1250],...
        'k','linewidth',3); 
plot([16.5 16.5],[scrap(13) max(scrap)],'k','linewidth',3);
plot([16 16.5],[scrap(13) scrap(13)],'k','linewidth',3); 
plot([16 16.5],[max(scrap) max(scrap)],'k','linewidth',3); 

% Labels 
text(6,(0.5*10^4)+1850,{'Inactive ','Bout Length'},'Fontsize',32,...
        'FontWeight','Bold','HorizontalAlignment','center'); 
text(13.5,nanmean(scrap(11:13)),{'Length','Mean','Variance','Total'},'Fontsize',32); 
text(12.25,max(scrap),'Maximum','Fontsize',32); 
text(13.25,scrap(13),'Minimum','Fontsize',32); 
text(18.25,nanmean(scrap(13):max(scrap)),{'Active Bout','Parameters'},'Fontsize',32,...
    'FontWeight','Bold','HorizontalAlignment','center'); 

% Nice Figure 
axis([1 size(scrap,1) ylim]); box off; set(gca,'Fontsize',32); 
xticks(2:2:size(scrap,1)); set(gca,'XTickLabel',{(1:2:size(scrap,1))/fps}); 
xlabel('Time (Seconds)'); 
yticks([]); ylabel('Delta Px'); 
