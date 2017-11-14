% Video Bout Features Figure 

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
figure; hold on; 
plot(scrap,'color',cmap(1,:),'linewidth',3); 

% Plotting - Marks 
scatter(example_frames,scrap(example_frames),72,...
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

% Labels 
text(5,(0.5*10^4)+1250,'Length','Fontsize',32); 
text(13.5,nanmean(scrap(11:13)),{'Length','Mean','Variance','Total'},'Fontsize',32); 
text(12.25,max(scrap),'Maximum','Fontsize',32); 
text(13.25,scrap(13),'Minimum','Fontsize',32); 

% Nice Figure 
axis([1 size(scrap,1) ylim]); box off; set(gca,'Fontsize',32); 
xticks(2:2:size(scrap,1)); set(gca,'XTickLabel',{(1:2:size(scrap,1))/fps}); 
xlabel('Time (seconds)'); 
yticks([]); ylabel('Delta Px (a.u)'); 

%% Old 
% Active 
subplot(1,2,1); title('Active Bout Features'); 
a = find(states(2,:) == 13,1,'last'); % Hard chosen bout 
hold on; clear scrap; 
scrap = raw_data(2,a-7:a+4); 
plot(scrap,'color',cmap(1,:),'linewidth',3); 
scatter(7,max(scrap),72,'k','filled'); 
text(7.5,double(max(scrap)),'Maximum','Fontsize',16); 
scatter(8,min(scrap(scrap>0)),72,'k','filled'); 
text(8.5,double(min(scrap(scrap>0))),'Minimum','Fontsize',16); 
plot([6 8], [mean(scrap(scrap > 0)) mean(scrap(scrap > 0))],'-k',...
    'linewidth',3); 
plot([6 6],[mean(scrap(scrap > 0))-0.5 mean(scrap(scrap > 0))+0.5],'k',...
    'linewidth',3); 
plot([8 8],[mean(scrap(scrap > 0))-0.5 mean(scrap(scrap > 0))+0.5],'k',...
    'linewidth',3); 
text(8.5,10.5,{'Length','Mean','Variance','Total'},'Fontsize',16); 
set(gca,'Fontsize',16); 
axis([1 12 0 15.5]); 
xticks(1:2:12); 
set(gca,'XTickLabel',{(1:2:12)/fps}); 
xlabel('Time (seconds)'); 
ylabel('Delta Px (a.u)'); 

% Inactive 
subplot(1,2,2); hold on; title('Inactive Bout Features'); set(gca,'Fontsize',16);  
a = find(states(2,:) == 2,1,'last'); % Hard chosen bout 
hold on; clear scrap;
scrap = raw_data(2,a-18:a+1);
plot(scrap,'color',cmap(1,:),'linewidth',3); 
axis([1 20 ylim]); 
plot([2 19],[7 7],'-k','linewidth',3); 
plot([2 2],[6.5 7.5],'k','linewidth',3); 
plot([19 19],[6.5 7.5],'k','linewidth',3); 
text(10,7.5,'Length','Fontsize',16); 
xticks(1:2:20); 
set(gca,'XTickLabel',{(1:2:20)/fps}); 
xlabel('Time (seconds)'); 
ylabel('Delta Px (a.u)'); 