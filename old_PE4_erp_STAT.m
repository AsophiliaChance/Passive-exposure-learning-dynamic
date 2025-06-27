%% load files
clc;clear
cd('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');
is100trial=0;
isremovebaseline=0;
ishighpass=0;
%%
if is100trial==1
    filler1='100trial_';
else
    filler1='';
end
if isremovebaseline==1
   filler2='debsln_';
else
     filler2='';
end
if ishighpass==0.1
   filler3='0_1';
else
     filler3='';
end
filename={[filler1,filler2,filler3,'day1_nsubavg.mat'],[filler1,filler2,filler3,'day2_nsubavg.mat'],...
   [[filler1,filler2,filler3],'day3_nsubavg.mat'],[filler1,filler2,filler3,'day4_nsubavg.mat']};
%%
for md=1:4
    load(filename{md});
    %reshape Ã¿¸ö±äÁ¿
    for isub=1:14
        for idev=1:2
    Diff{md,isub,idev}=Diff_avg{isub,idev};
        end
    end
end
%% statistics
% cfg_neighb        = [];
% cfg_neighb.method = 'distance';
% neighbours        = ft_prepare_neighbours(cfg_neighb,  STD2timelock4{1});
cfg = [];

layout = ft_prepare_layout(cfg,  STD_avg{1,1} );

cfg.method = 'triangulation';
% % cfg.neighbourdist = 0.7;
cfg.layout = layout;
cfg.feedback = 'yes'; %%with this you get the feedback plot
neighbours = ft_prepare_neighbours(cfg, STD_avg{1,1} );
%% P3

cfg=[];
cfg.latency          = [0.25 0.35];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = -1;
cfg.clustertail      = -1;
cfg.alpha            =0.1/6;
cfg.numrandomization = 'all';
cfg.minnbchan        = 2;

Nsubj  = 14;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

for idev=1:2
    i=1;
    for iday1=1:3
        for iday2=iday1+1:4
            P3_stat{i,idev}= ft_timelockstatistics(cfg, Diff{iday1,:,idev}, Diff{iday2,:,idev});
            i=i+1;
        end
    end
end
% MMN
cfg=[];
cfg.latency          = [0.15 0.25];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            =0.05/6;
cfg.numrandomization = 'all';
cfg.minnbchan        = 2;

Nsubj  = 14;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

for idev=1:2
    i=1;
    for iday1=1:3
        for iday2=iday1+1:4
            MMN_stat{i,idev}= ft_timelockstatistics(cfg, Diff{iday1,:,idev}, Diff{iday2,:,idev});
            i=i+1;
        end
    end
end
%%

%% effect size
% grand average for effect sizes
eeglab

cfg = [];
cfg.keepindividual = 'yes';
for md=1:4
    for idev=1:2
        Diff_gavg {md,idev}      = ft_timelockgrandaverage(cfg, Diff{md,:,idev});
    end
    
end
%
cfg = [];
cfg.latency = [0.25 0.35];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
roiFIC = ft_selectdata(cfg, Diff_gavg {1,2});
roiFC  = ft_selectdata(cfg, Diff_gavg {4,2});
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD

Nsubj  = 14;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvareffect_roi_unpaired   = 1;
cfg.ivar   = 2;
effect_roi_unpaired = ft_timelockstatistics(cfg, roiFIC, roiFC);
%% grand average for plotting
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
for md=1:4
    for idev=1:2
Diff_gavg {md,idev}      = ft_timelockgrandaverage(cfg, Diff{md,:,idev});
    end

end
%% plotting
figure
for idev=1
for iday=1:4
plot(Diff_gavg{1, 1}.time,Diff_gavg{iday, idev}.avg(3,:));hold on;
end
end
legend
%% Figure 1
%% plot waveform for 4days
chaninx=[4,8,5,1,3,2,6,9,7];

location_child= [0.148 0.85 0.086 0.09;...
    0.424898148148147 0.832145684877276 0.0889629629629622 0.0918448139350752;...
    0.705916666666664 0.832145684877276 0.0889629629629622 0.0918448139350752;...
    0.14375 0.535773555027712 0.0889629629629622 0.0918448139350752;...
    0.424898148148147 0.535773555027712 0.0889629629629622 0.0918448139350752;...
    0.705916666666664 0.535773555027712 0.0889629629629622 0.0918448139350752;...
    0.14375 0.232486144101346 0.0889629629629622 0.0918448139350752;...
    0.424898148148147 0.232486144101346 0.0889629629629622 0.0918448139350752;...
    0.705916666666664 0.232486144101346 0.0889629629629622 0.0918448139350752];

pos_main=[    0.1300    0.7093    0.2134    0.2157;0.4108    0.7093    0.2134    0.2157;...
    0.6916    0.7093    0.2134    0.2157;0.1300    0.4096    0.2134    0.2157;...
    0.4108    0.4096    0.2134    0.2157;0.6916    0.4096    0.2134    0.2157;...
    0.1300    0.1100    0.2134    0.2157;0.4108    0.1100    0.2134    0.2157;...
    0.6916    0.1100    0.2134    0.2157];
pos_main(:,3)=0.23;pos_main(:,4)=0.23;
w=0.27;
x1=0.13;x2=0.13+w;x3=x2+w;
ind2=find(pos_main==0.4108);
pos_main(ind2)=x2;
ind3=find(pos_main==0.6916);
pos_main(ind3)=x3;
location_child(:,1)=pos_main(:,1)+0.016;
location_child(:,2)=pos_main(:,2)+0.14;
h = figure('color','w');

for i=1:length(chaninx)
    
    data=[];
    subplot('Position',pos_main(i,:));
    
    
    axis([-0.1 0.6 -2.0 2.0]);
    plottime=Diff_gavg{1, 2}.time;
    
    % indices of electrode for plotting and ROI
    chanind = chaninx(i);
    channam = label(chanind);
    channalnam{i}=channam;
    single_chan =chanind;
    
    
    idev=2;
%four day data    
    plot(plottime, Diff_gavg{1, idev}.avg(single_chan,:), 'Color',([145,191,219]./255),'LineWidth',1);hold on;
    plot(plottime,Diff_gavg{2, idev}.avg(single_chan,:), 'Color',([254 224 144]./255),'LineWidth', 1);hold on;    
    plot(plottime,Diff_gavg{3, idev}.avg(single_chan,:), 'Color',[0.75 0.75 0.75],'LineWidth', 1);hold on; 
    plot(plottime,Diff_gavg{4, idev}.avg(single_chan,:), 'Color',([252,197,192]./255),'LineWidth', 1);hold on;

%mark

x=P3_stat{2,idev}.time;
for ii =1:length(x)-1
    if P3_stat{2,idev}.mask(single_chan,ii) == 1
        
        plot(x(ii:ii+1), [2.2 2.2],'Color',[215 48 39]./255,'LineWidth', 2);  
    end
        
        if P3_stat{3,idev}.mask(single_chan,ii) == 1
        plot(x(ii:ii+1), [2.3 2.3],'Color',[69 117 180]./255,'LineWidth', 2);  
        end
 end


box off;

graphname=fullfile([label{single_chan}]);

legend1='day1';
legend2='day2';
legend3='day3';
legend4='day4';
title(graphname,'Position', [-0.1, 2.3, 0],'fontweight','bold');

if i==7
    ylabel('Amplitude (\muV)');
    xlabel('Time (ms)');
    leg=legend(legend1,legend2,legend3,legend4,'Location','SouthOutside','Orientation', 'vertical','FontSize', 8,'fontweight','bold','box','off');
    leg.ItemTokenSize = [9,9];
end

set(gca,'FontSize',9,'fontweight','bold');

axis([-0.100 0.6 -2.4 2.3]);
set(h,'PaperPositionMode','auto');
end

%%
