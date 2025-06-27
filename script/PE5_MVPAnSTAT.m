%%%% 导入的数据应为单个被试keep trial的数据
clear;clc;

method = 'lda';

is100trial=0;
isremovebaseline=0;
ishighpass=0;
isresample=1;
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
if isresample==1
    filler4='resam_';
else
    filler4='';
end
filename={[filler1,filler2,filler3,'day1_nsubDEVnSTDmvpa.mat'],[filler1,filler2,filler3,'day2_nsubDEVnSTDmvpa.mat'],...
    [[filler1,filler2,filler3],'day3_nsubDEVnSTDmvpa.mat'],[filler1,filler2,filler3,'day4_nsubDEVnSTDmvpa.mat']};
eeglab
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');

trialcluster=1;%每trialcluster个trial平均
%%
for md=1:4
    
    load(filename{md});
    for subj=1:14
        
        for idev=1:2
            data4clssfy{md,subj, idev}=clssi_Diff{subj, idev};
            
            tmp=[];
            clss_ntrial=floor(size(data4clssfy{md,subj, idev}.sampleinfo,1)/trialcluster);
            tmp=data4clssfy{md,subj, idev}.trial;
            tmp=reshape(tmp(1:trialcluster.*clss_ntrial,:,:),trialcluster,clss_ntrial,9,175);
            tmp=squeeze(mean(tmp,1));
            
            cfg=[];
            cfg.trials=1:clss_ntrial;
            data4clssfy{md,subj, idev}=ft_selectdata(cfg,data4clssfy{md,subj, idev});
            
            data4clssfy{md,subj, idev}.trial = tmp;
            
            
        end
    end
end

if isresample==1
    for md=1:4
        for subj=1:14
            for idev=1:2
                cfg = [];
                cfg.resamplefs = 50;
                cfg.method='downsample';
                data4clssfy{md,subj, idev} = ft_resampledata(cfg, data4clssfy{md,subj, idev});
            end
        end
    end
end
%% classification using MVPA-LIGHT
stat=[];
for idev=1:2
    i=1;
    for i1stday=1:3
        for i2ndday=i1stday+1:4
            nnj=14;
            for nn=1:nnj
                i1stday
                i2ndday
                nn
                for iter=1
                    predata=[];
                    predata=data4clssfy{i1stday,nn,idev};postdata=data4clssfy{i2ndday,nn,idev};
                    npre=[];npost=[];
                    npre=size(predata.trial,1);
                    npost=size(postdata.trial,1);
                    
                    cfg = [];
                    %                 cfg.time            = 101:151;
                    cfg.classifier      = method;
                    cfg.metric          = 'auc';
                    cfg.preprocess      = 'demean';
                    % to be a bit faster, use 10-fold cross-validation with no extra repetitions
                    cfg.k               = 5;
                    %cfg.repeat          = 3;
                    
                    cfg.repeat          = 100;
                    
                    
                    design=[ones(npre,1); 2*ones(npost,1)];
                    data=cat(1,predata.trial, postdata.trial);
                    [cf_meeg{i,nn,idev}, result_meeg{i,nn,idev,iter}] = mv_classify_across_time(cfg, data,design);
                    idev
                    i1stday
                    i2ndday
                    nn
                   
                    %[cf_time{i,nn,idev}, result_time{i,nn,idev,iter}] = mv_classify_timextime(cfg, data,design);
                end
                
            end
          i=i+1; 
        end 
    end
   
end
save([filler1,filler2,filler3,filler4,method,'_eachtwoday_MVPA_result_meeg.mat'],'result_meeg');

%%
% stat_meeg=[];stat_time=[];
load('resam_lda_eachtwoday_MVPA_result_temporal.mat')
 load('resam_lda_eachtwoday_MVPA_result_meeg.mat')
 MMN_time=result_time;MMN_meeg=result_meeg; P3a_time=result_time;P3a_meeg=result_meeg;
 MMN_win=13:18;P3a_win=19:23;
 for i=1:6
    for idev=1:2
        for isub=1:14
            MMN_meeg{i,isub,idev}.perf=result_meeg{i,isub,idev}.perf(MMN_win,1);
            MMN_meeg{i,isub,idev}.perf_std=result_meeg{i,isub,idev}.perf_std(MMN_win,1);
            P3a_meeg{i,isub,idev}.perf=result_meeg{i,isub,idev}.perf(P3a_win,1);
            P3a_meeg{i,isub,idev}.perf_std=result_meeg{i,isub,idev}.perf_std(P3a_win,1);
            MMN_time{i,isub,idev}.perf=result_time{i,isub,idev}.perf(MMN_win,:);
            MMN_time{i,isub,idev}.perf_std=result_meeg{i,isub,idev}.perf_std(MMN_win,:);
            P3a_time{i,isub,idev}.perf=result_time{i,isub,idev}.perf(P3a_win,:);
            P3a_time{i,isub,idev}.perf_std=result_meeg{i,isub,idev}.perf_std(P3a_win,:);
        end
    end
 end
%%
for i=5:6
    for idev=1:2
%
cfg = [];
cfg.metric          = 'auc';
cfg.test            = 'permutation';
cfg.correctm        = 'cluster';
cfg.n_permutations  = 10000;
cfg.clustercritval  = 1.96;
cfg.alpha           =0.05/6;
% for alpha 0.05 = t-val is 1.96
% for alpha 0.01 = t-val is 2.58

% Level 2 stats settings
cfg.design          = 'within';
cfg.statistic       = 'ttest';
cfg.null            = 0.5;



        stat_meeg_MMN{i,idev} = mv_statistics(cfg,MMN_meeg(i,:,idev));
        stat_time_MMN{i,idev} = mv_statistics(cfg,MMN_time(i,:,idev));
        stat_meeg_P3a{i,idev} = mv_statistics(cfg,P3a_meeg(i,:,idev));
        stat_time_P3a{i,idev} = mv_statistics(cfg,P3a_time(i,:,idev));
        
    end
end
save([filler1,filler2,filler3,filler4,method,'_alltime_mvpa_permutation.mat'],'stat_time_MMN','stat_meeg_MMN','stat_time_P3a','stat_meeg_P3a');

%% plot
%% plot time generalization
load('resam_lda_eachtwoday_MVPA_result_temporal.mat')
 load('resam_lda_alltime_mvpa_permutation.mat')
mask=[];
mean_result=[];
plottitle={'Day1 vs Day2',' Day1 vs Day3',' Day1 vs Day4','Day2 vs Day3',' Day2 vs Day4',' Day3 vs Day4'};
results=[];
txt={'Small deviant','Large deviant'};
figure
for i=1:4
    a1=subplot(2,2,i);
    pos_main(i,:)=get(a1,'position');
end
pos_main(:,3)=0.3;pos_main(:,4)=0.3;
w=0.35;
x1=0.13;x2=0.13+w;x3=x2+w;
ind2=find(pos_main==pos_main(2,1));
pos_main(ind2)=x2;
h=0.32;
y1=0.5;y2=y1-h;
yind1=find(pos_main==0.5838);
pos_main(yind1)=y1;
yind2=find(pos_main==0.1100);
pos_main(yind2)=y2;
 iplt=4;
%subplot(2,2,3)
figure
for idev=[1,2]
    if idev==1
        ii=3;
    elseif idev==2
        ii=2:3;
    end
    for i=ii

        %subplot('Position',pos_main(iplt,:));
        subplot(2,3,iplt)
        iplt=iplt+1;
        results=result_time(i,:,idev);
        for isub=1:14
            results{isub}= mv_prepare_plot(results{isub});
        end
        mean_result{i} = mv_combine_results(results, 'average');
        
        mean_result{i}.plot{1}.title=[plottitle{i}];
        %mean_result{i}.plot{1}.title=['Day3 vs Day4'];
        
        legend_labels{i}=[plottitle{i}];
        
        mean_result{i}.metric='AUC';
        perf=cell2mat(mean_result{i}.perf(1));
        mask=zeros(35,35);
        mask(P3a_win,:)=double(stat_time_P3a{i, idev}.mask);
        mask=mask';
        perf = perf.*mask;
        perf(perf==0) = nan;
        imagesc([-0.1, 0.6],[-0.1, 0.6],perf');
        %mv_plot_result(mean_result{i}, predata.time, predata.time,'mask',stat_time{i, idev}.mask  )
       
        
        grid off
        text(0.59,-0.075,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',19,'FontWeight','bold')
        % 自定义白色到酒红色的渐变色
       n = 256; 
r = linspace(1, 0.5, n);
g = linspace(1, 0, n);  
b = linspace(1, 0, n);   

cmap = [r' g' b']; 
colormap(cmap);   

        title([plottitle{i}],'FontSize',19,'FontWeight','bold');
                ax = gca;
  
    pos = ax.Position;
    pos(1) = pos(1) * 0.9; 
    if iplt>3
    pos(2) = pos(2) + 0.1;  % 向上移动一些，减少间距
    end
    ax.Position = pos;
        
        if iplt==5
            xlabel('Testing time (s)')
            ylabel('Training time (s)')
        end
        if iplt==7
            colorbar;
            c = colorbar;
            c.Title.String = 'AUC';


% 获取当前颜色条的位置并右移
pos = c.Position; 
pos(1) = pos(1) + 0.1; % 将颜色条往右移 0.05 个单位
c.Position = pos; % 更新颜色条的位置
            
        end
        caxis([0.50 0.535]);
        set(gca, 'xtick', 0:0.1:0.6);
        set(gca, 'ytick', 0:0.1:0.6);

        set(gca,'FontSize',19,'FontWeight','bold')
        set(gca,'ticklength',[0.018 0.015]);
        % 反向 Y 轴
        set(gca,'YDir','normal');
             
      %  title('Time generalization','position',[-0.25 0.8])
               
        
    end
end
% mv_plot_result(mean_result, predata.time,'mask',mask','new_figure',0)

%% temporal decoding_small
legendnam=[];
load('resam_lda_eachtwoday_MVPA_result_temporal.mat')
 load('resam_lda_alltime_mvpa_permutation.mat')
 load('resam_lda_eachtwoday_MVPA_result_meeg.mat')
 plottitle={'Day1 vs Day2','Day1 vs Day3','Day1 vs Day4','Day2 vs Day3','Day2 vs Day4','Day3 vs Day4'};

%mycolor=[255,0,0;55,126,184;77,175,74;152,78,163;255,127,0;100,100,100]./255;
mycolor=[215,25,28;215,25,28;253,174,28;253,174,97;171,217,233;44,123,182]./255;

legend_labels=[];
figure
iplt=1;
%subplot(2,2,1)    
subplot(2,2,1)
ii=0;
for i=[2,3,5,6]

ii=ii+1;
    
    for idev=1
        
        
        results=result_meeg(i,:,idev);
        for isub=1:14
            results{isub}= mv_prepare_plot(results{isub});
        end
        mean_result{i} = mv_combine_results(results, 'average');
    end
        legendnam{ii}=plottitle{i};
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    predata=data4clssfy{1,1,1};
    xval=predata.time;
    
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    %legend_labels{i}=['Day1 vs',daynam{i}];
    % 绘制曲线
    plot(xval, dat ,'color',mycolor(i,:),'Linewidth',1);
    hold on;
    
end
    for i=[2:3,5:6]
    for idev=1
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    mask1=stat_meeg_MMN{i,idev}.mask;
    isig1=find(mask1==1)+12;
    mask2=stat_meeg_P3a{i,idev}.mask;
    isig2=find(mask2==1)+18;
    plot(xval(isig1), dat(isig1),'color',mycolor(i,:),'Linewidth',3.5);hold on
    plot(xval(isig2), dat(isig2),'color',mycolor(i,:),'Linewidth',3.5);hold on

    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    
    % 计算填充区域的上下界限
    upper_bound = dat + tmp;
    lower_bound = dat - tmp;
    
    % 填充标准差区域
    %fill([xval fliplr(xval)],[lower_bound fliplr(upper_bound)],mycolor(i,:),'edgealpha', '0', 'facealpha', '.1');hold on
    
    %text(0.59,0.485,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
   % text(0.59,0.485,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
   % text(-0.09,0.485,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
     
   
    end

    end
     text(0.18,0.535,'MMN','HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
      text(0.28,0.535,'P3a','HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
    line([0 0],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([-0.1 0.6],[0.5 0.5],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([0.15 0.15],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([0.25 0.25],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([0.35 0.35],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    ax = gca;


    pos = ax.Position;
    pos(3) = pos(3) * 0.8; 

    pos(2) = pos(2) - 0.1;  % 向上移动一些，减少间距

    ax.Position = pos;
    % legend_labels{5}=legend_labels{3};
    % legend_labels{3}=legend_labels{2};
    % legend_labels{2}=


        ylabel('AUC');
        xlabel('Time (s)');


    set(gca,'YLim',[0.48,0.54])
    set(gca,'XLim',[-0.1,0.6])
    %legend(legend_labels,'box','off')
    leg=legend(legendnam,'Location','northeast','Orientation', 'vertical','FontSize', 11,'fontweight','bold','box','off');
  % set(leg, 'LineWidth', 3); 
    leg.ItemTokenSize = [99,19];
   h=title('Small deviant','FontSize',11,'FontWeight','bold');
    h.Position(2) = h.Position(2) + 0.001; 
   xticks([0,0.1, 0.2,0.3,0.4,0.5,0.6]); % 设置 X 轴的刻度位置
        xticklabels({'0','0.1', '0.2','0.3','0.4','0.5','0.6'}); % 设置对应的刻度标签
    set(gca,'FontSize',15,'FontWeight','bold')
            set(gca,'ticklength',[0.018 0.015]);
    set(gca,'Box','off')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%temporal decoding_large
load('resam_lda_eachtwoday_MVPA_result_temporal.mat')
 load('resam_lda_alltime_mvpa_permutation.mat')
 load('resam_lda_eachtwoday_MVPA_result_meeg.mat')
 plottitle={'Day1 vs Day2','Day1 vs Day3','Day1 vs Day4','Day2 vs Day3','Day2 vs Day4','Day3 vs Day4'};

%mycolor=[228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;100,100,100]./255;

legend_labels=[];

iplt=1;
%subplot(2,2,1)    
subplot(2,2,2)
for i=[2,3,5,6]


    
    for idev=2
        
        
        results=result_meeg(i,:,idev);
        for isub=1:14
            results{isub}= mv_prepare_plot(results{isub});
        end
        mean_result{i} = mv_combine_results(results, 'average');
    end
    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    predata=data4clssfy{1,1,1};
    xval=predata.time;
    
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    %legend_labels{i}=['Day1 vs',daynam{i}];
    % 绘制曲线
    plot(xval, dat ,'color',mycolor(i,:),'Linewidth',1);
    hold on;
    
end
    for i=[2:3,5:6]
    for idev=2
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    mask1=stat_meeg_MMN{i,idev}.mask;
    isig1=find(mask1==1)+12;
    mask2=stat_meeg_P3a{i,idev}.mask;
    isig2=find(mask2==1)+18;
    plot(xval(isig1), dat(isig1),'color',mycolor(i,:),'Linewidth',3.5);hold on
    plot(xval(isig2), dat(isig2),'color',mycolor(i,:),'Linewidth',3.5);hold on

    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    
    % 计算填充区域的上下界限
    upper_bound = dat + tmp;
    lower_bound = dat - tmp;
    
    % 填充标准差区域
    %fill([xval fliplr(xval)],[lower_bound fliplr(upper_bound)],mycolor(i,:),'edgealpha', '0', 'facealpha', '.1');hold on
    
    %text(0.59,0.485,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
   % text(0.59,0.485,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
   % text(-0.09,0.485,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
     end
    end
     text(0.18,0.535,'MMN','HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
      text(0.28,0.535,'P3a','HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
    line([0 0],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([-0.1 0.6],[0.5 0.5],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([0.15 0.15],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([0.25 0.25],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([0.35 0.35],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    ax = gca;


    pos = ax.Position;
    pos(3) = pos(3) * 0.8; 

    pos(2) = pos(2) - 0.1;  % 向上移动一些，减少间距
 pos(1) = pos(1) - 0.11;  % 向上移动一些，减少间距
    ax.Position = pos;
    % legend_labels{5}=legend_labels{3};
    % legend_labels{3}=legend_labels{2};
    % legend_labels{2}=




    set(gca,'YLim',[0.48,0.54])
    set(gca,'XLim',[-0.1,0.6])
    %legend(legend_labels,'box','off')
   % leg=legend(plottitle,'Location','northeast','Orientation', 'vertical','FontSize', 10,'fontweight','bold','box','off');
    leg.ItemTokenSize = [10,10];
   h=title('Large deviant','FontSize',10,'FontWeight','bold');
    h.Position(2) = h.Position(2)+0.001; 
   xticks([0,0.1, 0.2,0.3,0.4,0.5,0.6]); % 设置 X 轴的刻度位置
        xticklabels({'0','0.1', '0.2','0.3','0.4','0.5','0.6'}); % 设置对应的刻度标签
    set(gca,'FontSize',15,'FontWeight','bold')
            set(gca,'ticklength',[0.018 0.015]);
    set(gca,'Box','off')
%% SMALL_DTA
load('resam_lda_eachtwoday_MVPA_result_temporal.mat')
 load('resam_lda_alltime_mvpa_permutation.mat')
 load('resam_lda_eachtwoday_MVPA_result_meeg.mat')
mycolor=[215,48,39;252,141,89;254,224,144;150,150,150;145,191,219;69,117,180]./255;

legend_labels=[];
figure
iplt=4;
%subplot(2,2,1)
for i=[6,5,3]
    subplot(2,3,iplt)
    iplt=iplt+1;
    
    for idev=1
        
        
        results=result_meeg(i,:,idev);
        for isub=1:14
            results{isub}= mv_prepare_plot(results{isub});
        end
        mean_result{i} = mv_combine_results(results, 'average');
    end
    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    predata=data4clssfy{1,1,1};
    xval=predata.time;
    
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    %legend_labels{i}=['Day1 vs',daynam{i}];
    % 绘制曲线
    plot(xval, dat ,'color',mycolor(i,:),'Linewidth',1.5);
    hold on;
    
    
    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    mask=stat_meeg{i,idev}.mask;
    isig=find(mask==1); 
    if i==3
    plot(xval(18:25), dat(18:25),'color',mycolor(i,:),'Linewidth',3.5);hold on
    plot(xval(27:30), dat(27:30),'color',mycolor(i,:),'Linewidth',3.5);hold on
    else
    plot(xval(isig), dat(isig),'color',mycolor(i,:),'Linewidth',3.5);hold on
   
    end
        
    line([0 0],[0.48 0.54],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([-0.1 0.6],[0.5 0.5],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
        ax = gca;
  
    pos = ax.Position;
    pos(4) = pos(4) * 0.8; 
    ax.Position = pos;
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    
    % 计算填充区域的上下界限
    upper_bound = dat + tmp;
    lower_bound = dat - tmp;
    
    % 填充标准差区域
    fill([xval fliplr(xval)],[lower_bound fliplr(upper_bound)],mycolor(i,:),'edgealpha', '0', 'facealpha', '.1');hold on
    
    %text(0.59,0.485,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
    text(0.59,0.485,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',8,'FontWeight','bold')
   % text(-0.09,0.485,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
    
    % legend_labels{5}=legend_labels{3};
    % legend_labels{3}=legend_labels{2};
    % legend_labels{2}=
    if iplt==5% || iplt==2

        ylabel('AUC');
    end
    set(gca,'YLim',[0.48,0.54])
    set(gca,'XLim',[-0.1,0.6])
    %legend(legend_labels,'box','off')
    title([plottitle{i}],'FontSize',8,'FontWeight','bold');
   xticks([0,0.1, 0.2,0.3,0.4,0.5,0.6]); % 设置 X 轴的刻度位置
        xticklabels({'0','0.1', '0.2','0.3','0.4','0.5','(s)'}); % 设置对应的刻度标签
        set(gca,'ticklength',[0.018 0.015]);
    set(gca,'FontSize',8,'FontWeight','bold')
    set(gca,'Box','off')
    
end
%% plot time generalization 
load('resam_lda_eachtwoday_MVPA_result_temporal.mat')
 load('resam_lda_alltime_mvpa_permutation.mat')
mask=[];
mean_result=[];
plottitle={'Day1 vs Day2',' Day1 vs Day3',' Day1 vs Day4','Day2 vs Day3',' Day2 vs Day4',' Day3 vs Day4'};
results=[];
txt={'Small deviant','Large deviant'};
figure
%subplot(2,2,3)
iplt=1;
for ii=1:3
    if ii==1
     idev=2;i=2;
    elseif ii==2
      idev=2;i=3;
    else ii==3
        idev=1;i=3;
    end
        subplot(2,3,iplt+3);
        %subplot(2,2,pos_main(iplt,:))
        iplt=iplt+1;
        results=result_time(i,:,idev);
        for isub=1:14
            results{isub}= mv_prepare_plot(results{isub});
        end
        mean_result{i} = mv_combine_results(results, 'average');
        
        mean_result{i}.plot{1}.title=[plottitle{i}];
        %mean_result{i}.plot{1}.title=['Day3 vs Day4'];
        
        legend_labels{i}=[plottitle{i}];
        
        mean_result{i}.metric='AUC';
        perf=cell2mat(mean_result{i}.perf(1));
        mask=double(stat_time{i, idev}.mask);
        perf = perf.*mask;
        perf(perf==0) = nan;
        imagesc([-0.1, 0.6],[-0.1, 0.6],perf');
        %mv_plot_result(mean_result{i}, predata.time, predata.time,'mask',stat_time{i, idev}.mask  )
       
        
        grid off
        text(0.59,-0.075,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',19,'FontWeight','bold')
       n = 256; 
r = linspace(1, 0.5, n);
g = linspace(1, 0, n);  
b = linspace(1, 0, n);   

cmap = [r' g' b']; 
colormap(cmap);    
        title([plottitle{i}],'FontSize',19,'FontWeight','bold');
        
        if iplt==2
            xlabel('Test time(s)')
            ylabel('Train time(s)')
        end
        if iplt==9
            colorbar;
            c = colorbar;
            c.Title.String = 'AUC';
            
        end
        caxis([0.50 0.535]);
        set(gca, 'xtick', 0:0.1:0.6);
        set(gca, 'ytick', 0:0.1:0.6);

        set(gca,'ticklength',[0.018 0.015]);
        set(gca,'FontSize',19,'FontWeight','bold')

        set(gca,'YDir','normal');
        % 设置刻度线只在左边和下边显示

        
        %title('Time generalization','position',[-0.25 0.8])
        
    end

% mv_plot_result(mean_result, predata.time,'mask',mask','new_figure',0)