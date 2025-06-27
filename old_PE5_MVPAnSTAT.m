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
cd('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');

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
%%
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
stat_meeg=[];stat_time=[];
for i=1:6
    for idev=2
%
cfg = [];
cfg.metric          = 'auc';
cfg.test            = 'permutation';
cfg.correctm        = 'cluster';
cfg.n_permutations  = 10000;
cfg.clustercritval  = 1.96;
cfg.alpha           =0.01;
% for alpha 0.05 = t-val is 1.96
% for alpha 0.01 = t-val is 2.58

% Level 2 stats settings
cfg.design          = 'within';
cfg.statistic       = 'ttest';
cfg.null            = 0.5;



        stat_meeg{i,idev} = mv_statistics(cfg,result_meeg(i,:,idev));
        stat_time{i,idev} = mv_statistics(cfg,result_time(i,:,idev));
        
    end
end
save([filler1,filler2,filler3,filler4,method,'_alltime_mvpa_permutation1x100_2.mat'],'stat_time','stat_meeg');

%% plot
%% plot time generalization
 load('debsln_0_1resam_lda_eachtwoday_MVPA_result_alltime1x100_2.mat')
 load('debsln_0_1resam_lda_alltime_mvpa_permutation1x100_2.mat')
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
 iplt=1;
%subplot(2,2,3)
for i=[1:3,5]
    for idev=2
        subplot('Position',pos_main(iplt,:));
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
        caxis([0.46 0.535]);
        
        grid off
        text(0.59,-0.075,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
        colormap(pink)
        title([plottitle{i}],'FontSize',10,'FontWeight','bold');
        
        if iplt==4
            xlabel('Test time(s)')
            ylabel('Train time(s)')
        end
        if iplt==8
            colorbar;
            c = colorbar;
            c.Title.String = 'AUC';
            
        end
        caxis([0.46 0.535]);
        set(gca, 'xtick', 0:0.1:0.5);
        set(gca, 'ytick', 0:0.1:0.5);
        set(gca,'FontSize',10,'FontWeight','bold')
        
        % 反向 Y 轴
        set(gca,'YDir','normal');
        
        %title('Time generalization','position',[-0.25 0.8])
        
    end
end
% mv_plot_result(mean_result, predata.time,'mask',mask','new_figure',0)

%%
mycolor=[215,48,39;252,141,89;254,224,144;150,150,150;145,191,219;69,117,180]./255;

legend_labels=[];
figure
iplt=1;
%subplot(2,2,1)
for i=[1:3,5]
    subplot(2,2,iplt)
    iplt=iplt+1;
    
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
    
    
    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    mask=stat_meeg{i,idev}.mask;
    isig=find(mask==1);
    plot(xval(isig), dat(isig),'color',mycolor(i,:),'Linewidth',4);hold on
    line([0 0],[0.46 0.56],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    line([-0.1 0.6],[0.5 0.5],'linestyle','--', 'Color',[0.5 0.5 0.5], 'LineWidth',1);hold on
    
    %[h.plt, h.patch] = boundedline(predata.time,mean_result{i}.perf,mean_result{i}.perf_std{1}, 'alpha','LineWidth',1.5,'Color',mycolor(i,:));hold on
    xval=predata.time;
    dat=(mean_result{i}.perf{1})';
    tmp=(mean_result{i}.perf_std{1})';
    
    % 计算填充区域的上下界限
    upper_bound = dat + tmp;
    lower_bound = dat - tmp;
    
    % 填充标准差区域
    fill([xval fliplr(xval)],[lower_bound fliplr(upper_bound)],mycolor(i,:),'edgealpha', '0', 'facealpha', '.1');hold on
    
    text(0.59,0.465,txt{idev},'HorizontalAlignment','right','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
    text(-0.09,0.465,[num2str(xval(isig(1))*1000),'-',num2str(xval(isig(end))*1000),'ms'],'HorizontalAlignment','left','color',[0.5 0.5 0.5],'FontSize',10,'FontWeight','bold')
    
    % legend_labels{5}=legend_labels{3};
    % legend_labels{3}=legend_labels{2};
    % legend_labels{2}=
    if iplt==4
        xlabel('Times');
        ylabel('AUC');
    end
    set(gca,'YLim',[0.46,0.56])
    set(gca,'XLim',[-0.1,0.6])
    %legend(legend_labels,'box','off')
    title([plottitle{i}],'FontSize',10,'FontWeight','bold');
    set(gca, 'xtick', 0:0.1:0.5);
    set(gca,'FontSize',10,'FontWeight','bold')
    set(gca,'Box','off')
    
end
