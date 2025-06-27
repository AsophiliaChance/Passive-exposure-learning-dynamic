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
filename={[filler1,filler2,filler3,'day1_nsubDEVnSTDmvpa_re.mat'],[filler1,filler2,filler3,'day2_nsubDEVnSTDmvpa_re.mat'],...
    [[filler1,filler2,filler3],'day3_nsubDEVnSTDmvpa_re.mat'],[filler1,filler2,filler3,'day4_nsubDEVnSTDmvpa_re.mat']};
eeglab
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');

trialcluster=1;%每trialcluster个trial平均
%%
for md=1:4
    
    load(filename{md});
    for subj=1:14
        
        for idev=1:2
            data=DEV_keeptrials_nsub{subj, idev};
            cfg = [];
            %cfg.channel ={'C3','C4','Cz','F3','F4','P3','P4','Fz','Pz'};
            cfg.channel ={'F3','F4','Fz',};
            data_selected = ft_selectdata(cfg, data);
            
            data4clssfy{md,subj, idev}=data_selected;
            
            tmp=[];
            clss_ntrial=floor(size(data4clssfy{md,subj, idev}.sampleinfo,1)/trialcluster);
            tmp=data4clssfy{md,subj, idev}.trial;
            tmp=reshape(tmp(1:trialcluster.*clss_ntrial,:,:),trialcluster,clss_ntrial,3,175);
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
time=data4clssfy{md,subj, idev}.time;
%% classification using MVPA-LIGHT
stat=[];
for idev=2
    i=1;
    for i1stday=3
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
                    [cf_meeg, result_meeg{nn}] = mv_classify_across_time(cfg, data,design);
                    idev
                    i1stday
                    i2ndday
                    nn
                    
                    [cf_time, result_time{nn}] = mv_classify_timextime(cfg, data,design);
                end
                
            end
            i=i+1;
            save([method,'dev_',num2str(idev),num2str(i1stday),num2str(i2ndday),'_eachtwoday_MVPA_result_meeg_dev_re3.mat'],'result_meeg','result_time');
            
        end
        
    end
    
end
%% temporal decoding单样本t检验 % 16-18 20-23
win1=15:18; win2=20:23;
for idev=1:2
    i=1;
    for i1stday=1:3
        for i2ndday=i1stday+1:4
            data=[];mean_result=[];
            load([method,'dev_',num2str(idev),num2str(i1stday),num2str(i2ndday),'_eachtwoday_MVPA_result_meeg_std_re3.mat']);
            results=result_meeg;
            for isub=1:14
                dataMMN(isub)=mean( results{isub}.perf(win1));
                dataP3a(isub)=mean( results{isub}.perf(win2));
            end
            mean_result= mv_combine_results(results, 'average');
            %mean_auc(i,idev,:)=mean_result.perf{1, 1}  ;
            [h, pMMN(i,idev), ci, stats] = ttest(dataMMN, 0.5,'Tail', 'right');
            [h, pP3a(i,idev), ci, stats] = ttest(dataP3a, 0.5,'Tail', 'right');
            i=i+1;
            
        end
    end
end
save([filler1,filler2,filler3,filler4,method,'_alltime_mvpa_permutation_std3_ttest.mat'],'pMMN','pP3a');
% waveform
win1=16:18; win2=20:23;
for idev=1:2
    figure
    set(gcf, 'Position', [100, 100, 800, 200]);
    i=1;
    %%
    for i1stday=1:3
        for i2ndday=i1stday+1:4
            subplot(1,6,i)
            data=[];mean_result=[];
            load([method,'dev_',num2str(idev),num2str(i1stday),num2str(i2ndday),'_eachtwoday_MVPA_result_meeg_std_re3.mat']);
            results=result_meeg;
            
            mean_result= mv_combine_results(results, 'average');
            mean_auc=mean_result.perf{1, 1}  ;
            %%
            X=[];Y=[];upper=[];lower=[];
            y=mean_auc;
            err=(mean_result.perf_std{1, 1})/14;
            x=data4clssfy{md,subj, idev}.time;
            upper = y + err';
            lower = y - err';          
            X = [x, fliplr(x)];
            Y = [upper', fliplr(lower')];
            fill(X, Y, [166,189,219]/255, 'EdgeColor', 'none'); hold on;
            %%
            time=data4clssfy{md,subj, idev}.time;
            plot(time, mean_auc, 'LineWidth',1,'Color','b'); hold on;
            idx = (x >= 0.188) & (x <= 0.248);
            plot(time(idx), mean_auc(idx), 'LineWidth', 1.5, 'Color', 'k');hold on;
            idx = (x >= 0.288) & (x <= 0.348);
            plot(time(idx), mean_auc(idx), 'LineWidth', 1.5, 'Color', [202,0,32]./255);
            %}
            %%
            if pMMN(i,idev) < 0.05/6
                text(0.208, 0.535, '*', 'Color', 'k', 'FontSize', 12, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
            if pP3a(i,idev) < 0.05/6
                text(0.308, 0.535, '*', 'Color', [202,0,32]./255, 'FontSize', 12, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end            
            %%
            graphname = {'Day1 vs Day2',' Day1 vs Day3',' Day1 vs Day4','Day2 vs Day3',' Day2 vs Day4',' Day3 vs Day4'};
            title( graphname{i}, 'fontweight', 'bold');                    
            if i==1
                ylabel('AUC');
                xlabel('Time (s)');
            end
            
            i=i+1;
            set(gca, 'FontSize', 8, 'fontweight', 'bold');
            axis([-0.100 0.6 0.47 0.555]);
            xticks([0,0.2, 0.4, 0.6]);
            
            % 设置 y 轴刻度位置
            yticks([0.48, 0.50, 0.52,0.54,0.56,0.58,0.60]);ytickformat('%.2f');  % 格式化为保留两位小数
            set(gca, 'box', 'off');
            axis square;
            set(gca, 'TickLength', [0.03, 0.05]);
            hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.075, 0.5, 'Deviant', 'Units', 'normalized', 'Rotation', 90, ...
                'FontSize', 10,'fontweight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            change={'Small change','Large change'};
            annotation('textbox', [0.3, 0.85, 0.4, 0.05], ... 
                'String',change{idev}, ...
                'EdgeColor', 'none', ...  % 取消边框
                'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'FontSize', 10);
        end
    end
    ax = gca;
    set(ax, 'LooseInset', [0, 0, 0, 0]);
    filename = sprintf('MVPAdev%d.tiff', idev);
    print(filename, '-dtiff', '-r300');
end

%% temporal generalization单样本t检验 % 15-18 20-23
dataMMN=[];dataP3a=[];
win1=15:18; win2=20:23;
for idev=1:2
    i=1;
    for i1stday=1:3
        for i2ndday=i1stday+1:4
            data=[];mean_result=[];
            load([method,'dev_',num2str(idev),num2str(i1stday),num2str(i2ndday),'_eachtwoday_MVPA_result_meeg_std_re3.mat']);
            results=result_time;
            for isub=1:14
                dataMMN(isub)=mean(mean( results{isub}.perf(win1,win1)));
                dataP3a(isub)=mean(mean( results{isub}.perf(win2,win2)));
            end
            mean_result= mv_combine_results(results, 'average');
            %mean_auc(i,idev,:)=mean_result.perf{1, 1}  ;
            [h, pMMN(i,idev), ci, stats] = ttest(dataMMN, 0.5,'Tail', 'right');
            [h, pP3a(i,idev), ci, stats] = ttest(dataP3a, 0.5,'Tail', 'right');
            i=i+1;
            
        end
    end
end
save([filler1,filler2,filler3,filler4,method,'_alltime_mvpa_permutation_std3_timettest.mat'],'pMMN','pP3a');
% temporal generalization 表格
% 创建一个 5x8 的 cell 数组，用于存放表格数据
tableData = cell(5,8);

% 第一行第三至第八列：设置比较组标题
headers = {'Day1 vs Day2', 'Day1 vs Day3', 'Day1 vs Day4', 'Day2 vs Day3', 'Day2 vs Day4', 'Day3 vs Day4'};
tableData(1,3:8) = headers;

% 第一列：第2-3行为 'Small change'，第4-5行为 'Large change'
tableData{2,1} = 'Small change';
tableData{3,1} = 'Small change';
tableData{4,1} = 'Large change';
tableData{5,1} = 'Large change';

% 第二列：第2行和第4行为 'MMN'，第3行和第5行为 'P3a'
tableData{2,2} = 'MMN';
tableData{3,2} = 'P3a';
tableData{4,2} = 'MMN';
tableData{5,2} = 'P3a';

% 填充数据部分（第2-5行，第三至第八列），此处以随机数为例
data = [pMMN';pP3a'];  % 生成 4x6 的随机数据
temp = data(2,:);
data(2,:) = data(4,:);
data(4,:) = temp;

for i = 1:4
    for j = 1:6
        tableData{i+1, j+2} = data(i,j);
    end
end

% 将数据中小于 0.0083 的数值替换为星号 '*'
for i = 2:5         % 行：数据从第2行开始（不包含标题行）
    for j = 3:8     % 列：第三列到第八列为数据
        % 判断当前单元格是否为数值且小于 0.0083
        if isnumeric(tableData{i,j}) && tableData{i,j} < 0.0083
            tableData{i,j} = '*';
        else
            tableData{i,j} = [];

            
        end
    end
end

% 显示表格
disp(tableData);



