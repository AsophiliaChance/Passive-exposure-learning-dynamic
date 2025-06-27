clc;clear;

basedir = {'E:\passiveexp22222\data\4thanalysis\day1','E:\passiveexp22222\data\4thanalysis\day2',...
    'E:\passiveexp22222\data\4thanalysis\day3','E:\passiveexp22222\data\4thanalysis\day4'};
filt ='*_ref4.set';
filt1='*_nodelete1.set';
deviant_types = {'2','3'};
standard_type = '1';
art_thresh = 120; % +/- uV
twin = [-0.1 0.6];

%%%%%%%%%%%%%%
eeglab;
close(gcf);

%% %%%%%%%%
for md=1:4
    
    cd(basedir{md});
   
    outputdir = basedir{md};
    files = dir(filt);
    files1 = dir(filt1);
    ICtotal_brain_good=[];
    
    
    %%
    for curfile =1:length(files)
        file = files(curfile).name;
        file1 = files1(curfile).name;
        EEG = pop_loadset(file,pwd);

        EEG= pop_select( EEG, 'channel',{'C3','C4','Cz','F3','F4','P3','P4','Fz','Pz'});
        
        %%
        ind_boundary=[];
        for i=1:length(EEG.event)
            if ismember(str2num(EEG.event(i).type),[1,4,16,64])
                EEG.event(i).type='1';
            elseif ismember(str2num(EEG.event(i).type),[3,12,48,192])
                EEG.event(i).type='3';
            elseif ismember(str2num(EEG.event(i).type),[2,8,32,128])
                EEG.event(i).type='2';
            elseif ismember('boundary',EEG.event(i).type)
                ind_boundary=[ind_boundary;i];
            end
        end
        %%
        EEG0 = pop_loadset(file1,pwd);
        ind_boundary=[];
        for i=1:length(EEG0.event)
            if ismember(str2num(EEG0.event(i).type),[1,4,16,64])
                EEG0.event(i).type='1';
            elseif ismember(str2num(EEG0.event(i).type),[3,12,48,192])
                EEG0.event(i).type='3';
            elseif ismember(str2num(EEG0.event(i).type),[2,8,32,128])
                EEG0.event(i).type='2';
            elseif ismember('boundary',EEG0.event(i).type)
                ind_boundary=[ind_boundary;i];
            end
        end
        %% 去掉不必要的event
        tmpEEG0 = pop_selectevent( EEG0, 'type',{'2'},'deleteevents','on');
        
        %% epoch
        for idev=1
            uniqueToDEV =[];uniqueToSTD = [];indices1=[];indices2=[];indices3=[];indices4=[];
            EEG = pop_selectevent( EEG, 'type',{'1','2'},'deleteevents','on');
            EEG = pop_epoch( EEG, { deviant_types{idev}  }, [-1         0.6]);
            % 找到唯一值和它们的出现次数
            data = [EEG.event.epoch].';
            [unique_values, ~, indices] = unique(data);
            counts = histc(indices, 1:numel(unique_values));

            % 找到出现次数超过两次的值
            values_exceeding_two = unique_values(counts > 2);
              EEG = pop_selectevent( EEG, 'omitepoch', values_exceeding_two ,'deleteevents','off','deleteepochs','on','invertepochs','off');

            tmpEEG = pop_selectevent( EEG, 'type',{'2'},'deleteevents','on');
            data1 = [tmpEEG.event.urevent].';
            data2 = [tmpEEG0.event.urevent].';
            epochID{md,curfile} = cell2mat(arrayfun(@(x) find(data2 == x, 1), data1, 'UniformOutput', false));
            length(unique(epochID{md,curfile}(:) ))
            %             for ievent=1:length(urevent)
            %                 epochID{md,curfile}(ievent)=find(urevent(ievent)==urevent0);
            %             end
        end
    end
end

save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\epochID_small.mat','epochID')



%%
load('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\epochID_small.mat')
slide_data=[];
pari=[1,2;1,3;1,4;2,3;2,4;3,4];

for curfile=1:14
    %%
    for ipari=1:6
        imd=0;
    for md=pari(ipari,:)
        epochtosave=[];
        com=[];
        imd=imd+1;
        

        cd(basedir{md})
        files = dir(filt);
        %倒入数据并epoch改mark reref
        file = files(curfile).name;
        EEG = pop_loadset('filename',file,'filepath',basedir{md});

        EEG= pop_select( EEG, 'channel',{'C3','C4','Cz','F3','F4','P3','P4','Fz','Pz'});

        %%
        ind_boundary=[];
        for i=1:length(EEG.event)
            if ismember(str2num(EEG.event(i).type),[1,4,16,64])
                EEG.event(i).type='1';
            elseif ismember(str2num(EEG.event(i).type),[3,12,48,192])
                EEG.event(i).type='3';
            elseif ismember(str2num(EEG.event(i).type),[2,8,32,128])
                EEG.event(i).type='2';
            elseif ismember('boundary',EEG.event(i).type)
                ind_boundary=[ind_boundary;i];
            end
        end
                EEG = pop_selectevent( EEG, 'type',{'1','2'},'deleteevents','on');
        %% epoch
        idev=1;
        
        com = intersect(epochID{pari(ipari,1),curfile}, epochID{pari(ipari,2),curfile});   
        tmp=epochID{md,curfile};
        data1=com;data2=tmp;
        epochtosave = cell2mat(arrayfun(@(x) find(data2 == x, 1), data1, 'UniformOutput', false));
        
        
        EEG = pop_epoch( EEG, { deviant_types{idev}}, [-1         0.60]);
                    % 找到唯一值和它们的出现次数
            data = [EEG.event.epoch].';
            [unique_values, ~, indices] = unique(data);
            counts = histc(indices, 1:numel(unique_values));

            % 找到出现次数超过两次的值
            values_exceeding_two = unique_values(counts > 2);
               EEG = pop_selectevent( EEG, 'omitepoch', values_exceeding_two ,'deleteevents','off','deleteepochs','on','invertepochs','off');


       [ EEG,indices0] = pop_selectevent( EEG, 'epoch',epochtosave' ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        
        [DEV1,indices1] = pop_epoch( EEG, { deviant_types{idev} }, [-0.1         0.60]);
        [STD1,indices2] = pop_epoch( EEG, { '1' }, [-0.1         0.60]);
         A_only = setdiff(indices1, indices2);
         B_only = setdiff(indices2, indices1);
        
        DEV = pop_rmbase( DEV1, [-100 0] ,[]);
        STD = pop_rmbase( STD1, [-100 0] ,[]);
        
        Diff=DEV;
        Diff.data=DEV.data-STD.data;
        %plotData=Diff.data;%250-400ms
        P3aData=squeeze(mean(Diff.data(:,88:113,:),2));%250-350ms
         MMNData=squeeze(mean(Diff.data(:,64:88,:),2));%150-250ms
        slide_data_tmp=[];
        for x=1:Diff.trials-99
            slide_P3a_tmp(:,x)=mean(P3aData(:,x:x+99),2);
            slide_MMN_tmp(:,x)=mean(MMNData(:,x:x+99),2);
        end
        slide_MMN{ipari,imd,curfile}=slide_MMN_tmp;
        slide_P3a{ipari,imd,curfile}=slide_P3a_tmp;
        
    end       
    end
end
chanlocs=EEG.chanlocs;
time=EEG.times;
save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intraday_slidedata_small.mat','slide_MMN','slide_P3a','chanlocs','time')

%% PLOT within-sub correlation

load('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intraday_slidedata_small.mat')
cof_XC=[];cof_lag=[]; data1=[]; data2=[];is_significant =[];p=[];non_significant=[];
for P3a=1:2
    if P3a==1
        slide_data=slide_P3a;
    else
        slide_data=slide_MMN;
    end
for ichan=1:9%{'C3';'C4';'Cz';'F3';'F4';'P3';'P4';'Fz';'Pz'}
for ipari=1:6
%     figure
for curfile=1:14
%     subplot(4,4,curfile)
    for imd=1:2
%         plot(1:length(slide_data{ipari,imd,curfile}(1,:)),slide_data{ipari,imd,curfile}(1,:));hold on
        
    end

    data1=slide_data{ipari,1,curfile}(ichan,:);data2=slide_data{ipari,2,curfile}(ichan,:);
    [xc, lags] = xcorr(data1, data2, 'coeff');%cof1(ipari,curfile) = corr(data1', data2'); 
    [cof_XC{ichan,P3a}(ipari,curfile), idx] = max(xc);
    cof_lag{ichan,P3a}(ipari,curfile)  = lags(idx);
    % 定义显著性检验函数


% 显著性水平
alpha = 0.05;
r=cof_XC{ichan,P3a}(ipari,curfile);
n=length(data1);
% 初始化显著性结果存储
    t = r .* sqrt((n - 2) ./ (1 - r.^2));
    p{ichan,P3a}(ipari,curfile) = 2 * (1 - tcdf(abs(t), n - 2));

is_significant {ichan,P3a}(ipari,curfile)= p{ichan,P3a}(ipari,curfile) < alpha;

end
end
non_significant(ichan,P3a)= 6*14-sum(sum(is_significant{ichan,P3a}));
end
end
within_pvalue=p;
save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intraday_xcorr_small.mat','cof_XC','is_significant','non_significant')

%% 比较六组intraday_xcorr是否有统计学差异，采用非参数方差分析
% stats=[];
load('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intraday_xcorr_small.mat')
for P3a=2
for ichan=1:9
data =cof_XC{ichan,P3a};
data = data(:);
mu0 =0.5; % 总体均值

% 执行单样本t检验
%[h, p, ci, stats] = ttest(data, mu0,'Tail', 'right');

% 创建分组变量
group = [repmat({'A'}, 14, 1); ...
         repmat({'B'}, 14, 1); ...
         repmat({'C'}, 14, 1); ...
         repmat({'D'}, 14, 1); ...
         repmat({'E'}, 14, 1); ...
         repmat({'F'}, 14, 1)];

% 进行 Kruskal-Wallis 检验
 [p, tbl, stats{ichan,P3a}] = kruskalwallis(data, group);

% 显示结果
%disp(['Kruskal-Wallis 检验 p 值: ', num2str(p)]);

% 如果检验结果显著，进行多重比较
 if p < 0.05
%     disp('检验结果显著，进行多重比较...');

  results_bonferroni {ichan,P3a}= multcompare(stats{ichan,P3a}, 'CType', 'bonferroni','display','on');
else
    disp('检验结果不显著，组间没有显著差异。');
end
end
end
save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intraday_statistics_small.mat','results_bonferroni','stats','chanlocs','time')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intrasubject  1
pari=[];

slide_data=[];
imd=0;
for curfile1=1:13
    
    for curfile2=curfile1+1:14
        imd=imd+1;
        pari(imd,1)=curfile1;
        pari(imd,2)=curfile2;
    end
end
%% intrasub analysis 2
slide_MMN=[];slide_P3a=[];
for md=1:4
    com=[];        
    cd(basedir{md})
    files = dir(filt);
   
    %%
    for ipari=1:size(pari,1)
        imd=0;
    for curfile=pari(ipari,:)
        imd=imd+1; % 1/2
        epochtosave=[];
        com=[];
        %倒入数据并epoch改mark reref
        file = files(curfile).name;
        EEG = pop_loadset('filename',file,'filepath',basedir{md});

        EEG= pop_select( EEG, 'channel',{'C3','C4','Cz','F3','F4','P3','P4','Fz','Pz'});

        %%
        ind_boundary=[];
        for i=1:length(EEG.event)
            if ismember(str2num(EEG.event(i).type),[1,4,16,64])
                EEG.event(i).type='1';
            elseif ismember(str2num(EEG.event(i).type),[3,12,48,192])
                EEG.event(i).type='3';
            elseif ismember(str2num(EEG.event(i).type),[2,8,32,128])
                EEG.event(i).type='2';
            elseif ismember('boundary',EEG.event(i).type)
                ind_boundary=[ind_boundary;i];
            end
        end
                EEG = pop_selectevent( EEG, 'type',{'1','2'},'deleteevents','on');
        %% epoch
        idev=1;
        
        com = intersect(epochID{md,pari(ipari,1)}, epochID{md,pari(ipari,2)});   
        tmp=epochID{md,curfile};
        data1=com;data2=tmp;
        epochtosave = cell2mat(arrayfun(@(x) find(data2 == x, 1), data1, 'UniformOutput', false));
        
        
        EEG = pop_epoch( EEG, { deviant_types{idev}}, [-1         0.60]);
                    % 找到唯一值和它们的出现次数
            data = [EEG.event.epoch].';
            [unique_values, ~, indices] = unique(data);
            counts = histc(indices, 1:numel(unique_values));

            % 找到出现次数超过两次的值
            values_exceeding_two = unique_values(counts > 2);
              EEG = pop_selectevent( EEG, 'omitepoch', values_exceeding_two ,'deleteevents','off','deleteepochs','on','invertepochs','off');


       [ EEG,indices0] = pop_selectevent( EEG, 'epoch',epochtosave' ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        
        [DEV1,indices1] = pop_epoch( EEG, { deviant_types{idev} }, [-0.1         0.60]);
        [STD1,indices2] = pop_epoch( EEG, { '1' }, [-0.1         0.60]);
         A_only = setdiff(indices1, indices2);
         B_only = setdiff(indices2, indices1);
        
        DEV = pop_rmbase( DEV1, [-100 0] ,[]);
        STD = pop_rmbase( STD1, [-100 0] ,[]);
        
        Diff=DEV;
        Diff.data=DEV.data-STD.data;
         P3aData=squeeze(mean(Diff.data(:,88:113,:),2));%250-350ms
         MMNData=squeeze(mean(Diff.data(:,64:88,:),2));%150-250ms
        slide_P3a_tmp=[];slide_MMN_tmp=[];
        for x=1:Diff.trials-99
            slide_P3a_tmp(:,x)=mean(P3aData(:,x:x+99),2);
            slide_MMN_tmp(:,x)=mean(MMNData(:,x:x+99),2);
        end
        slide_MMN{md,ipari,imd}=slide_MMN_tmp;
        slide_P3a{md,ipari,imd}=slide_P3a_tmp;
       
        
    end       
    end
end

chanlocs=EEG.chanlocs;
time=EEG.times;
intrasub_MMN=slide_MMN;
intrasub_P3a=slide_P3a;
save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intrasub_slidedata_small.mat','intrasub_MMN','intrasub_P3a','chanlocs','time')


%% intrasub statistics 3
load('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intrasub_slidedata_small.mat')
cof_XC=[];cof_lag=[];
is_significant=[];non_significant=[];
imd=0;
bad_pari=[];
cof_XC=[];cof_lag=[]; data1=[]; data2=[];is_significant =[];
for P3a=1:2
    if P3a==1
        slide_data_1=intrasub_P3a;
    else
        slide_data_1=intrasub_MMN;
    end
    for ichan=1:9
for md=1:4
for ipari=1:size(pari,1)

    data1=slide_data_1{md,ipari,1}(ichan,:);data2=slide_data_1{md,ipari,2}(ichan,:);
    if length(data1)==length(data2)
    [xc, lags] = xcorr(data1, data2, 'coeff');%cof1(ipari,curfile) = corr(data1', data2'); 
    [cof_XC{ichan,P3a}(md,ipari), idx] = max(xc);
    cof_lag{ichan,P3a}(md,ipari)  = lags(idx);
    % 定义显著性检验函数


% 显著性水平
alpha = 0.05;
r=cof_XC{ichan,P3a}(md,ipari);
n=length(data1);
% 初始化显著性结果存储
    t = r .* sqrt((n - 2) ./ (1 - r.^2));
   intrasubject_pvalue {ichan}(md,ipari) = 2 * (1 - tcdf(abs(t), n - 2));

is_significant {ichan,P3a}(md,ipari)= intrasubject_pvalue {ichan}(md,ipari) < alpha;

    else 
        imd=imd+1;
        bad_pari(imd,:)=pari(ipari,:);
        bad_slidata(imd,:)=slide_data_1(md,ipari,:);
    end

end
end
non_significant(ichan,P3a)= 4*196-sum(sum(is_significant{ichan,P3a}));

    end
end
a=1-sum(non_significant(:,1))/(9*4*196);b=1-sum(non_significant(:,2))/(9*4*196);
save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intrasub_xcorr_small.mat','cof_XC','non_significant','is_significant','chanlocs','time')

%% 比较四天intrasub_xcorr是否有统计学差异，采用非参数方差分析
stats=[];data=[];
for ichan=1:9
    for P3a=1:2
data =cof_XC{ichan,P3a};
data = data(:);
mu0 =0.5; % 总体均值

% 执行单样本t检验
%[h, p, ci, stats] = ttest(data, mu0,'Tail', 'right');

% 创建分组变量
group = [repmat({'A'}, 91, 1); ...
         repmat({'B'}, 91, 1); ...
         repmat({'C'}, 91, 1); ...
         repmat({'D'}, 91, 1)];

% 进行 Kruskal-Wallis 检验
[p, tbl, stats{ichan,P3a}] = kruskalwallis(data, group);
    end
end
%% 多重比较
    for P3a=2
for ichan=1:9

figure
   results_bonferroni {ichan,P3a}= multcompare(stats{ichan,P3a}, 'CType', 'bonferroni');

    end
end
%save('E:\passiveexp22222\data\4thanalysis\MVPA_dataNresult\intrasub_statistics_small.mat','results_bonferroni','stats')


