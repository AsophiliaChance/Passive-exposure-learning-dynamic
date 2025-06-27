clc;clear
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');

filename={'day1_nsubavg_re.mat','day2_nsubavg_re.mat','day3_nsubavg_re.mat','day4_nsubavg_re.mat'};
%%
channelsToAvg = {'F3', 'F4', 'Fz'};

for md = 1:4
    load(filename{md});channel_idx = find(ismember( Diff_avg{1,1} .label, channelsToAvg));
    % 对每个受试者和设备进行数据重排及平均
    for isub = 1:14
        for idev = 1:2
            % 对 Diff_avg 中 F3, F4, Fz 通道数据求平均
            avg_diff=Diff_avg{isub,idev};
            avg_diff.avg = mean(Diff_avg{isub,idev}.avg(channel_idx,:));
            Diff{md,isub,idev} = avg_diff;
            
            % 同理，对 STD_avg 与 DEV_avg 也可以计算平均值
            avg_std=STD_avg{isub,idev};
            avg_std.avg  = mean(STD_avg{isub,idev}.avg (channel_idx,:));
            STD{md,isub,idev} = avg_std;
            
            avg_dev=DEV_avg{isub,idev};
            avg_dev.avg  = mean(DEV_avg{isub,idev}.avg (channel_idx,:));
            DEV{md,isub,idev}  = avg_dev;
        end
    end
end


%% grand average for plotting
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
for md=1:4
    for idev=1:2
        STD_gavg {md,idev}      = ft_timelockgrandaverage(cfg, STD{md,:,idev});
        DEV_gavg {md,idev}      = ft_timelockgrandaverage(cfg, DEV{md,:,idev});
        Diff_gavg {md,idev}      = ft_timelockgrandaverage(cfg, Diff{md,:,idev});
    end
    
end

% channelsToAvg = {'F3', 'F4', 'Fz'};
% chanIdx = find(ismember( Diff_gavg{1,1}  .label, channelsToAvg));
%% define window for 2 deviant
minTime=[];maxTime=[];
for idev = 1:2
    data=[];
for md = 1:4
    
        data(md,:) = Diff_gavg{md, idev}.avg;        

                

end
data=mean(data,1);
time=Diff_gavg{md, idev}.time;
idx_min = find(Diff_gavg{md, idev}.time >= 0.19 & Diff_gavg{md, idev}.time <= 0.26);
        [minVal, minIdx] = min(data(idx_min));
        minTime(idev) = time(idx_min(minIdx));
        idx_max = find(time >= 0.29 & time <= 0.34);
        [maxVal, maxIdx] = max(data(idx_max));
        maxTime(idev) = time(idx_max(maxIdx));
end

win=[];
win.mmn=minTime;win.p3=maxTime;
 %% define window jari
% win.mmn1=0.19; win.mmn2=0.24;
% win.p31=0.25; win.p32=0.3;
%%
%% define window combine deviant
winlenght=[0.025,0.03,0.035,0.04,0.05];
clc
diary('output.txt');
for iwin=3%:length(winlenght)
Diff_ga=[];
for idev=1:2   
    for md=1:4
    Diff_ga (idev,md,:)      =  Diff_gavg {md,idev}.avg;
    end
end

 Diff_ga =squeeze(mean(mean(Diff_ga,2),1)); 

time=Diff_gavg{md, idev}.time;
idx_min = find(Diff_gavg{md, idev}.time >= 0.19 & Diff_gavg{md, idev}.time <= 0.26);
        [minVal, minIdx] = min(Diff_ga(idx_min));
        minTime = time(idx_min(minIdx));
        idx_max = find(time >= 0.25 & time <= 0.35);
        [maxVal, maxIdx] = max(Diff_ga(idx_max));
        maxTime= time(idx_max(maxIdx));
        win=[];
win.mmn1=minTime-winlenght(iwin); win.mmn2=minTime+winlenght(iwin);
win.p31=maxTime-winlenght(iwin); win.p32=maxTime+winlenght(iwin);
%% Peak detection
latency=[]; amplitude=[];
for md = 1:4
    for idev = 1:2
        
        for isub = 1:size(Diff,2)
            data = Diff{md, isub, idev};
            idx = find(data.time >= win.mmn1 & data.time <= win.mmn2);
%             amplitude_diff.mmn(md, isub, idev) = mean(data.avg(idx));
            [~, peakIdx] = min(data.avg(idx));
            latency_diff.mmn(md, isub, idev) = data.time(idx(peakIdx));
             idx = find(data.time >= data.time(idx(peakIdx))-0.025 & data.time <= data.time(idx(peakIdx))+0.025 );
            amplitude_diff.mmn(md, isub, idev) = mean(data.avg(idx));
            
            
            idx = find(data.time >= win.p31 & data.time <= win.p32);
           % amplitude_diff.p3(md, isub, idev) = mean(data.avg(idx));
            [~, peakIdx] = max(data.avg(idx));
            latency_diff.p3(md, isub, idev) = data.time(idx(peakIdx));
             idx = find(data.time >= data.time(idx(peakIdx))-0.025 & data.time <= data.time(idx(peakIdx))+0.025 );
            amplitude_diff.p3(md, isub, idev) = mean(data.avg(idx));
        end
    end
end
%%
% 假设你已有这两个变量：amplitude_diff 和 latency_diff
% 其中每个字段是 4×14×2 数组（day × subject × condition）

% 单独提取字段
amplitude_mmn = amplitude_diff.mmn;   % 4×14×2
amplitude_p3  = amplitude_diff.p3;

latency_mmn = latency_diff.mmn;
latency_p3  = latency_diff.p3;

% 保存成独立变量，R中易读取（避免嵌套结构体）
save('erp_statisticsdata_simple.mat', ...
     'amplitude_mmn', 'amplitude_p3', ...
     'latency_mmn', 'latency_p3', '-v7');

save('erp_statisticsdata.mat','amplitude_diff','latency_diff');
%%
%% 1. 准备数据（模拟）
% 对应的Day
set(groot,'defaultAxesTickLength',[0.03 0.06])   % 放在脚本最开头

clc
data = cat(3, amplitude_diff.mmn, amplitude_diff.p3);
figure('Position', [100 100 800 230])    % 改宽一点，合理布局
componentTitles = {'Small MMN','Large MMN','Small P3a','Large P3a'};

for i = 1:4
    subplot(1,4,i)
    hold on

    Days = repmat((1:4)',14,1);     
    Subjects = repelem((1:14)',4);
    Y = reshape(data(:,:,i),[],1); 
    n = 14;    % 每组14个受试者，注意加！

    T = table(Subjects,Days,Y);
    T.Subjects = categorical(T.Subjects);
    T.Days = categorical(T.Days);

    lme = fitlme(T, 'Y ~ Days + (Days|Subjects)');
    %gee = fitgee(T, 'Y ~ Days', 'Subjects', 'Correlation','AR(1)');

    %disp(lme)

    if iscategorical(T.Days)
        T.DayNum = double(T.Days);
    else
        T.DayNum = T.Days;
    end

    days_list = unique(T.DayNum);
    mean_Y = groupsummary(T,"DayNum","mean","Y");
    sem_Y  = groupsummary(T,"DayNum","std","Y");
    sem_Y.std_Y  = sem_Y.std_Y / sqrt(n);

    % ---- 只画柱形图 + errorbar ----
    bar(days_list, mean_Y.mean_Y, 0.5, 'FaceAlpha',0.5, 'EdgeColor','none');
    errorbar(days_list, mean_Y.mean_Y, sem_Y.std_Y, 'k.', 'LineWidth',1, 'CapSize',5);

ax = gca;                        % 当前 subplot 句柄

if i <= 2                       % -------- 前两幅（MMN）--------
       % 负值朝上
    ylim(ax,[-3 0])             % 0 到 –2 ?V
    ax.YTick    = [-3 -2 -1 0]; 
    ax.YDir     = 'reverse'; % 可选：刻度
else                            % -------- 后两幅（P3a）--------
    ax.YDir     = 'normal';     % 正值朝上（默认）
    ylim(ax,[0 3])              % 0 到 3 ?V
    ax.YTick    = 0:1:3;        % 可选：刻度
end
    

    xticks(days_list)
   % xlim([0.5 4.5])
   % grid on   % 打开网格更专业
set(gca, 'FontSize', 10);
if i==1
xlabel('Day', 'FontSize', 15,'FontWeight', 'bold');
ylabel('Amplitude (μV)', 'FontSize', 15,'FontWeight', 'bold');
end
title(componentTitles{i}, 'FontSize', 15)
set(findall(gcf,'-property','FontWeight'), 'FontWeight', 'bold')
% --- SIGNIFICANCE: Day 2–4 vs Day 1 -----------------------------
 % --- SIGNIFICANCE: Day 2–4 vs Day 1 -----------------------------
    [~,~,FE] = fixedEffects(lme,'DFMethod','Satterthwaite');
   amplitudeFE{i}= FE;
   disp(FE)
    pvals = FE.pValue;                    % (Intercept, Days_2, Days_3, Days_4)

    p2star = @(p) repmat('*',1,(p<0.05)+(p<0.01)+(p<0.001));

    % Determine offsets -------------------------------------------------------
    ySpan        = range(ylim(ax));           % full axis span
    lineOffset   = 0.15 * ySpan;             % distance from bar top to line (was 0.08)
    starOffset   = 0.02 * ySpan;             % distance star above line for clarity
    if strcmp(ax.YDir,'reverse')             % reverse axis => invert offsets
        lineOffset = -lineOffset;
        starOffset = -starOffset;
    end

    nStar = 0;   % counter to stack multiple lines vertically if needed

    for d = 2:4
        if pvals(d) < 0.05
            nStar = nStar + 1;
            % Height of this significance line (stack if multiple) ---
            barTops = mean_Y.mean_Y([1 d]);
            yBase   = max(barTops);
            if strcmp(ax.YDir,'reverse'), yBase = min(barTops)-1; end
            yLine   = yBase + (nStar-1)*0.15+0.45;    % higher line level

            % Draw horizontal line connecting Day1 and Day d ----------
            line([1 d],[yLine yLine], 'Parent',ax, ...
                 'Color','k','LineWidth',1.2');

            % Place star(s) slightly above line ----------------------
            xMid = mean([1 d]);
            yStar = yLine + starOffset;            % closer to the line
            text(xMid, yStar, p2star(pvals(d)), ...
                 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle', ...
                 'FontSize',12);
        end
    end


    hold off
end
sgtitle('Amplitude', 'FontWeight','bold', 'FontSize',20);
print(gcf,'-dtiff','-r400','bar_amplitude.tif');  % saves 400?dpi TIFF in current folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 准备数据（模拟）
% 对应的Day
set(groot,'defaultAxesTickLength',[0.03 0.06])   % 放在脚本最开头


data = cat(3, latency_diff.mmn, latency_diff.p3);
figure('Position', [100 100 800 230])    % 改宽一点，合理布局
componentTitles = {'Small MMN','Large MMN','Small P3a','Large P3a'};

for i = 1:4
    subplot(1,4,i)
    hold on

    Days = repmat((1:4)',14,1);     
    Subjects = repelem((1:14)',4);
    Y = reshape(data(:,:,i),[],1); 
    n = 14;    % 每组14个受试者，注意加！

    T = table(Subjects,Days,Y);
    T.Subjects = categorical(T.Subjects);
    T.Days = categorical(T.Days);

    lme = fitlme(T, 'Y ~ Days + (Days|Subjects)');
    %disp(lme)

    if iscategorical(T.Days)
        T.DayNum = double(T.Days);
    else
        T.DayNum = T.Days;
    end

    days_list = unique(T.DayNum);
    mean_Y = groupsummary(T,"DayNum","mean","Y");
    sem_Y  = groupsummary(T,"DayNum","std","Y");
    sem_Y.std_Y  = sem_Y.std_Y / sqrt(n);

    % ---- 只画柱形图 + errorbar ----
    bar(days_list, 1000*(mean_Y.mean_Y), 0.5, 'FaceAlpha',0.5, 'EdgeColor','none');
    errorbar(days_list, 1000*(mean_Y.mean_Y), 1000*(sem_Y.std_Y), 'k.', 'LineWidth',1, 'CapSize',5);

ax = gca;                        % 当前 subplot 句柄

if i <= 2                       % -------- 前两幅（MMN）--------
       % 负值朝上
    ylim(ax,[190 260])             % 0 到 –2 ?V
    ax.YTick    = 190:20:260; 
    %ax.YDir     = 'reverse'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else                            % -------- 后两幅（P3a）--------
    ax.YDir     = 'normal';     % 正值朝上（默认）
    ylim(ax,[280 350])              % 0 到 3 ?V
    ax.YTick    = 280:20:350;        % 可选：刻度
end
%     title(componentTitles{i})

    xticks(days_list)
   % xlim([0.5 4.5])
   % grid on   % 打开网格更专业
set(gca, 'FontSize', 10);
if i==1
    xlabel('Day', 'FontSize', 15,'FontWeight', 'bold');
    ylabel('Latency (ms)', 'FontSize', 15,'FontWeight', 'bold');
end
title(componentTitles{i}, 'FontSize', 15)
set(findall(gcf,'-property','FontWeight'), 'FontWeight', 'bold')
% --- SIGNIFICANCE: Day 2–4 vs Day 1 -----------------------------
 % --- SIGNIFICANCE: Day 2–4 vs Day 1 -----------------------------
    [~,~,FE] = fixedEffects(lme,'DFMethod','Satterthwaite');
   latencyFE{i}= FE;
   disp(FE)
    pvals = FE.pValue;                    % (Intercept, Days_2, Days_3, Days_4)

    p2star = @(p) repmat('*',1,(p<0.05)+(p<0.01)+(p<0.001));

    % Determine offsets -------------------------------------------------------
    ySpan        = range(ylim(ax));           % full axis span
    lineOffset   = 0.015 * ySpan;             % distance from bar top to line (was 0.08)
    starOffset   = 0.02 * ySpan;             % distance star above line for clarity
    if strcmp(ax.YDir,'reverse')             % reverse axis => invert offsets
        lineOffset = -lineOffset;
        starOffset = -starOffset;
    end

    nStar = 0;   % counter to stack multiple lines vertically if needed

    for d = 2:4
        if pvals(d) < 0.05
            nStar = nStar + 1;
            % Height of this significance line (stack if multiple) ---
            barTops = mean_Y.mean_Y([1 d]);
            yBase   = max(barTops);
            if strcmp(ax.YDir,'reverse'), yBase = min(barTops)-1; end
            yLine   = (yBase+ (nStar-1)*0.1+0.01)*1000;    % higher line level

            % Draw horizontal line connecting Day1 and Day d ----------
            line([1 d],[yLine yLine], 'Parent',ax, ...
                 'Color','k','LineWidth',1.2');

            % Place star(s) slightly above line ----------------------
            xMid = mean([1 d]);
            yStar = yLine + starOffset;            % closer to the line
            text(xMid, yStar, p2star(pvals(d)), ...
                 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle', ...
                 'FontSize',12);
        end
    end


    hold off
end
sgtitle('Latency', 'FontWeight','bold', 'FontSize',20);
print(gcf,'-dtiff','-r400','bar_latency.tif');  % saves 400?dpi TIFF in current folder


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for md = 1:4
    for idev = 1:2
        
        for isub = 1:size(Diff,2)
            data = STD{md, isub, idev};
            idx = find(data.time >= win.mmn1 & data.time <= win.mmn2);
            %amplitude_std.mmn(md, isub, idev) = mean(data.avg(idx));
            [~, peakIdx] = min(data.avg(idx));
            latency_std.mmn(md, isub, idev) = data.time(idx(peakIdx));
                         idx = find(data.time >= data.time(idx(peakIdx))-0.025 & data.time <= data.time(idx(peakIdx))+0.025 );
            amplitude_std.mmn(md, isub, idev) = mean(data.avg(idx));
            
            idx = find(data.time >= win.p31 & data.time <= win.p32);
            %amplitude_std.p3(md, isub, idev) = mean(data.avg(idx));
            [~, peakIdx] = max(data.avg(idx));
            latency_std.p3(md, isub, idev) = data.time(idx(peakIdx));
            idx = find(data.time >= data.time(idx(peakIdx))-0.025 & data.time <= data.time(idx(peakIdx))+0.025 );
            amplitude_std.p3(md, isub, idev) = mean(data.avg(idx));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for md = 1:4
    for idev = 1:2
        
        for isub = 1:size(Diff,2)
            data = DEV{md, isub, idev};
            idx = find(data.time >= win.mmn1 & data.time <= win.mmn2);
            %amplitude_dev.mmn(md, isub, idev) = mean(data.avg(idx));
            [~, peakIdx] = min(data.avg(idx));
            latency_dev.mmn(md, isub, idev) = data.time(idx(peakIdx));
             idx = find(data.time >= data.time(idx(peakIdx))-0.025 & data.time <= data.time(idx(peakIdx))+0.025 );
            amplitude_dev.mmn(md, isub, idev) = mean(data.avg(idx));           
            
            idx = find(data.time >= win.p31 & data.time <= win.p32);
            %amplitude_dev.p3(md, isub, idev) = mean(data.avg(idx));
            [~, peakIdx] = max(data.avg(idx));
            latency_dev.p3(md, isub, idev) = data.time(idx(peakIdx));
             idx = find(data.time >= data.time(idx(peakIdx))-0.025 & data.time <= data.time(idx(peakIdx))+0.025 );
            amplitude_dev.p3(md, isub, idev) = mean(data.avg(idx)); 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latency_stim1.mmn=cat(3, latency_std.mmn(:, :, 1), latency_dev.mmn(:, :,1));
latency_stim1.p3=cat(3, latency_std.p3(:, :, 1), latency_dev.p3(:, :,1));
latency_stim2.mmn=cat(3, latency_std.mmn(:, :, 2), latency_dev.mmn(:, :,2));
latency_stim2.p3=cat(3, latency_std.p3(:, :, 2), latency_dev.p3(:, :,2));
amplitude_stim1.mmn=cat(3, amplitude_std.mmn(:, :, 1), amplitude_dev.mmn(:, :,1));
amplitude_stim1.p3=cat(3, amplitude_std.p3(:, :, 1), amplitude_dev.p3(:, :,1));
amplitude_stim2.mmn=cat(3, amplitude_std.mmn(:, :, 2), amplitude_dev.mmn(:, :,2));
amplitude_stim2.p3=cat(3, amplitude_std.p3(:, :, 2), amplitude_dev.p3(:, :,2));
%% anova
%% MATLAB 脚本：对4×14×2×6数据进行重复测量方差分析


data1{1}=latency_diff.mmn;
data1{2}=latency_diff.p3;
data1{3}=latency_dev.mmn;
data1{4}=latency_dev.p3;
data1{5}=amplitude_diff.mmn;
data1{6}=amplitude_diff.p3;
data1{7}=amplitude_dev.mmn;
data1{8}=amplitude_dev.p3;
data1{9}=latency_stim1.mmn;
data1{10}=latency_stim1.p3;
data1{11}=latency_stim2.mmn;
data1{12}=latency_stim2.p3;
data1{13}=amplitude_stim1.mmn;
data1{14}=amplitude_stim1.p3;
data1{15}=amplitude_stim2.mmn;
data1{16}=amplitude_stim2.p3;
datanam{1}  = 'Latency difference wave MMN';
datanam{2}  = 'Latency difference wave P3';
datanam{3}  = 'Latency deviant MMN';
datanam{4}  = 'Latency deviant P3';
datanam{5}  = 'Amplitude difference wave MMN';
datanam{6}  = 'Amplitude difference wave P3';
datanam{7}  = 'Amplitude deviant MMN';
datanam{8}  = 'Amplitude deviant P3';
datanam{9}  = 'Latency small stimulus MMN';
datanam{10} = 'Latency small stimulus P3';
datanam{11} = 'Latency large stimulus MMN';
datanam{12} = 'Latency large stimulus P3';
datanam{13} = 'Amplitude small stimulus MMN';
datanam{14} = 'Amplitude small stimulus P3';
datanam{15} = 'Amplitude large stimulus MMN';
datanam{16} = 'Amplitude large stimulus P3';

%%


for idata =  [1,2,5,6] %1:16%
    data = data1{idata};
    fprintf('%s\n', repmat('*', 1, 100));

fprintf('%s\n', datanam{idata});
fprintf('window length: %.3f\n', winlenght(iwin));

  

    % Data reshaping
    data_2D = reshape(data, 14, []);  % Resulting matrix is 14×48
%{
    %% Normality Test (Lilliefors Test)
    fprintf('Normality Test (Lilliefors Test):\n');
    for i = 1:size(data_2D,2)
        [h, p] = lillietest(data_2D(:,i));
        if h == 0
            %fprintf('Column %d: Failed to reject the null hypothesis (p = %.4f)\n', i, p);
        else
           fprintf('%s: Rejected the null hypothesis (p = %.4f)\n', varNames{i}, p);
        end
    end

    % Levene Test for Homogeneity of Variance (small difference wave)
    [p_levene, stats_levene] = vartestn(data_2D(:,1:4), [], 'TestType', 'LeveneAbsolute', 'Display', 'off');
    fprintf('\nLevene Test for Homogeneity of Variance (small difference wave):\n');
    fprintf('Test Statistic: %.4f, p-value: %.4f\n', stats_levene.fstat, p_levene);

    % Levene Test for Homogeneity of Variance (Columns 5-8)
    [p_levene, stats_levene] = vartestn(data_2D(:,5:8), [], 'TestType', 'LeveneAbsolute', 'Display', 'off');
    fprintf('\nLevene Test for Homogeneity of Variance (large difference wave):\n');
    fprintf('Test Statistic: %.4f, p-value: %.4f\n', stats_levene.fstat, p_levene);

    % Mauchly's Test for Sphericity (small difference wave)
    S = cov(data_2D(:,1:4));
    p_var = size(S,1);
    n = size(data_2D(:,1:4),1);
    W = det(S) / ((trace(S) / p_var)^p_var);
    C = 1 - (2*p_var^2 + p_var + 2) / (6*(p_var - 1)*p_var);
    chi2_stat = - (n - 1) * C * log(W);
    df = (p_var*(p_var - 1))/2 - 1;
    p_value_sphericity = 1 - chi2cdf(chi2_stat, df);
    fprintf('\nMauchly''s Test for Sphericity (small difference wave):\n');
    fprintf('W = %.4f\n', W);
    fprintf('Chi-square Statistic: %.4f, Degrees of Freedom: %d, p-value: %.4f\n', chi2_stat, df, p_value_sphericity);

    % Mauchly's Test for Sphericity (Columns 5-8)
    S = cov(data_2D(:,5:8));
    p_var = size(S,1);
    n = size(data_2D(:,1:4),1); % Using the same n as above
    W = det(S) / ((trace(S) / p_var)^p_var);
    C = 1 - (2*p_var^2 + p_var + 2) / (6*(p_var - 1)*p_var);
    chi2_stat = - (n - 1) * C * log(W);
    df = (p_var*(p_var - 1))/2 - 1;
    p_value_sphericity = 1 - chi2cdf(chi2_stat, df);
    fprintf('\nMauchly''s Test for Sphericity (large difference wave):\n');
    fprintf('W = %.4f\n', W);
    fprintf('Chi-square Statistic: %.4f, Degrees of Freedom: %d, p-value: %.4f\n', chi2_stat, df, p_value_sphericity);

    %% Friedman Test and Pairwise Comparisons (Holm Correction) for small difference wave
    dataa = data_2D(:,1:4);
    [p_friedman, tbl, stats] = friedman(dataa, 1, 'off');
    disp('Friedman Test Result (small difference wave):');
    disp(tbl);

    k = size(dataa, 2);
    pvals = nan(k, k);
    for i = 1:k-1
        for j = i+1:k
            pvals(i, j) = signrank(dataa(:, i), dataa(:, j));
            pvals(j, i) = pvals(i, j);
        end
    end

    idx = find(triu(ones(k), 1));
    p_vector = pvals(idx);
    m = length(p_vector);
    [sorted_p, sort_index] = sort(p_vector);
    adjusted_p = zeros(size(sorted_p));
    for i = 1:m
        adjusted_p(i) = min(sorted_p(i) * (m - i + 1), 1);
    end
    for i = 2:m
        adjusted_p(i) = max(adjusted_p(i-1), adjusted_p(i));
    end
    adjusted_p_original_order = zeros(size(p_vector));
    adjusted_p_original_order(sort_index) = adjusted_p;
    pvals_holm = nan(k, k);
    pvals_holm(idx) = adjusted_p_original_order;
    for i = 1:k-1
        for j = i+1:k
            pvals_holm(j, i) = pvals_holm(i, j);
        end
    end

    disp('Pairwise Comparisons (p-values after Holm Correction) for small difference wave:');
    disp(pvals_holm);

    %% Friedman Test and Pairwise Comparisons (Holm Correction) for Columns 5-8
    dataa = data_2D(:,5:8);
    [p_friedman, tbl, stats] = friedman(dataa, 1, 'off');
    disp('Friedman Test Result (large difference wave):');
    disp(tbl);

    k = size(dataa, 2);
    pvals = nan(k, k);
    for i = 1:k-1
        for j = i+1:k
            pvals(i, j) = signrank(dataa(:, i), dataa(:, j));
            pvals(j, i) = pvals(i, j);
        end
    end

    idx = find(triu(ones(k), 1));
    p_vector = pvals(idx);
    m = length(p_vector);
    [sorted_p, sort_index] = sort(p_vector);
    adjusted_p = zeros(size(sorted_p));
    for i = 1:m
        adjusted_p(i) = min(sorted_p(i) * (m - i + 1), 1);
    end
    for i = 2:m
        adjusted_p(i) = max(adjusted_p(i-1), adjusted_p(i));
    end
    adjusted_p_original_order = zeros(size(p_vector));
    adjusted_p_original_order(sort_index) = adjusted_p;
    pvals_holm = nan(k, k);
    pvals_holm(idx) = adjusted_p_original_order;
    for i = 1:k-1
        for j = i+1:k
            pvals_holm(j, i) = pvals_holm(i, j);
        end
    end

    disp('Pairwise Comparisons (p-values after Holm Correction) for large difference wave:');
    disp(pvals_holm);
%}
    
    % The following section converts data_2D to a table and performs repeated measures ANOVA.
    % Construct variable names
    varNames = cell(1, 8);
    idxVar = 1;
    for stim = 1:2
        for day = 1:4
            varNames{idxVar} = sprintf('S%d_Day%d', stim, day);
            idxVar = idxVar + 1;
        end
    end

    % Convert data_2D to a table and set variable names
    T = array2table(data_2D, 'VariableNames', varNames);

    % Construct within-design table
    Stimulus = [ones(4,1); 2*ones(4,1)];            
    Day = [repelem((1:4)', 1); repelem((1:4)', 1)];      
    WithinDesign = table(Stimulus, Day);

    % Fit repeated measures model using fitrm (48 measurements as dependent variables without covariates)
   %rm = fitrm(T, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WithinDesign);
T_sub = T(:, 5:8);
WithinDesign_sub = WithinDesign(5:8, :);
rm = fitrm(T_sub, sprintf('%s-%s ~ 1', varNames{5}, varNames{8}), 'WithinDesign', WithinDesign_sub);

    % Repeated measures ANOVA (including factors: Stimulus, Day, and their interaction)
    ranovatbl = ranova(rm, 'WithinModel', 'Stimulus*Day');

    % Display the ANOVA results
    % disp(ranovatbl);

    % Pairwise comparisons for the Day factor using the 'multcompare' function
    %compDay = multcompare(rm, 'Day', 'By', 'Stimulus');
   % compDay = multcompare(rm, 'Day','bonferroni');
    compDay = multcompare(rm,'Day','ComparisonType','bonferroni');


%     raw_p = compDay.pValue;
%     adj_p = mafdr(raw_p, 'BHFDR', true);
%     compDay.FDR_p = adj_p;
    % Display only the rows with p-values less than 0.05
    sigRows = compDay.pValue < 0.05;
    disp(compDay(sigRows, :));
    %}
end
end

diary off;
%% each participant waveform
for idev = 1:2
    figure;  % 为每个 idev 创建一个新图
    for md = 1:4
        subplot(2,2,md);
        hold on;
        for isub = 1:14
            plot(time,STD{md,isub,idev}.avg, 'LineWidth', 1);
        end
        title(['day' num2str(md)]);
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
             axis([-0.100 0.6 -2.4 5]);
        hold off;
                   set(gca, 'XTick', -0.1:0.05:0.6);
       grid on;
    end
    % 新版本 MATLAB 可用 sgtitle 添加整体标题
    sgtitle(['Difference wave ' num2str(idev)]);

end

for idev = 1:2
    figure;  % 为每个 idev 创建一个新图
    for md = 1:4
        subplot(2,2,md);
        hold on;
        for isub = 1:14
            plot(time,STD{md,isub,idev}.avg, 'LineWidth', 1);
        end
        title(['day' num2str(md)]);
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
             axis([-0.100 0.6 -2.4 5]);
        hold off;         
        set(gca, 'XTick', -0.1:0.05:0.6);
       grid on;
    end
    % 新版本 MATLAB 可用 sgtitle 添加整体标题
    sgtitle(['STD ' num2str(idev)]);
 
end

for idev = 1:2
    figure;  % 为每个 idev 创建一个新图
    for md = 1:4
        subplot(2,2,md);
        hold on;
        for isub = 1:14
            plot(time,DEV{md,isub,idev}.avg, 'LineWidth', 1);
        end
        title(['day' num2str(md)]);
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
             axis([-0.100 0.6 -2.4 5]);
        hold off;
                  set(gca, 'XTick', -0.1:0.05:0.6);
       grid on; 
    end
    % 新版本 MATLAB 可用 sgtitle 添加整体标题
    sgtitle(['DEV ' num2str(idev)]);

end



 
 %%
 figure
 
 for idev=1
 
 
       % subplot(2,2,i(idev))
        
            plot(Diff_gavg{1, 1}.time, Diff_ga(:,:), 'LineWidth',1); hold on;
        
                
        %%
        graphname = {'day1'; 'day2'; 'day3'; 'day4'};
        title(  deviant_type{idev}, 'fontweight', 'bold');
        
        
        if iday==3
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
            leg = legend('F3-F4-Fz', 'Location', 'Northeast', 'Orientation', 'vertical', 'FontSize', 12, 'fontweight', 'bold', 'box', 'off');
            leg.ItemTokenSize = [9,9];

        end
                    % 在坐标 (0, -2) 处添加文本 "190-260ms"

        
        set(gca, 'FontSize', 10, 'fontweight', 'bold');
        axis([-0.100 0.6 -2.4 2.3]);
        
        % 设置 x 轴刻度，每隔 50ms (0.05秒)一条
       set(gca, 'XTick', -0.1:0.05:0.6);
       grid on;
    end



%% waveform
chan_nam={'F3','Fz','F4','C3','Cz','C4'}; chan_idx=[8,16,9,1,7,2];
deviant_type={'Small deviant','Large deviant'}; i=[1,2,4,5,];
for idev=1:2
    figure
    for iday=1:4
        subplot(2,3,i(iday))
        
            plot(Diff_gavg{1, 1}.time, Diff_gavg{iday, idev}.avg(1,:), 'LineWidth',1); hold on;
        
                
        %%
        graphname = {'day1'; 'day2'; 'day3'; 'day4'};
        title( graphname{iday}, 'fontweight', 'bold');
        
        
        if iday==3
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
            leg = legend('F3-F4-Fz', 'Location', 'Northeast', 'Orientation', 'vertical', 'FontSize', 12, 'fontweight', 'bold', 'box', 'off');
            leg.ItemTokenSize = [9,9];

        end
                    % 在坐标 (0, -2) 处添加文本 "190-260ms"

        
        set(gca, 'FontSize', 10, 'fontweight', 'bold');
        axis([-0.100 0.6 -2.4 2.3]);
        
        % 设置 x 轴刻度，每隔 50ms (0.05秒)一条
       set(gca, 'XTick', -0.1:0.05:0.6);
       grid on;
    end
hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.06, 0.5, deviant_type{idev}, 'Units', 'normalized', 'Rotation', 90, ...
     'FontSize', 20,'fontweight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
