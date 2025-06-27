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





