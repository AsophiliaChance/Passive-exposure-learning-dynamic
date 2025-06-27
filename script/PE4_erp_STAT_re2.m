clc;clear
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');

filename={'day1_nsubavg_re.mat','day2_nsubavg_re.mat','day3_nsubavg_re.mat','day4_nsubavg_re.mat'};
%%
for md=1:4
    load(filename{md});
    %reshape ????��???
    for isub=1:14
        for idev=1:2
            Diff{md,isub,idev}=Diff_avg{isub,idev};
            STD{md,isub,idev}=STD_avg{isub,idev};
            DEV{md,isub,idev}=DEV_avg{isub,idev};
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
channelsToAvg = {'F3', 'F4', 'Fz'};
for idev = 1:2
for md = 1:4
    
        data = Diff_gavg{md, idev};
        chanIdx = find(ismember(data.label, channelsToAvg));
        avgData = mean(data.avg(chanIdx, :), 1);
        data.avg = avgData;
        data.label = {'F3_F4_Fz'};
%         Diff_gavg{md, idev} = data;
                idx_min = find(data.time >= 0.19 & data.time <= 0.26);
        [minVal, minIdx] = min(data.avg(idx_min));
        minTime(md,idev) = data.time(idx_min(minIdx));
        idx_max = find(data.time >= 0.29 & data.time <= 0.34);
        [maxVal, maxIdx] = max(data.avg(idx_max));
        maxTime(md,idev) = data.time(idx_max(maxIdx));

    end
end

win=[];
win.mmn1=minTime-0.035; win.mmn2=minTime+0.035;
win.p31=maxTime-0.035; win.p32=maxTime+0.035;
win.mmn=minTime;win.p3=maxTime;
%% Peak detection
data=[];
latency=[]; amplitude=[];
for md = 1:4
    for idev = 1:2
        
        for isub = 1:size(Diff,2)
            data = Diff{md, isub, idev};
            idx = find(data.time >= win.mmn1(md,idev) & data.time <= win.mmn2(md,idev));
            amplitude.mmn(md, isub, idev) = mean(data.avg(idx));
            %[amplitude.mmn(md, isub, idev), peakIdx] = min(data.avg(idx));
            [~, peakIdx] = min(data.avg(idx));
            latency.mmn(md, isub, idev) = data.time(idx(peakIdx));
            
            idx = find(data.time >= win.p31(md,idev) & data.time <= win.p32(md,idev));
            amplitude.p3(md, isub, idev) = mean(data.avg(idx));
            %[amplitude.p3(md, isub, idev), peakIdx] = max(data.avg(idx));
             [~, peakIdx] = max(data.avg(idx));
            latency.p3(md, isub, idev) = data.time(idx(peakIdx));
        end
    end
end
%% ����ͼ
n = 14;
data=[];
data = cat(3, amplitude_diff.mmn, amplitude_diff.p3);
means = mean(data,2);
sem = std(data,0,2)./sqrt(n);

figure(1); clf
groupLabels = {'day1','day2','day3','day4'};
componentTitles = {'MMNs','MMNl','P3s','P3l'};

for c = 1:4
    subplot(2,2,c)
    b = bar(means(:,c));
    hold on
    errorbar(1:4,means(:,c),sem(:,c),'k','LineStyle','none','LineWidth',1.4,'CapSize',10);
    b.FaceColor = [0.2 0.6 0.8];
    b.EdgeColor = 'none';
    %ylim padded
    set(gca,'XTick',1:4,'XTickLabel',groupLabels,'Box','off')
    title(componentTitles{c})
    ylabel('Value')
    hold off
end

sgtitle('Amplitude (Mean �� SEM)')

%% anova
%% MATLAB �ű�����4��14��2��6���ݽ����ظ������������
data = latency.p3;


% ��̬�Լ���
% ����һ����ÿ�������������һ�У�����̬�Խ��м���
[h_norm1, p_norm1] = lillietest(data(:,1));
fprintf('��̬�Լ��� p ֵ: %f\n', p_norm1);

% ������������������������Ƿ������̬�ֲ�����������ƽ��������
[h_norm_all, p_norm_all] = lillietest(data(:));
fprintf('��������任����̬�Լ��� p ֵ: %f\n', p_norm_all);

% �������Լ���
% ���� data_log ��ÿһ�д���һ���������
% ���������
[nRows, nCols] = size(data);
group = repmat(1:nCols, nRows, 1);  % ÿ��Ϊһ��
group = group(:);
data_vector = data(:);

% ʹ�� vartestn ���� Levene ���飨'LeveneAbsolute'��
[p_var, stats] = vartestn(data_vector, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
fprintf('�������Լ��飨Levene���飩 p ֵ: %f\n', p_var);
% ��������
data_2D = reshape(data, 14, []);  % ���Ϊ 14��48 ����
% ���μ���
n = size(data_2D, 1);
k = size(data_2D, 2);
S = cov(data_2D);
W = det(S) / ((trace(S)/k)^k);
c = 1 - (2*k^2 + k + 2) / (6*(k-1)*(n-1));
chi2_stat = -(n-1)*c*log(W);
df = (k*(k-1)/2) - 1;
p = 1 - chi2cdf(chi2_stat, df);
fprintf('Mauchly''s W = %.4f\n', W);
fprintf('chi-square = %.4f, df = %d, p = %.4f\n', chi2_stat, df, p);

% �����������
% ��������˳�����Ϊ Stimulus���ڲ�����Ϊ Day��Channel�������������
varNames = cell(1, 8);
idx = 1;
for stim = 1:2
    for day = 1:4
        
            varNames{idx} = sprintf('S%d_D%d', stim, day);
            idx = idx + 1;
        
    end
end

% �� data_2D ת��Ϊ table�����趨������
T = array2table(data_2D, 'VariableNames', varNames);

% ���� within-design ��

Stimulus = [ones(4,1); 2*ones(4,1)];            
Day = [repelem((1:4)', 1); repelem((1:4)', 1)];      


WithinDesign = table(Stimulus, Day);

% ����ظ�����ģ��
% ʹ�� fitrm �� 48 ������ֵ��Ϊ�����������Э����������ָ�� within-design ��Ϣ
rm = fitrm(T, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WithinDesign);

% �ظ������������
% ָ��ģ���п��� Stimulus��Day��Channel ���佻������
ranovatbl = ranova(rm, 'WithinModel', 'Stimulus*Day');

% ��ʾ����������
disp(ranovatbl);
% �� Day ���ؽ��������Ƚ�
% ��ʽһ��ֱ�Ӵ�����������

compDay = multcompare(rm, 'Day', 'By', 'Stimulus');
raw_p = compDay.pValue;
adj_p = mafdr(raw_p, 'BHFDR', true);
compDay.FDR_p = adj_p;
disp(compDay);
disp(compDay);

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
                    % ������ (0, -2) ������ı� "190-260ms"

        
        set(gca, 'FontSize', 10, 'fontweight', 'bold');
        axis([-0.100 0.6 -2.4 2.3]);
        
        % ���� x ��̶ȣ�ÿ�� 50ms (0.05��)һ��
       set(gca, 'XTick', -0.1:0.05:0.6);
       grid on;
    end
hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.06, 0.5, deviant_type{idev}, 'Units', 'normalized', 'Rotation', 90, ...
     'FontSize', 20,'fontweight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
