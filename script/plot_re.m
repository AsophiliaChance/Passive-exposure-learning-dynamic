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
%%
%% topo
chan_nam={'F3','Fz','F4','C3','Cz','C4'}; chan_idx=[8,16,9,1,7,2];
% ���û�ͼ����
figure
for idev = 1:2
    
    for md = 1:4
        subplot(4,4,i)
        
%         t = 100; % �޸�Ϊ����Ҫ��ʱ�������
% dataVector = EEG.data(:, t); % ��ȡ��ʱ������е缫����ֵ
% topoplot(dataVector, EEG.chanlocs, 'maplimits', 'absmax'); 
% colorbar;  % ��ѡ��Ϊͼ�������ɫ��

        
        cfg = [];
        cfg.layout    = 'elec1020.lay';   % �������ʵ��ѡ����ʵ� layout �ļ�
        cfg.xlim      = [0.281 0.351];% [0.189 0.259];%[0.27 0.38];      % ������Ҫ��ʾ��ʱ�䴰[0.19 0.26]
        cfg.zlim      = [-4 4];             % �������ݵ���ɫ��Χ
        cfg.comment   = 'no';               % ȡ��ͼ��Ĭ�ϵ�ע��
        cfg.marker    = 'labels';           % ��ʾ�缫���Ʊ�ǩ
         cfg.figure   = 'no';
        
        % ���� ERP ����ͼ
        ft_topoplotER(cfg,Diff_gavg{md,idev});
        

    end
            % ��ӱ��⣬��ע��ǰ�� day ���
        title(['Day ' num2str(md)], 'FontWeight', 'bold', 'FontSize', 12);
                filename = sprintf('ERP_P3_Topo_Day%d_Diff%d.tif', md, idev);
        print(gcf, filename, '-dtiff', '-r300');
        colorbar;
end
%% topo
chan_nam={'F3','Fz','F4','C3','Cz','C4'}; chan_idx=[8,16,9,1,7,2];
% ���û�ͼ����
for idev = 1:2
    for md = 1:4
        
        cfg = [];
        cfg.layout    = 'elec1020.lay';   % �������ʵ��ѡ����ʵ� layout �ļ�
        cfg.xlim      = [0.281 0.351];% [0.189 0.259];%[0.27 0.38];      % ������Ҫ��ʾ��ʱ�䴰[0.19 0.26]
        cfg.zlim      = [-4 4];             % �������ݵ���ɫ��Χ
        cfg.comment   = 'no';               % ȡ��ͼ��Ĭ�ϵ�ע��
        cfg.marker    = 'labels';           % ��ʾ�缫���Ʊ�ǩ
        
        % ���� ERP ����ͼ
        ft_topoplotER(cfg,Diff_gavg{md,idev});
        
        % ��ӱ��⣬��ע��ǰ�� day ���
        title(['Day ' num2str(md)], 'FontWeight', 'bold', 'FontSize', 12);
                filename = sprintf('ERP_P3_Topo_Day%d_Diff%d.tif', md, idev);
        print(gcf, filename, '-dtiff', '-r300');
        colorbar;
    end
end


%%
channelsToAvg = {'F3', 'F4', 'Fz'};
for md = 1:4
    for idev = 1:2
        data = Diff_gavg{md, idev};
        chanIdx = find(ismember(data.label, channelsToAvg));
        avgData = mean(data.avg(chanIdx, :), 1);
        data.avg = avgData;
        data.label = {'F3_F4_Fz'};
        Diff_gavg{md, idev} = data;
                idx_min = find(data.time >= 0.19 & data.time <= 0.26);
        [minVal, minIdx] = min(data.avg(idx_min));
        minTime(md,idev) = data.time(idx_min(minIdx));
        idx_max = find(data.time >= 0.29 & data.time <= 0.34);
        [maxVal, maxIdx] = max(data.avg(idx_max));
        maxTime(md,idev) = data.time(idx_max(maxIdx));

    end
end
%%
% waveform
deviant_type={'Small deviant','Large deviant'}; i=[1,2,4,5,];
for idev=1:2
    figure
    for iday=1:4
        subplot(2,3,i(iday))
        for  ichan=1:6
            plot(Diff_gavg{1, 1}.time, Diff_gavg{iday, idev}.avg, 'LineWidth',1); hold on;
        end
        x1 = 0.180;
               if idev==2
        y1 = -2.25;
        else
            y1=-1.3;
        end
        width1 = 0.080;
        height1 = 1.5;
        rectangle('Position', [x1, y1, width1, height1], 'EdgeColor', [0.3 0.3 0.3], 'LineStyle', '--','LineWidth', 1);
        
        % �ڶ�������
        % x ��Χ��270 �� 380����� = 380 - 270 = 110
        % y ��Χ��0.5 �� 2.5���ײ�Ϊ 0.5���߶� = 2.5 - 0.5 = 2.0
        x2 = 0.270;
        y2 = 0.25;
        width2 = 0.110;
        height2 = 2.0;
        rectangle('Position', [x2, y2, width2, height2], 'EdgeColor', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);
                %%
                x1 = 0.190; 
        if idev==2
        y1 = -2.25;
        else
            y1=-1.3;
        end
        width1 = 0.050;
        height1 = 1.5;
       
        rectangle('Position', [x1, y1, width1, height1], 'EdgeColor', 'r', 'LineStyle', '--','LineWidth', 1);
        
        
        % �ڶ�������
        % x ��Χ��270 �� 380����� = 380 - 270 = 110
        % y ��Χ��0.5 �� 2.5���ײ�Ϊ 0.5���߶� = 2.5 - 0.5 = 2.0
        x2 = 0.250;
        y2 = 0.25;
        width2 = 0.050;
        height2 = 2.0;
        rectangle('Position', [x2, y2, width2, height2], 'EdgeColor','r', 'LineStyle', '--', 'LineWidth', 1);
        %%
        graphname = {'day1'; 'day2'; 'day3'; 'day4'};
        title( graphname{iday}, 'fontweight', 'bold');
        
        
        if iday==3
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
            leg = legend(chan_nam, 'Location', 'Northeast', 'Orientation', 'vertical', 'FontSize', 12, 'fontweight', 'bold', 'box', 'off');
            leg.ItemTokenSize = [9,9];

        end
                    % ������ (0, -2) ������ı� "190-260ms"

if idev==1
text(0.265, -1.3, '180-260ms', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
text(0.265, -1.1, '190-240ms', 'FontSize', 10, 'Color','r');%
else
text(0.265, -2, '180-260ms', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
text(0.265, -1.8, '190-240ms', 'FontSize', 10, 'Color','r');
end
% ������ (0.1, 2) ������ı� "270-380ms"
text(0.1, 2, '270-380ms', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
text(0.1, 1.8, '250-300ms', 'FontSize', 10, 'Color','r');
        
        set(gca, 'FontSize', 17, 'fontweight', 'bold');
        axis([-0.100 0.6 -2.4 2.3]);
        
        % ���� x ��̶ȣ�ÿ�� 50ms (0.05��)һ��
       % set(gca, 'XTick', -0.1:0.05:0.6);
       % grid on;
    end
hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.06, 0.5, deviant_type{idev}, 'Units', 'normalized', 'Rotation', 90, ...
     'FontSize', 20,'fontweight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
%%
deviant_type={'Small deviant','Large deviant'};
for idev=1:2
    figure
    for ichan=1:6
        subplot(2,3,ichan)
        for  iday=1:4
            plot(Diff_gavg{1, 1}.time, Diff_gavg{iday, idev}.avg, 'LineWidth',1); hold on;
        end
        x1 = 0.180; 
        if idev==2
        y1 = -2.25;
        else
            y1=-1.3;
        end
        width1 = 0.080;
        height1 = 1.5;
       
        rectangle('Position', [x1, y1, width1, height1], 'EdgeColor', [0.3 0.3 0.3], 'LineStyle', '--','LineWidth', 1);
        
        
        % �ڶ�������
        % x ��Χ��270 �� 380����� = 380 - 270 = 110
        % y ��Χ��0.5 �� 2.5���ײ�Ϊ 0.5���߶� = 2.5 - 0.5 = 2.0
        x2 = 0.270;
        y2 = 0.25;
        width2 = 0.110;
        height2 = 2.0;
        rectangle('Position', [x2, y2, width2, height2], 'EdgeColor', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);
        %%
                x1 = 0.190; 
        if idev==2
        y1 = -2.25;
        else
            y1=-1.3;
        end
        width1 = 0.050;
        height1 = 1.5;
       
        rectangle('Position', [x1, y1, width1, height1], 'EdgeColor', 'r', 'LineStyle', '--','LineWidth', 1);
        
        
        % �ڶ�������
        % x ��Χ��270 �� 380����� = 380 - 270 = 110
        % y ��Χ��0.5 �� 2.5���ײ�Ϊ 0.5���߶� = 2.5 - 0.5 = 2.0
        x2 = 0.250;
        y2 = 0.25;
        width2 = 0.050;
        height2 = 2.0;
        rectangle('Position', [x2, y2, width2, height2], 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1);
        %%
        graphname = {'day1'; 'day2'; 'day3'; 'day4'};
        title(chan_nam{ichan}, 'fontweight', 'bold');
        
        
        if ichan==4
            ylabel('Amplitude (\muV)');
            xlabel('Time (ms)');
            leg = legend(graphname, 'Location', 'Northeast', 'Orientation', 'vertical', 'FontSize', 12, 'fontweight', 'bold', 'box', 'off');
            leg.ItemTokenSize = [9,9];

        end
                    % ������ (0, -2) ������ı� "190-260ms"

if idev==1
text(0.265, -1.3, '180-260ms', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
text(0.265, -1.1, '190-240ms', 'FontSize', 10, 'Color','r');%text(0.265, -1.8, '190-240ms', 'FontSize', 10, 'Color','r');
else
text(0.265, -2, '180-260ms', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
text(0.265, -1.8, '190-240ms', 'FontSize', 10, 'Color','r');
end

% ������ (0.1, 2) ������ı� "270-380ms"
text(0.1, 2, '270-380ms', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
text(0.1, 1.8, '250-300ms', 'FontSize', 10, 'Color','r');
        
        set(gca, 'FontSize', 17, 'fontweight', 'bold');
        axis([-0.100 0.6 -2.4 2.3]);
        
        % ���� x ��̶ȣ�ÿ�� 50ms (0.05��)һ��
       % set(gca, 'XTick', -0.1:0.05:0.6);
       % grid on;
    end
hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.06, 0.5, deviant_type{idev}, 'Units', 'normalized', 'Rotation', 90, ...
     'FontSize', 20,'fontweight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
%%
channelsToAvg = {'F3', 'F4', 'Fz'};
for md = 1:4
    for idev = 1:2
        data = Diff_gavg{md, idev};
        chanIdx = find(ismember(data.label, channelsToAvg));
        avgData = mean(data.avg(chanIdx, :), 1);
        data.avg = avgData;
        data.label = {'F3_F4_Fz'};
        Diff_gavg{md, idev} = data;
                idx_min = find(data.time >= 0.19 & data.time <= 0.26);
        [minVal, minIdx] = min(data.avg(idx_min));
        minTime(md,idev) = data.time(idx_min(minIdx));
        idx_max = find(data.time >= 0.29 & data.time <= 0.34);
        [maxVal, maxIdx] = max(data.avg(idx_max));
        maxTime(md,idev) = data.time(idx_max(maxIdx));

    end
end

%%
win.mmn1=minTime-0.03; win.mmn2=minTime+0.03;
win.p31=maxTime-0.03; win.P31=maxTime+0.03;
