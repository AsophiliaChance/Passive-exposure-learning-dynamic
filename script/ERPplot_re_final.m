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
eeglab;
%% grand average for plotting

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual =  'no';
for md=1:4
    for idev=1:2
        a=STD(md,:,idev);
        STD_gavg {md,idev}      = ft_timelockgrandaverage(cfg, a{1}, a{2}, a{3}, a{4}, a{5}, a{6}, a{7}, a{8}, a{9}, a{10}, a{11}, a{12}, a{13}, a{14});
        a=DEV(md,:,idev);
        DEV_gavg {md,idev}      = ft_timelockgrandaverage(cfg, a{1}, a{2}, a{3}, a{4}, a{5}, a{6}, a{7}, a{8}, a{9}, a{10}, a{11}, a{12}, a{13}, a{14});
        a=Diff(md,:,idev);
        Diff_gavg {md,idev}      = ft_timelockgrandaverage(cfg, a{1}, a{2}, a{3}, a{4}, a{5}, a{6}, a{7}, a{8}, a{9}, a{10}, a{11}, a{12}, a{13}, a{14});
    end
    
end

%% topo
chan_nam={'F3','Fz','F4','C3','Cz','C4'}; chan_idx=[8,16,9,1,7,2];
% ���û�ͼ����
change={'Small change','Large change'};
for idev = 1:2
    % ������Ϊ��λ����ͼ�δ��ڵ�λ�úͳߴ�
    fig = figure('Position', [100, 100, 800, 250]);
    i=1;
    for md = 1:4
        subplot(1,4,i)
        cfg = [];
        cfg.layout    = 'elec1020.lay';   % �������ʵ��ѡ����ʵ� layout �ļ�
        cfg.xlim      =  [0.189 0.259];% [0.281 0.351];%    % ������Ҫ��ʾ��ʱ�䴰[0.19 0.26]
        cfg.zlim      = [-4 4];             % �������ݵ���ɫ��Χ
        cfg.comment   = 'no';               % ȡ��ͼ��Ĭ�ϵ�ע��
        cfg.marker    = 'labels';           % ��ʾ�缫���Ʊ�ǩ
        ax = subplot(1,4,i);
        
        pos = get(ax, 'Position');
        % ���磬�����ұ߾���С��������󣻾�����ֵ�ɸ�����Ҫ΢��
        newPos = [pos(1)-0.02*(md-1), pos(2)-0.2, pos(3)+0.01, pos(4)+0.01];
        cfg.figure = subplot('Position', newPos); % ���� subplot
        ft_topoplotER(cfg,Diff_gavg{md,idev});
        hTitle = title(['Day ' num2str(md)], 'FontWeight', 'bold', 'FontSize', 25); % ���ñ���
        
        % ���ñ���λ�ã����Ʊ�����ͼ��֮��ľ���
        set(hTitle, 'Units', 'normalized'); % ʹ���������
        titlePos = get(hTitle, 'Position');
        if md==1
        titlePos(2) =0.84;  % �޸� y ֵ��>1 �����ƣ�<1 ������
        set(hTitle, 'Position', titlePos);
        end
        
        % ��Ϊ�������ͼʱ�������ת�ı� "deviant"
        if i == 1
            pos = get(gca, 'Position');  % ��ȡ�������ͼ��λ�� [x, y, width, height]
            % ����ͼ��ഴ��һ���µġ�͸�������������ڷ�����ת�ı�
            hTextAxes = axes('Position', [pos(1)-0.05, pos(2)-0.01, 0.04, pos(4)]);
            text(0.5, 0.5,'MMN', 'FontWeight', 'bold', 'FontSize', 25, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
            set(hTextAxes, 'Visible', 'off', 'HitTest', 'off');  % ��������������̶�
        end
        
        axis square;
        
        i=i+1;
    end
    % colorbar('LOCATION','EastOutside');
    
        annotation('textbox', [0.3, 0.93, 0.4, 0.05], ... % ����Ը�����Ҫ������Щ��ֵ
            'String',change{idev}, ...
            'EdgeColor', 'none', ...  % ȡ���߿�
            'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 30);
    filename = sprintf('topoMMN%d.tiff', idev);
    print(filename, '-dtiff', '-r350');
end
%colorbar

%%
channelsToAvg = {'F3', 'F4', 'Fz'};
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

for md = 1:4
    for idev = 1:2
        data = Diff_gavg{md, idev};
        channelsToAvg = {'F3', 'F4', 'Fz'};
        chanIdx = find(ismember(data.label, channelsToAvg));
        avgData = mean(data.avg(chanIdx, :), 1);
        data.avg = avgData;
        data.label = {'F3_F4_Fz'};
        Diff_gavg{md, idev} = data;
        
        data = DEV_gavg{md, idev};
        chanIdx = find(ismember(data.label, channelsToAvg));
        avgData = mean(data.avg(chanIdx, :), 1);
        data.avg = avgData;
        data.label = {'F3_F4_Fz'};
        DEV_gavg{md, idev} = data;
        
        data = STD_gavg{md, idev};
        chanIdx = find(ismember(data.label, channelsToAvg));
        avgData = mean(data.avg(chanIdx, :), 1);
        data.avg = avgData;
        data.label = {'F3_F4_Fz'};
        STD_gavg{md, idev} = data;
        
        
    end
end

%%
% waveform
clc
deviant_type={'Small deviant','Large deviant'};% i=[1,2,3,4];
ax=[];
dataplot{3}=Diff_gavg;dataplot{1}=STD_gavg;dataplot{2}=DEV_gavg;
for idev=1:2
    ax=[];
    figure('Position', [100, 100, 800, 300]);
    
    for idata=1:3
        data=dataplot{idata};
        
        ax{idata}= subplot(1,3,idata);
        
        for iday=1:4
            plot(STD_gavg{1, 1}.time, data{iday, idev}.avg, 'LineWidth',2); hold on;
            %             plot(Diff_gavg{1, 1}.time, STD_gavg{iday, idev}.avg, 'LineWidth',2); hold on;
            %             plot(Diff_gavg{1, 1}.time, DEV_gavg{iday, idev}.avg, 'LineWidth',2); hold on;
            %
            
            x1 = 0.189;
            y1=-4;
            width1 = 0.07;
            height1 = 8;
            rectangle('Position', [x1, y1, width1, height1], 'EdgeColor', [0.3 0.3 0.3], 'LineStyle', '--','LineWidth', 1);
            
            % �ڶ�������
            % x ��Χ��270 �� 380����� = 380 - 270 = 110
            % y ��Χ��0.5 �� 2.5���ײ�Ϊ 0.5���߶� = 2.5 - 0.5 = 2.0
            x2 = 0.281;
            y2 = -4;
            width2 = 0.07;
            height2 = 8;
            rectangle('Position', [x2, y2, width2, height2], 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1);
            %%
            graphname = { 'Standard'; 'Deviant';'Differential'};
            title( graphname{idata}, 'fontweight', 'bold');
            
            
            if idata==1
                ylabel('Amplitude (\muV)');
                xlabel('Time (s)');
                leg = legend({'Day1'; 'Day2'; 'Day3'; 'Day4'}, 'Location', 'Northeast', 'Orientation', 'vertical', 'FontSize', 10, 'fontweight', 'bold', 'box', 'off');
                leg.ItemTokenSize = [9,9];
                set(leg, 'Units', 'normalized');  % ʹ�ù�һ����λ���㶨λ
                pos = get(leg, 'Position');  % pos = [left, bottom, width, height]
                
                % ��ͼ�����������ƽ�ư��ͼ�����
               % pos(2)=pos(2)+0.04;
                pos(1) = pos(1) + pos(3)/4;
                set(leg, 'Position', pos);
                
            end
            % ������ (0, -2) ������ı� "190-260ms"
            
            
            text(0.215, -3.9, 'MMN', 'FontSize', 11, 'Color',[0.3 0.3 0.3], 'Rotation', 90);
            text(0.315, -3.9, 'P3a', 'FontSize', 11, 'Color','r', 'Rotation', 90);
            
            
            
            set(gca, 'FontSize', 12, 'fontweight', 'bold');
            axis([-0.100 0.6 -4 4]);
            xticks([0,0.2, 0.4, 0.6]);
            
            % ���� y ��̶�λ��
            yticks([-4, -2, 0, 2, 4]);
            set(gca, 'box', 'off');
            axis square;
            set(gca, 'TickLength', [0.02, 0.025]);
            
            
            
            % ���� x ��̶ȣ�ÿ�� 50ms (0.05��)һ��
            % set(gca, 'XTick', -0.1:0.05:0.6);
            % grid on;
        end
        % ����ÿ����ͼ��λ�úʹ�С�����ټ��
        
        
        % hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
        %text(0.06, 0.5, deviant_type{idev}, 'Units', 'normalized', 'Rotation', 90, ...
        %      'FontSize', 20,'fontweight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    for i = 1:length(ax)
        ax1=ax;
        pos = get(ax{i}, 'Position');
        pos(1) = pos(1) - 0.03*(i-1);  % �����ƶ�һ��
        pos(2) = pos(2) - 0.02;  % �����ƶ�һ��
         pos(3) = pos(3) - 0.025;
        pos(4) = pos(4) - 0.02 ;
        set(ax{i}, 'Position', pos);
    end
    
    change={'Small change','Large change'};
    annotation('textbox', [0.3, 0.94, 0.4, 0.05], ... % ����Ը�����Ҫ������Щ��ֵ
        'String',change{idev}, ...
        'EdgeColor', 'none', ...  % ȡ���߿�
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 20);
    filename = sprintf('ERP%d.tiff', idev);
    print(filename, '-dtiff', '-r600');
    
end
