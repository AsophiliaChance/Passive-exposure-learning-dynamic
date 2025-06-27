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
%% 
%% each participant waveform
for isub = 1:14
    for idev = 2
        figure;  
        for md = 1:4
            subplot(2,2,md);
            hold on;
            time=Diff_avg{md, idev}.time;
            plot(time,Diff{md,isub,idev}.avg, 'LineWidth', 1);        
            set(gca, 'XTick', -0.1:0.05:0.6);
        grid on;   sgtitle(['Sub',num2str(isub), 'Difference wave ' num2str(idev)]);axis([-0.100 0.6 -3 3]);
        end
        title(['day' num2str(md)]);
        ylabel('Amplitude (\muV)');
        xlabel('Time (ms)');
        
        hold off;

       % print(gcf,['Sub',num2str(isub), 'Difference wave ' num2str(idev) '.tif'], '-dtiffn', '-r600');
        filename = fullfile(pwd, sprintf('Sub%dDifference wave %d.tif', isub, idev));
if exist(filename, 'file')
    delete(filename);
end
print(gcf, filename, '-dtiffn', '-r600');
    end
      
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for isub = 1:14
    for idev = 1:2
        figure;  % 为每个 idev 创建一个新图
        for md = 1:4
            subplot(2,2,md);
            hold on;
            
            plot(time,STD{md,isub,idev}.avg, 'LineWidth', 1);
        end
        title(['day' num2str(md)]);
        ylabel('Amplitude (\muV)');
        xlabel('Time (ms)');
        axis([-0.100 0.6 -2.4 5]);
        hold off;
        set(gca, 'XTick', -0.1:0.05:0.6);
        grid on;sgtitle(['Sub',num2str(isub),' STD ' num2str(idev)]);
        print(gcf, ['Sub',num2str(isub), 'STD ' num2str(idev) '.tif'], '-dtiffn', '-r600');
    end
 
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for isub = 1:14
    for idev = 1:2
        figure;  % 为每个 idev 创建一个新图
        for md = 1:4
            subplot(2,2,md);
            hold on;
            
            plot(time,DEV{md,isub,idev}.avg, 'LineWidth', 1);
        end
        title(['day' num2str(md)]);
        ylabel('Amplitude (\muV)');
        xlabel('Time (ms)');
        axis([-0.100 0.6 -2.4 5]);
        hold off;
        set(gca, 'XTick', -0.1:0.05:0.6);
        grid on;
        sgtitle(['Sub',num2str(isub),'DEV ' num2str(idev)]);
        print(gcf, ['Sub',num2str(isub), 'DEV' num2str(idev) '.tif'], '-dtiffn', '-r600');
    end
   
end



