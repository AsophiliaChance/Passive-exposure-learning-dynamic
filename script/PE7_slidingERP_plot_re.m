clear;clc;

filename={'day1_nsubDEVnSTDmvpa_re.mat','day2_nsubDEVnSTDmvpa_re.mat',...
    'day3_nsubDEVnSTDmvpa_re.mat','day4_nsubDEVnSTDmvpa_re.mat'};
eeglab
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');
channelsToAvg = {'F3', 'F4', 'Fz'};
% chanIdx = find(ismember(DEV_keeptrials_nsub{1,1}.label, channelsToAvg));


for md=1:4
    load(filename{md});
    chanIdx = find(ismember(DEV_keeptrials_nsub{1,1}.label, channelsToAvg));
    %reshape 每个变量
    for isub=1:14
        for idev=1:2
            DEV{md,isub,idev}=squeeze(mean(DEV_keeptrials_nsub{isub,idev}.trial(:,chanIdx, :),2));
            STD{md,isub,idev}=squeeze(mean(STD_keeptrials_nsub{isub,idev}.trial(:,chanIdx, :),2));
                        
            Diff{md,isub,idev}=DEV{md,isub,idev}(1:502,:)- STD{md,isub,idev}(1:502,:);
        end
    end
end

%%
plotData=[]; slide_data=[];slide_all_std=[];slide_data1=[];
time=DEV_keeptrials_nsub{1,1}.time;ax=[];
i=1;figure('Position', [20, 20, 1100, 350]);

for P3a=[2,1]
    for idev=1:2
        if P3a==1
            
            plotwin= find(time >= 0.281 & time <= 0.351 );%***********************************************************
        else
            plotwin=find(time >= 0.189 & time <= 0.259 );
        end
        
        
        for md=1:4
            ax{i}=subplot(2,8,i);
    pos = get(ax{i}, 'Position');
    % 根据需要手动调整 pos 中的 [left, bottom, width, height]
    % 例如减少左右间距
    %pos(1) = pos(1) + 0.02;
    pos(3) = pos(3) + 0.005;
    set(ax{i}, 'Position', pos);
            if P3a==2
                pos(2)=pos(2)-0.065;
                set(ax{i}, 'Position', pos);
            
            end
                if idev==2
                pos(1)=pos(1)+0.01;
                set(ax{i}, 'Position', pos);
            
        end

i=i+1;
            
            for isub=1:14
                % plot
                
                plotData{md,isub,idev}=squeeze(mean(STD{md,isub,idev}(:,plotwin),2));
                % slide_data=[];
                for x=1:400%length(plotData{md,isub,idev})-99
                    slide_data(isub,x)=mean(plotData{md,isub,idev}(x:x+99));
                    %slide_data1(x,:)=mean(plotData1{md,isub,idev}((x:x+99),:),1);
                end
                % ntrial(md,isub)=x;
            end
            
            slide_data_gavg=mean(slide_data,1);
           
            slide_std=std(slide_data,0,1);
            
            
            X=[];Y=[];upper=[];lower=[];
            % 计算上下边界
            y=slide_data_gavg;
            err= slide_std/14;
            x=1:400;
            upper = y + err;
            lower = y - err;
            
            
            X = [x, fliplr(x)];
            Y = [upper, fliplr(lower)];
            fill(X, Y,[166,189,219]/255,'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold on;
            plot(x, slide_data_gavg, 'color','b','linewidth',1);hold on;
            axis square;
            
            title(['Day',num2str(md)],'fontsize',9,'fontweight','bold');
            if i==2||i==10
                ylabel('Amplitude (\muV)');
                xlabel('Datapoints');
                
            end

                
           
            if P3a==2
                axis([0 400 -3 4]);
            else
                axis([0 400 -3 4]);
            end
            set(gca,'fontsize',9,'LineWidth',1.5,'fontweight','bold');
            box off
            set(gca,'XTick',0:100:400);
            %pos = get(gca, 'Position');  % 获取第五个子图的位置 [x, y, width, height]
            % 在子图左侧创建一个新的、透明的坐标轴用于放置旋转文本
            if i==2
            hTextAxes = axes('Position', [pos(1)-0.055, pos(2)-0.01, 0.04, pos(4)]);
            text(0.5, 0.5, 'MMN', 'FontWeight', 'bold', 'FontSize', 9, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
            set(hTextAxes, 'Visible', 'off', 'HitTest', 'off');  % 隐藏坐标轴线与刻度
            end
             if i==10
            hTextAxes = axes('Position', [pos(1)-0.057, pos(2)-0.01, 0.04, pos(4)]);
            text(0.5, 0.5, 'P3a', 'FontWeight', 'bold', 'FontSize', 9, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
            set(hTextAxes, 'Visible', 'off', 'HitTest', 'off');  % 隐藏坐标轴线与刻度
             end

        end
    
    end
end
% 在整个 figure 上创建一个隐藏的坐标轴
axes('Position',[0 0 1 1],'Visible','off');
% 在归一化坐标中添加旋转90°的文本
text(0.08, 0.5, 'Standard', 'Units', 'normalized', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');

pos1 = get(ax{1}, 'Position');  % [left, bottom, width, height]
pos4 = get(ax{4}, 'Position');
centerX = (pos1(1) + (pos4(1) + pos4(3)))/2;
topY = pos1(2) + pos1(4);
annotation('textbox', [centerX-0.1, topY+0.02, 0.2, 0.1], ...
    'String', 'Small change', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

pos1 = get(ax{5}, 'Position');  % [left, bottom, width, height]
pos4 = get(ax{8}, 'Position');
centerX = (pos1(1) + (pos4(1) + pos4(3)))/2;
topY = pos1(2) + pos1(4);
annotation('textbox', [centerX-0.1, topY+0.02, 0.2, 0.1], ...
    'String', 'Large change', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

filename = sprintf('sliding_std.tiff');
print(filename, '-dtiff', '-r300');
