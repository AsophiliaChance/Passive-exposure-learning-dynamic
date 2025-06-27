clear;clc;

filename={'day1_nsubDEVnSTDmvpa_re2.mat','day2_nsubDEVnSTDmvpa_re2.mat',...
    'day3_nsubDEVnSTDmvpa_re2.mat','day4_nsubDEVnSTDmvpa_re2.mat'};
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
            Diff{md,isub,idev}=squeeze(mean(clssi_Diff{isub,idev}.trial(:,chanIdx, :),2));
        end
    end
end

%%
plotData=[]; slide_data=[];slide_all_std=[];slide_data1=[];
time=DEV_keeptrials_nsub{1,1}.time;ax=[];



    for idev=1:2
        i=1;figure('Position', [20, 20, 550, 350]);
        for P3a=[2,1]
        if P3a==1
            
            plotwin= find(time >= 0.281 & time <= 0.351 );%***********************************************************
        else
            plotwin=find(time >= 0.189 & time <= 0.259 );
        end
        
        
        for md=1:4
            ax{i}=subplot(2,4,i);
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
%                 if idev==2
%                 pos(1)=pos(1)+0.01;
%                 set(ax{i}, 'Position', pos);
%             
%         end

i=i+1;
            
            for isub=1:14
                % plot
                
                plotData{md,isub,idev}=squeeze(mean(Diff{md,isub,idev}(:,plotwin),2));hold on
                % slide_data=[];
                for x=1:400%length(plotData{md,isub,idev})-99
                    slide_data(isub,x)=mean(plotData{md,isub,idev}(x:x+99));
                    %slide_data1(x,:)=mean(plotData1{md,isub,idev}((x:x+99),:),1);
                end
                % ntrial(md,isub)=x;
            end
            
            slide_data_gavg=mean(slide_data,1);
           
            slide_std=std(slide_data,0,1);
            mdl = fitlm(1:400, slide_data_gavg);
            plot(mdl.Fitted,'color','k','linewidth',1);hold on
            R2 = mdl.Rsquared.Ordinary;
            xLimits = [1,400];
            if P3a==2
                yLimits=([ -3 4]);
            else
                yLimits=([-1 5]);
            end
text(xLimits(1) + 0.05*(xLimits(2)-xLimits(1)), yLimits(2) - 0.05*(yLimits(2)-yLimits(1)), ...
    sprintf('R^2 = %.2f', R2), 'FontSize', 8, 'FontWeight', 'bold');


if idev==2 && P3a==1 && md==3
    text(xLimits(1) + 0.05*(xLimits(2)-xLimits(1)), yLimits(2) - 0.2*(yLimits(2)-yLimits(1)), ...
    sprintf('slope = %.4f ', mdl.Coefficients.Estimate(2)), 'FontSize', 8, 'FontWeight', 'bold');
else
    text(xLimits(1) + 0.05*(xLimits(2)-xLimits(1)), yLimits(2) - 0.2*(yLimits(2)-yLimits(1)), ...
    sprintf('slope = %.3f ', mdl.Coefficients.Estimate(2)), 'FontSize', 8, 'FontWeight', 'bold');
end
%text(xLimits(1) + 0.05*(xLimits(2)-xLimits(1)), yLimits(2) - 0.35*(yLimits(2)-yLimits(1)), ...
%    sprintf('*'), 'FontSize', 8, 'FontWeight', 'bold');

            %
            X=[];Y=[];upper=[];lower=[];
            % 计算上下边界
            y=slide_data_gavg;
            err= slide_std/sqrt(14);
            x=1:400;
            upper = y + err;
            lower = y - err;
            
            
            X = [x, fliplr(x)];
            Y = [upper, fliplr(lower)];
 fill(X, Y, [166,189,219]/255,'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold on;
 plot(x, slide_data_gavg, 'color','b','linewidth',1);hold on;
            axis square;
            
            title(['Day',num2str(md)],'fontsize',8,'fontweight','bold');
            if i==2||i==6
                ylabel('Amplitude (\muV)');
                xlabel('Sliding trials');
                
            end

                
           
            if P3a==2
                axis([1 400 -3 4]);
            else
                axis([1 400 -1 5]);
            end
            set(gca,'fontsize',8,'LineWidth',1.5,'fontweight','bold');
            box off
            set(gca,'XTick',[1,100,200,300,400]);
            %pos = get(gca, 'Position');  % 获取第五个子图的位置 [x, y, width, height]
            % 在子图左侧创建一个新的、透明的坐标轴用于放置旋转文本
            if i==2
            hTextAxes = axes('Position', [pos(1)-0.095, pos(2)-0.01, 0.04, pos(4)]);
            text(0.5, 0.5, 'MMN', 'FontWeight', 'bold', 'FontSize', 9, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
            set(hTextAxes, 'Visible', 'off', 'HitTest', 'off');  % 隐藏坐标轴线与刻度
            end
             if i==6
            hTextAxes = axes('Position', [pos(1)-0.09, pos(2)-0.01, 0.04, pos(4)]);
            text(0.5, 0.5, 'P3a', 'FontWeight', 'bold', 'FontSize', 9, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
            set(hTextAxes, 'Visible', 'off', 'HitTest', 'off');  % 隐藏坐标轴线与刻度
             end

        end
    
    end

% 在整个 figure 上创建一个隐藏的坐标轴
axes('Position',[0 0 1 1],'Visible','off');
% 在归一化坐标中添加旋转90°的文本
% text(0.08, 0.5, 'Standard', 'Units', 'normalized', ...
%     'Rotation', 90, 'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');

pos1 = get(ax{1}, 'Position');  % [left, bottom, width, height]
pos4 = get(ax{4}, 'Position');
centerX = (pos1(1) + (pos4(1) + pos4(3)))/2;
topY = pos1(2) + pos1(4);
            
change={'Small change','Large change'};
annotation('textbox', [centerX-0.2, topY+0.01, 0.4, 0.1], ...
    'String', change{idev}, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

% pos1 = get(ax{5}, 'Position');  % [left, bottom, width, height]
% pos4 = get(ax{8}, 'Position');
% centerX = (pos1(1) + (pos4(1) + pos4(3)))/2;
% topY = pos1(2) + pos1(4);
% annotation('textbox', [centerX-0.1, topY+0.02, 0.2, 0.1], ...
%     'String', 'Large change', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'FontSize', 12, ...
%     'FontWeight', 'bold');
   
filename = sprintf('sliding_diff%d.tiff', idev);
print(filename, '-dtiff', '-r400');
end
