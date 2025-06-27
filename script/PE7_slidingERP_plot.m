%% load files
clc;clear
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult')
is100trial=0;
isremovebaseline=0;
ishighpass=0;
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
filename={[filler1,filler2,filler3,'day1_nsubDEVnSTDmvpa.mat'],[filler1,filler2,filler3,'day2_nsubDEVnSTDmvpa.mat'],...
    [[filler1,filler2,filler3],'day3_nsubDEVnSTDmvpa.mat'],[filler1,filler2,filler3,'day4_nsubDEVnSTDmvpa.mat']};
for md=1:4
    load(filename{md});
    %reshape 每个变量
    for isub=1:14
        for idev=1:2
            Diff{md,isub,idev}=clssi_Diff{isub,idev}.trial;
        end
    end
end

%%

load('G:\passiveexp22222\data\3rdanalysis\MVPA_dataNresult\debsln_0_1_erpstat.mat')

time=ERP_stat{1, 1}.time;
label=ERP_stat{1, 2}.label;

chaninx=[1,3,2,4,8,5,6,9,7];
for i=1:length(chaninx)
    siglabel{i}=ERP_stat{1, 2}.label(chaninx(i));
end

%%
plotData=[]; slide_data=[];slide_all_std=[];slide_data1=[];
for P3a=1:2
    for idev=1:2
    if P3a==1
        
        plotwin=89:114;
    else
        plotwin=63:88;
    end
    figure
    for i=1:4
        a1=subplot(2,2,i);
        pos_main(i,:)=get(a1,'position');
    end
    pos_main(:,3)=0.3;pos_main(:,4)=0.3;
    w=0.35;
    x1=0.13;x2=0.13+w;x3=x2+w;
    ind2=find(pos_main==pos_main(2,1));
    pos_main(ind2)=x2;
    h=0.35;
    y1=0.5;y2=y1-h;
    yind1=find(pos_main==0.5838);
    pos_main(yind1)=y1;
    yind2=find(pos_main==0.1100);
    pos_main(yind2)=y2;
    % ind3=find(pos_main==0.6916);
    % pos_main(ind3)=x3;
    chan_oi=find(ismember(label,'Cz'));
    for md=1:4
        subplot('Position',pos_main(md,:));
        for isub=1:14
            % plot
            
            plotData{md,isub,idev}=squeeze(mean(mean(Diff{md,isub,idev}(:,:,plotwin),3),2));
            plotData1{md,isub,idev}=squeeze(mean(Diff{md,isub,idev}(:,:,plotwin),3));
            for x=1:length(plotData{md,isub,idev})-99
                slide_data(x)=mean(plotData{md,isub,idev}(x:x+99));
                slide_data1(x,:)=mean(plotData1{md,isub,idev}((x:x+99),:),1);
            end
            if md==4 && isub==9
                x=492;
            end
           slide_all{P3a,md,isub,idev}=slide_data1; 
           slide_all_std{P3a,idev}(isub,md,:)=std(slide_data1,0,1);
            plot(1:1:x, slide_data1(1:x,chan_oi), 'linewidth',1.5);hold on;
            
            
            
        end
        title(['Day',num2str(md)],'fontsize',9,'fontweight','bold');
        if md==3
            ylabel('Amplitude (\muV)');
            xlabel('Trials');
        end
        
        axis([0 900 -8 8]);
        set(gca,'fontsize',9,'LineWidth',1.5,'fontweight','bold');
        box off
        set(gca,'XTick',0:200:900);
        set(gca,'YTick',-6:2:8);
    end
    end
end
%%
%{
p=[];h=[];
alpha = 0.05 /144;

for idev=1:2       
    for channel=1:9
            for P3a=1:2
    i=1;
    for md=1:3
        for md2=md+1:4

            
[h{idev,P3a}(channel,i), p{idev,P3a}(channel,i)] = ttest(slide_all_std{P3a,idev}(:,md,channel), slide_all_std{P3a,idev,:}(:,md2,channel),'Tail', 'right','Alpha', alpha);
i=i+1;
        end
        end
    end
    end
end
%}
%%
slope_p=[];
for P3a=1:2
    for md=1:4
        for isub=1:14
            for idev=1:2
                for ichan=1:9
                x = 1:length(slide_all{P3a,md,isub,idev}(:,ichan)); % 自变量
                y = slide_all{P3a,md,isub,idev}(:,ichan); % 随机数据
                
                % 构建线性回归模型
                lm = fitlm(x, y);
                
                % 提取回归系数表
                coefficients = lm.Coefficients;
                
                % 提取斜率的p值
                slope{P3a,idev}(md,isub,ichan) = coefficients.Estimate(2);
                slope_p{P3a,idev}(md,isub,ichan) = coefficients.pValue(2);  % 第二行是斜率的p值
                % 获取95%置信区间
                ci = coefCI(lm, 0.05); % 默认95%置信区间
                % ci是一个矩阵，每行对应一个系数 [下界，上界]
                slope_ci = ci(2,:); % 第2行对应斜率系数
                if slope_p{P3a,idev}(md,isub,ichan) < 0.05  && (slope_ci(1)*slope_ci(2) > 0)
                    slope_sig{P3a,idev}(md,isub,ichan)=slope{P3a,idev}(md,isub,ichan);
                else
                    slope_sig{P3a,idev}(md,isub,ichan)=0;
                end
                end
               
            end
        end
    end
end
slope_ori=slope_sig;
tmp=slope_ori{2,1};
%%
irow=1; perct=[];slope_CI=[];
for P3a=1:2        
    for idev=1:2
    if P3a==2
        slope_sig{P3a,idev}(:,:,:)=-(slope_ori{P3a,idev}(:,:,:));       
    end
    for md=1:4        
icol=1;
            for ichan=1:9

                data = slope_sig{P3a, idev}(md, :, ichan);

                % 将数据转换为一维向量
                data = squeeze(data);
                
                % 统计大于零的个数
                num_positive = sum(data > 0);
                
                % 统计小于零的个数
                num_negative = sum(data < 0);
                
                % 总的数据点数
                total_count = numel(data);
                
                % 计算大于零的百分比
                percent_positive{P3a,idev}(md,ichan) = (num_positive / total_count) * 100;
                perct(irow,icol)=num_positive;
                perct_plt(irow,icol)=num_positive;
                icol=icol+1;
               
                % 计算小于零的百分比
                percent_negative{P3a,idev}(md,ichan) = (num_negative / total_count) * 100; 
                perct(irow,icol)=num_negative;
                 perct_plt(irow,icol)=-num_negative;
                icol=icol+1;
                slope_CI(P3a,idev,md,ichan,1) = min(slope{P3a,idev}(md,:,ichan));
                slope_CI(P3a,idev,md,ichan,2) = max(slope{P3a,idev}(md,:,ichan));
            end
            irow=irow+1;
        end       
    end
end

min(slope_CI,[],'all')
max(slope_CI,[],'all')

%%
% 假设 perct_plt 为您的数据矩阵
figure;
imagesc(perct_plt);

% 应用自定义 colormap
color1 = [103,169,207]/255; % 蓝
color2 = [1,1,1];           % 白
color3 = [239,138,98]/255;  % 红
numColors = 200;
cmap = [linspace(color1(1), color2(1), numColors/2)', linspace(color1(2), color2(2), numColors/2)', linspace(color1(3), color2(3), numColors/2)';
        linspace(color2(1), color3(1), numColors/2)', linspace(color2(2), color3(2), numColors/2)', linspace(color2(3), color3(3), numColors/2)'];
colormap(cmap);

% 设置颜色范围，0在中间
caxis([-14 14]);
colorbar;

% 在图上叠加文字显示绝对值
[numRows, numCols] = size(perct_plt);
for i = 1:numRows
    for j = 1:numCols
        val = perct_plt(i,j);
        text(j, i, num2str(abs(val),'%.f'), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k','FontSize',17,'FontName','Times New Roman','FontWeight','bold');
    end
end



%%
%{
figure
slo=[];

for P3a=[2,1]
    for idev=1:2
        if P3a==2
            slope_sig{P3a,idev}(:,:)=-(slope_sig{P3a,idev}(:,:));
        end
            
slo=[slo;slope_sig{P3a,idev}(:,:)];
end
end
h=heatmap(slo);
% 定义自定义颜色映射
color1 = [107,174,214] / 255; % 负值端颜色：#bdd7e7 (蓝色)
color2 = [1, 1, 1];              % 中间颜色：白色
color3 = [253,174,97] / 255;  % 正值端颜色：#feebe2 (红色)

% 生成颜色映射矩阵，线性插值从负值颜色到中间颜色再到正值颜色
numColors = 200;
cmap = [linspace(color1(1), color2(1), numColors/2)', linspace(color1(2), color2(2), numColors/2)', linspace(color1(3), color2(3), numColors/2)';
        linspace(color2(1), color3(1), numColors/2)', linspace(color2(2), color3(2), numColors/2)', linspace(color2(3), color3(3), numColors/2)'];

% 应用自定义颜色映射
colormap(h, cmap);

% 设置颜色范围，确保0在中间
h.ColorLimits = [-max(abs(slo(:))), max(abs(slo(:)))];

%}







