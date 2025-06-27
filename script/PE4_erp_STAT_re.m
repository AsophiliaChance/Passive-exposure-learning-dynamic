%% load files
clc;clear
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');
%%
load('peakdetection.mat')
%%
%% MATLAB 脚本：对4×14×2×6数据进行重复测量方差分析
data = amplitude.mmn;

% 检查数据是否存在零或负值
if any(data(:) <= 0)
    data_shifted = data - min(data(:)) + 1;
else
    data_shifted = data;
end
% 对数变换
data_log = log(data_shifted);
%% 正态性检验
% 方法一：对每个变量（例如第一列）的正态性进行检验
[h_norm1, p_norm1] = lillietest(data_log(:,1));
fprintf('变量1对数变换后的正态性检验 p 值: %f\n', p_norm1);

% 方法二：整体检验所有数据是否符合正态分布（将矩阵拉平成向量）
[h_norm_all, p_norm_all] = lillietest(data_log(:));
fprintf('整体对数变换后正态性检验 p 值: %f\n', p_norm_all);

%% 方差齐性检验
% 假设 data_log 的每一列代表一个组的数据
% 构造组变量
[nRows, nCols] = size(data_log);
group = repmat(1:nCols, nRows, 1);  % 每列为一组
group = group(:);
data_vector = data_log(:);

% 使用 vartestn 进行 Levene 检验（'LeveneAbsolute'）
[p_var, stats] = vartestn(data_vector, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
fprintf('方差齐性检验（Levene检验） p 值: %f\n', p_var);
data=data_log;
%% 数据重排
data_2D = reshape(data, 14, []);  % 结果为 14×48 矩阵

%% 构造变量名称
% 根据条件顺序（外层为 Stimulus，内层依次为 Day、Channel）构造变量名称
varNames = cell(1, 48);
idx = 1;
for stim = 1:2
    for day = 1:4
        for ch = 1:6
            varNames{idx} = sprintf('S%d_D%d_C%d', stim, day, ch);
            idx = idx + 1;
        end
    end
end

% 将 data_2D 转换为 table，并设定变量名
T = array2table(data_2D, 'VariableNames', varNames);

%% 构造 within-design 表
% 48 个测量条件对应 3 个因素：
%   Stimulus：前24个测量为刺激1，后24个为刺激2
%   Day：每个刺激下，每6个测量为一天（共4天）
%   Channel：每一天内按顺序为通道1～6
Stimulus = [ones(24,1); 2*ones(24,1)];            % 24 个 1 后接 24 个 2
Day = [repelem((1:4)', 6); repelem((1:4)', 6)];      % 对每个刺激组，每天 6 个测量
Channel = repmat((1:6)', 8, 1);                      % 8 个区块（4天×2刺激），每个区块 6 个通道

WithinDesign = table(Stimulus, Day, Channel);

%% 拟合重复测量模型
% 使用 fitrm 将 48 个测量值作为因变量（不含协变量），并指定 within-design 信息
rm = fitrm(T, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WithinDesign);

%% 重复测量方差分析
% 指定模型中考虑 Stimulus、Day、Channel 及其交互作用
ranovatbl = ranova(rm, 'WithinModel', 'Stimulus*Day*Channel');

% 显示方差分析结果
disp(ranovatbl);
% 对 Day 因素进行两两比较
% 方式一：直接传入因子名称
compDay = multcompare(rm, 'Day');
disp(compDay);

%%
for stim = 1:2
    fprintf('Analyzing Day × Channel at Stimulus %d\n', stim);
    
    % 提取当前 Stimulus 的数据（14 被试 × 4 天 × 6 通道）
    tempData = squeeze(data(:, :,:, stim));  % 结果应为 14×4×6
    
    % 重新排列为 14×24 矩阵（14 行，被试 × 24 列，4×6 组合）
    tempData_2D = reshape(tempData, 14, []);  

    % 生成变量名称
    tempVarNames = cell(1, 24);
    idx = 1;
    for day = 1:4
        for ch = 1:6
            tempVarNames{idx} = sprintf('D%d_C%d', day, ch);
            idx = idx + 1;
        end
    end
    
    % 转换为 table
    tempT = array2table(tempData_2D, 'VariableNames', tempVarNames);

    % 构造 WithinDesign（Day × Channel 设计）
    Day = repelem((1:4)', 6);
    Channel = repmat((1:6)', 4, 1);
    WithinDesign = table(Day, Channel);

    % 适配重复测量方差分析模型
    tempRM = fitrm(tempT, sprintf('%s-%s ~ 1', tempVarNames{1}, tempVarNames{end}), 'WithinDesign', WithinDesign);

    % 运行 ANOVA 以分析 Day × Channel 交互作用
    tempRanova = ranova(tempRM, 'WithinModel', 'Day*Channel');
    disp(tempRanova);
    
    % 事后检验（如果交互作用显著）
    fprintf('Post-hoc test for Day × Channel at Stimulus %d\n', stim);
    compDayChan = multcompare(tempRM, 'Day', 'By', 'Channel');
    disp(compDayChan);

    % 绘制交互作用图

    % 确保 Day 和 Channel 为数值向量或 categorical 变量
Day = categorical(repelem((1:4)', 6));  
Channel = categorical(repmat((1:6)', 4, 1));

% 计算均值作为 y
y = mean(tempData_2D, 1)';

% 调用 interactionplot
figure;
interactionplot(y, {Day, Channel}, 'varnames', {'Day', 'Channel'});
title('Day × Channel Interaction');

    title(sprintf('Stimulus %d: Day × Channel Interaction', stim));
end

%{
%% 假设 T 是一个宽格式表，每行对应一个被试，每列对应一种测量条件
% 变量名形如 S1_D1_C1, S1_D1_C2, ... S2_D4_C6

for stim = 1:2
    % 1. 提取 Stimulus = stim 的所有列 (共4天×6通道=24列)
    stimPrefix = sprintf('S%d_', stim);  % 例如 'S1_' 或 'S2_'
    idxStim = startsWith(T.Properties.VariableNames, stimPrefix);
    T_stim = T(:, idxStim);
    varNames_stim = T_stim.Properties.VariableNames;

    % 2. 构造 WithinDesign 表
    %    Day: 1~4，每个通道重复一次 -> 4×6=24行
    %    Channel: 1~6，每天重复一次
    Day = repelem((1:4)', 6);         % 4天×6通道 = 24
    Channel = repmat((1:6)', 4, 1);   % 6通道×4天 = 24
    WithinDesign_stim = table(Day, Channel);

    % 3. 拟合重复测量模型
    rm_stim = fitrm(T_stim, ...
        sprintf('%s-%s ~ 1', varNames_stim{1}, varNames_stim{end}), ...
        'WithinDesign', WithinDesign_stim);

    % 4. 只针对 Day 进行多重比较
    compDay = multcompare(rm_stim, 'Day');

    % 输出结果
    fprintf('\n===== Stimulus = %d: Day multcompare result =====\n', stim);
    disp(compDay);
end


%%
%{ 
%假设 T 是完整数据表，变量名格式为 'Sx_Dy_Cz'
% 循环遍历所有通道（1到6）和刺激（1到2）
for ch = 1:6
    for stim = 1:2
        % 构造当前刺激和通道下，4天的数据变量名
        varNames = cell(1, 4);
        for d = 1:4
            varNames{d} = sprintf('S%d_D%d_C%d', stim, d, ch);
        end

        % 提取当前组合的数据子集
        T_subset = T(:, varNames);

        % 构造仅包含 Day 因素的 within-design 表
        Day = (1:4)';  % 四天
        WithinDesign = table(Day);

        % 建立重复测量模型，模型公式中使用提取的4个变量
        rm_temp = fitrm(T_subset, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WithinDesign);

        % 对 Day 因素进行多重比较
        compDay = multcompare(rm_temp, 'Day');

        % 显示当前通道和刺激组合的比较结果
        fprintf('====================\n');
        fprintf('Channel %d, Stimulus %d 的 Day compare result：\n', ch, stim);
        disp(compDay);
    end
end
%}
%%
%% 正态性检验结果
%{
% 对 T 表中每个变量（测量条件）进行 Lilliefors 正态性检验
varNames = T.Properties.VariableNames;
fprintf('正态性检验结果（Lilliefors 检验）：\n');
for i = 1:length(varNames)
    data = T.(varNames{i});
    [h, p] = lillietest(data);
    fprintf('%s: h = %d, p = %.4f\n', varNames{i}, h, p);
end


%% 2. 方差齐性检验（仅适用于独立样本设计）
% 如果你有独立样本数据，可以使用 vartestn 检验不同组之间的方差是否相等。
% 例如，假设 data 是一列数据，group 是对应的分组变量：
% data = [...];   % 数据向量
% group = [...];  % 分组变量，类别信息
% [p_var, stats] = vartestn(data, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
% fprintf('Levene 检验的 p 值为: %.4f\n', p_var);

%% 3. 重复测量设计的球形性检验（Sphericity）
% 对于重复测量 ANOVA，需要检验各水平间的差值协方差是否满足球形性假设。
% 使用 mauchly 函数对 rm 模型进行检验。
mauchlyResults = mauchly(rm);
disp('Mauchly 球形性检验结果：');
disp(mauchlyResults);
% 结果中通常会显示 W 统计量和对应的 p 值，
% 若 p < 0.05 则拒绝球形性假设，此时建议采用如 Greenhouse-Geisser 校正。

%% 4. 独立性检验
% 假设 T 是宽格式数据表，每行代表一个被试，各列为不同条件的测量值
% 如果 T 中没有 Subject 变量，则添加一个
if ~ismember('Subject', T.Properties.VariableNames)
    T.Subject = (1:height(T))';
end

% 获取所有测量变量名（排除 Subject）
measureVars = setdiff(T.Properties.VariableNames, {'Subject'});

% 将 T 转换为长格式
T_long = stack(T, measureVars, 'NewDataVariable', 'Y', 'IndexVariable', 'Condition');

% 建立线性混合效应模型（只含随机被试效应）
lme = fitlme(T_long, 'Y ~ 1 + (1|Subject)');

% 提取模型残差
residuals = lme.Residuals.Raw;

% 绘制残差的自相关图
figure;
autocorr(residuals);
title('线性混合效应模型残差自相关图');

%}

%}
