%% load files
clc;clear
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');
%%
load('peakdetection.mat')
%%
%% MATLAB �ű�����4��14��2��6���ݽ����ظ������������
data = amplitude.mmn;

% ��������Ƿ�������ֵ
if any(data(:) <= 0)
    data_shifted = data - min(data(:)) + 1;
else
    data_shifted = data;
end
% �����任
data_log = log(data_shifted);
%% ��̬�Լ���
% ����һ����ÿ�������������һ�У�����̬�Խ��м���
[h_norm1, p_norm1] = lillietest(data_log(:,1));
fprintf('����1�����任�����̬�Լ��� p ֵ: %f\n', p_norm1);

% ������������������������Ƿ������̬�ֲ�����������ƽ��������
[h_norm_all, p_norm_all] = lillietest(data_log(:));
fprintf('��������任����̬�Լ��� p ֵ: %f\n', p_norm_all);

%% �������Լ���
% ���� data_log ��ÿһ�д���һ���������
% ���������
[nRows, nCols] = size(data_log);
group = repmat(1:nCols, nRows, 1);  % ÿ��Ϊһ��
group = group(:);
data_vector = data_log(:);

% ʹ�� vartestn ���� Levene ���飨'LeveneAbsolute'��
[p_var, stats] = vartestn(data_vector, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
fprintf('�������Լ��飨Levene���飩 p ֵ: %f\n', p_var);
data=data_log;
%% ��������
data_2D = reshape(data, 14, []);  % ���Ϊ 14��48 ����

%% �����������
% ��������˳�����Ϊ Stimulus���ڲ�����Ϊ Day��Channel�������������
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

% �� data_2D ת��Ϊ table�����趨������
T = array2table(data_2D, 'VariableNames', varNames);

%% ���� within-design ��
% 48 ������������Ӧ 3 �����أ�
%   Stimulus��ǰ24������Ϊ�̼�1����24��Ϊ�̼�2
%   Day��ÿ���̼��£�ÿ6������Ϊһ�죨��4�죩
%   Channel��ÿһ���ڰ�˳��Ϊͨ��1��6
Stimulus = [ones(24,1); 2*ones(24,1)];            % 24 �� 1 ��� 24 �� 2
Day = [repelem((1:4)', 6); repelem((1:4)', 6)];      % ��ÿ���̼��飬ÿ�� 6 ������
Channel = repmat((1:6)', 8, 1);                      % 8 �����飨4���2�̼�����ÿ������ 6 ��ͨ��

WithinDesign = table(Stimulus, Day, Channel);

%% ����ظ�����ģ��
% ʹ�� fitrm �� 48 ������ֵ��Ϊ�����������Э����������ָ�� within-design ��Ϣ
rm = fitrm(T, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WithinDesign);

%% �ظ������������
% ָ��ģ���п��� Stimulus��Day��Channel ���佻������
ranovatbl = ranova(rm, 'WithinModel', 'Stimulus*Day*Channel');

% ��ʾ����������
disp(ranovatbl);
% �� Day ���ؽ��������Ƚ�
% ��ʽһ��ֱ�Ӵ�����������
compDay = multcompare(rm, 'Day');
disp(compDay);

%%
for stim = 1:2
    fprintf('Analyzing Day �� Channel at Stimulus %d\n', stim);
    
    % ��ȡ��ǰ Stimulus �����ݣ�14 ���� �� 4 �� �� 6 ͨ����
    tempData = squeeze(data(:, :,:, stim));  % ���ӦΪ 14��4��6
    
    % ��������Ϊ 14��24 ����14 �У����� �� 24 �У�4��6 ��ϣ�
    tempData_2D = reshape(tempData, 14, []);  

    % ���ɱ�������
    tempVarNames = cell(1, 24);
    idx = 1;
    for day = 1:4
        for ch = 1:6
            tempVarNames{idx} = sprintf('D%d_C%d', day, ch);
            idx = idx + 1;
        end
    end
    
    % ת��Ϊ table
    tempT = array2table(tempData_2D, 'VariableNames', tempVarNames);

    % ���� WithinDesign��Day �� Channel ��ƣ�
    Day = repelem((1:4)', 6);
    Channel = repmat((1:6)', 4, 1);
    WithinDesign = table(Day, Channel);

    % �����ظ������������ģ��
    tempRM = fitrm(tempT, sprintf('%s-%s ~ 1', tempVarNames{1}, tempVarNames{end}), 'WithinDesign', WithinDesign);

    % ���� ANOVA �Է��� Day �� Channel ��������
    tempRanova = ranova(tempRM, 'WithinModel', 'Day*Channel');
    disp(tempRanova);
    
    % �º���飨�����������������
    fprintf('Post-hoc test for Day �� Channel at Stimulus %d\n', stim);
    compDayChan = multcompare(tempRM, 'Day', 'By', 'Channel');
    disp(compDayChan);

    % ���ƽ�������ͼ

    % ȷ�� Day �� Channel Ϊ��ֵ������ categorical ����
Day = categorical(repelem((1:4)', 6));  
Channel = categorical(repmat((1:6)', 4, 1));

% �����ֵ��Ϊ y
y = mean(tempData_2D, 1)';

% ���� interactionplot
figure;
interactionplot(y, {Day, Channel}, 'varnames', {'Day', 'Channel'});
title('Day �� Channel Interaction');

    title(sprintf('Stimulus %d: Day �� Channel Interaction', stim));
end

%{
%% ���� T ��һ�����ʽ��ÿ�ж�Ӧһ�����ԣ�ÿ�ж�Ӧһ�ֲ�������
% ���������� S1_D1_C1, S1_D1_C2, ... S2_D4_C6

for stim = 1:2
    % 1. ��ȡ Stimulus = stim �������� (��4���6ͨ��=24��)
    stimPrefix = sprintf('S%d_', stim);  % ���� 'S1_' �� 'S2_'
    idxStim = startsWith(T.Properties.VariableNames, stimPrefix);
    T_stim = T(:, idxStim);
    varNames_stim = T_stim.Properties.VariableNames;

    % 2. ���� WithinDesign ��
    %    Day: 1~4��ÿ��ͨ���ظ�һ�� -> 4��6=24��
    %    Channel: 1~6��ÿ���ظ�һ��
    Day = repelem((1:4)', 6);         % 4���6ͨ�� = 24
    Channel = repmat((1:6)', 4, 1);   % 6ͨ����4�� = 24
    WithinDesign_stim = table(Day, Channel);

    % 3. ����ظ�����ģ��
    rm_stim = fitrm(T_stim, ...
        sprintf('%s-%s ~ 1', varNames_stim{1}, varNames_stim{end}), ...
        'WithinDesign', WithinDesign_stim);

    % 4. ֻ��� Day ���ж��رȽ�
    compDay = multcompare(rm_stim, 'Day');

    % ������
    fprintf('\n===== Stimulus = %d: Day multcompare result =====\n', stim);
    disp(compDay);
end


%%
%{ 
%���� T ���������ݱ���������ʽΪ 'Sx_Dy_Cz'
% ѭ����������ͨ����1��6���ʹ̼���1��2��
for ch = 1:6
    for stim = 1:2
        % ���쵱ǰ�̼���ͨ���£�4������ݱ�����
        varNames = cell(1, 4);
        for d = 1:4
            varNames{d} = sprintf('S%d_D%d_C%d', stim, d, ch);
        end

        % ��ȡ��ǰ��ϵ������Ӽ�
        T_subset = T(:, varNames);

        % ��������� Day ���ص� within-design ��
        Day = (1:4)';  % ����
        WithinDesign = table(Day);

        % �����ظ�����ģ�ͣ�ģ�͹�ʽ��ʹ����ȡ��4������
        rm_temp = fitrm(T_subset, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WithinDesign);

        % �� Day ���ؽ��ж��رȽ�
        compDay = multcompare(rm_temp, 'Day');

        % ��ʾ��ǰͨ���ʹ̼���ϵıȽϽ��
        fprintf('====================\n');
        fprintf('Channel %d, Stimulus %d �� Day compare result��\n', ch, stim);
        disp(compDay);
    end
end
%}
%%
%% ��̬�Լ�����
%{
% �� T ����ÿ���������������������� Lilliefors ��̬�Լ���
varNames = T.Properties.VariableNames;
fprintf('��̬�Լ�������Lilliefors ���飩��\n');
for i = 1:length(varNames)
    data = T.(varNames{i});
    [h, p] = lillietest(data);
    fprintf('%s: h = %d, p = %.4f\n', varNames{i}, h, p);
end


%% 2. �������Լ��飨�������ڶ���������ƣ�
% ������ж����������ݣ�����ʹ�� vartestn ���鲻ͬ��֮��ķ����Ƿ���ȡ�
% ���磬���� data ��һ�����ݣ�group �Ƕ�Ӧ�ķ��������
% data = [...];   % ��������
% group = [...];  % ��������������Ϣ
% [p_var, stats] = vartestn(data, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
% fprintf('Levene ����� p ֵΪ: %.4f\n', p_var);

%% 3. �ظ�������Ƶ������Լ��飨Sphericity��
% �����ظ����� ANOVA����Ҫ�����ˮƽ��Ĳ�ֵЭ�����Ƿ����������Լ��衣
% ʹ�� mauchly ������ rm ģ�ͽ��м��顣
mauchlyResults = mauchly(rm);
disp('Mauchly �����Լ�������');
disp(mauchlyResults);
% �����ͨ������ʾ W ͳ�����Ͷ�Ӧ�� p ֵ��
% �� p < 0.05 ��ܾ������Լ��裬��ʱ��������� Greenhouse-Geisser У����

%% 4. �����Լ���
% ���� T �ǿ��ʽ���ݱ�ÿ�д���һ�����ԣ�����Ϊ��ͬ�����Ĳ���ֵ
% ��� T ��û�� Subject �����������һ��
if ~ismember('Subject', T.Properties.VariableNames)
    T.Subject = (1:height(T))';
end

% ��ȡ���в������������ų� Subject��
measureVars = setdiff(T.Properties.VariableNames, {'Subject'});

% �� T ת��Ϊ����ʽ
T_long = stack(T, measureVars, 'NewDataVariable', 'Y', 'IndexVariable', 'Condition');

% �������Ի��ЧӦģ�ͣ�ֻ���������ЧӦ��
lme = fitlme(T_long, 'Y ~ 1 + (1|Subject)');

% ��ȡģ�Ͳв�
residuals = lme.Residuals.Raw;

% ���Ʋв�������ͼ
figure;
autocorr(residuals);
title('���Ի��ЧӦģ�Ͳв������ͼ');

%}

%}
