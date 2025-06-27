%% ׼�����ݲ�����LMM�����Ͳв���ϣ�Amplitude & Latency��
clc;
components = {'Small MMN','Large MMN','Small P3a','Large P3a'};

fprintf('--- Linear Mixed Model Diagnostics ---\n\n');

for dataType = {amplitude_diff}%, latency_diff}
    if isequal(dataType{1}, amplitude_diff)
        fprintf('\n====== Amplitude ======\n');
        data = cat(3, amplitude_diff.mmn, amplitude_diff.p3);
    else
        fprintf('\n====== Latency ======\n');
        data = cat(3, latency_diff.mmn, latency_diff.p3);
    end
    
    for i = 3%1:4
        Days = repmat((1:4)',14,1);
        Subjects = repelem((1:14)',4);
        Y = reshape(data(:,:,i),[],1);
        

        T = table(categorical(Subjects),categorical(Days),Y,...
            'VariableNames',{'Subjects','Days','Y'});
       
        lme = fitlme(T, 'Y ~ Days + (Days|Subjects)');
        % ��ȡ�в�
        res = residuals(lme);
        
        fprintf('\nComponent: %s\n', components{i});
        
        %% 1. �в���̬�Լ���
        %fprintf('- �в���̬�Լ���(Shapiro-Wilk)��\n');
        % �滻Ϊlillietest (Lilliefors����, ����Shapiro-Wilk)
        [H_lillie,p_lillie] = lillietest(res);
        
        if H_lillie == 0
            %fprintf('  ��̬������ (p = %.4f)\n', p_lillie);
        else
            fprintf('  ����̬�в� (p = %.4f)\n', p_lillie);
        end
        
        
        %% 2. �������Լ��飨Levene��
        fprintf('- �������Լ���(Levene)��\n');
        group = double(T.Days);
        p_levene = vartestn(res, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
        
        if p_levene > 0.05
             fprintf('  �������Գ��� (p = %.4f)\n', p_levene);
        else
            fprintf('  �����췽�� (p = %.4f)\n', p_levene);
        end
        
        %% 3. �����Լ��飨Durbin-Watson��
        fprintf('- �����Լ���(Durbin-Watson)��\n');
        % ����ƾ�����Ϊ������
        X = ones(size(res));
        
        % ��ȡ�����������Ȱ� days �����ٰ� subjects ����
        [~, sort_idx] = sortrows([T.Days, T.Subjects]);
        
        % ���Ųв�
        res_sorted = res(sort_idx);
        
        % ��ȷ���� dwtest����ȡ p ֵ �� DW ͳ����
        [pval, dwStat] = dwtest(res_sorted, X);
        
        % ���
        fprintf('  DWͳ���� = %.4f\n', dwStat);
        fprintf('  pֵ = %.4f\n', pval);
        
        % �ж��Ƿ���������
        if dwStat > 1.5 && dwStat < 2.5
            fprintf('  �в�δ�������������\n');
        else
            fprintf('  ���ܴ����������\n');
        end
        
        
    end
end
