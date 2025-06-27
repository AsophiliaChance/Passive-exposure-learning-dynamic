%% 准备数据并进行LMM分析和残差诊断（Amplitude & Latency）
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
        % 提取残差
        res = residuals(lme);
        
        fprintf('\nComponent: %s\n', components{i});
        
        %% 1. 残差正态性检验
        %fprintf('- 残差正态性检验(Shapiro-Wilk)：\n');
        % 替换为lillietest (Lilliefors检验, 类似Shapiro-Wilk)
        [H_lillie,p_lillie] = lillietest(res);
        
        if H_lillie == 0
            %fprintf('  正态性良好 (p = %.4f)\n', p_lillie);
        else
            fprintf('  非正态残差 (p = %.4f)\n', p_lillie);
        end
        
        
        %% 2. 方差齐性检验（Levene）
        fprintf('- 方差齐性检验(Levene)：\n');
        group = double(T.Days);
        p_levene = vartestn(res, group, 'TestType', 'LeveneAbsolute', 'Display', 'off');
        
        if p_levene > 0.05
             fprintf('  方差齐性成立 (p = %.4f)\n', p_levene);
        else
            fprintf('  存在异方差 (p = %.4f)\n', p_levene);
        end
        
        %% 3. 独立性检验（Durbin-Watson）
        fprintf('- 独立性检验(Durbin-Watson)：\n');
        % 将设计矩阵设为常数项
        X = ones(size(res));
        
        % 获取排序索引：先按 days 升序，再按 subjects 升序
        [~, sort_idx] = sortrows([T.Days, T.Subjects]);
        
        % 重排残差
        res_sorted = res(sort_idx);
        
        % 正确调用 dwtest：获取 p 值 和 DW 统计量
        [pval, dwStat] = dwtest(res_sorted, X);
        
        % 输出
        fprintf('  DW统计量 = %.4f\n', dwStat);
        fprintf('  p值 = %.4f\n', pval);
        
        % 判断是否存在自相关
        if dwStat > 1.5 && dwStat < 2.5
            fprintf('  残差未发现明显自相关\n');
        else
            fprintf('  可能存在自相关性\n');
        end
        
        
    end
end
