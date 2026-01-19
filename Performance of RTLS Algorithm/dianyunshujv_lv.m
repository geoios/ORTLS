%% RTLS_Eqn 实测数据测试脚本
% 从Excel文件读取实测数据并进行RTLS估计

clear; clc; close all;

%% 从Excel文件读取数据
% Excel文件路径（请根据实际情况修改）
excel_file = 'C:\Users\24159\Desktop\synthetic_data.xlsx';

% 检查文件是否存在
if ~exist(excel_file, 'file')
    error('文件不存在: %s\n请检查文件路径！', excel_file);
end

% 读取Excel文件
try
    % 读取工作表信息
    [~, sheet_names] = xlsfinfo(excel_file);
    
    % 读取第一个工作表（数据）
    if exist('readtable', 'file')
        data = readtable(excel_file, 'Sheet', 1);
        data_matrix = table2array(data);
    else
        [data_matrix, ~, ~] = xlsread(excel_file, 1);
        data_matrix = data_matrix(~all(isnan(data_matrix), 2), ~all(isnan(data_matrix), 1));
    end
    
    [num_rows, num_cols] = size(data_matrix);
    
    % 解析数据格式
    if num_cols == 2
        % 线性回归: y = ax + b
        x_data = data_matrix(:, 1);
        y = data_matrix(:, 2);
        A = [x_data, ones(num_rows, 1)];
    elseif num_cols >= 3
        % 多变量回归
        x_data = data_matrix(:, 1:end-1);
        y = data_matrix(:, end);
        A = [x_data, ones(num_rows, 1)];
    else
        error('数据列数不足，至少需要2列');
    end
    
    [m, n] = size(A);
    
    % 尝试读取Q_c矩阵（第2个工作表）
    Q_c = [];
    if length(sheet_names) >= 2
        try
            if exist('readtable', 'file')
                Q_c_data = readtable(excel_file, 'Sheet', 2);
                Q_c = table2array(Q_c_data);
            else
                [Q_c, ~, ~] = xlsread(excel_file, 2);
                Q_c = Q_c(~all(isnan(Q_c), 2), ~all(isnan(Q_c), 1));
            end
        catch
        end
    end
    
    % 如果没有Q_c或维度不匹配，使用单位矩阵
    if isempty(Q_c) || size(Q_c, 1) ~= m || size(Q_c, 2) ~= m
        Q_c = eye(m);
    end
    
catch ME
    error('读取Excel文件失败: %s', ME.message);
end

%% RTLS估计
% 设置参数
options.k0 = 1.0;
options.k1 = 1.8;
options.max_iter = 20;
options.tol = 1e-4;
options.max_inner_iter = 3;

% 开始计算
time_start = tic;

try
    [x_est, sigma0_sq, Q_x, iter_info] = RTLS_Eqn(A, y, Q_c, [], options);
    time_rtls = toc(time_start);
    
    %% 显示结果
    fprintf('\n========== RTLS估计结果 ==========\n\n');
    
    % 1. 参数估计值及不确定度
    fprintf('参数估计值及不确定度:\n');
    param_names = cell(length(x_est), 1);
    for i = 1:length(x_est)
        % 标准不确定度 (k=1)
        u_standard = sqrt(sigma0_sq * Q_x(i, i));
        % 扩展不确定度 (k=2, 约95%置信水平)
        u_expanded = 2 * u_standard;
        
        if i == length(x_est)
            param_names{i} = '截距 b';
            fprintf('  截距 b = %.8f ± %.8f (k=1)\n', x_est(i), u_standard);
            fprintf('           %.8f ± %.8f (k=2, 95%%置信)\n', x_est(i), u_expanded);
        else
            param_names{i} = sprintf('斜率 k_%d', i);
            fprintf('  斜率 k = %.8f ± %.8f (k=1)\n', x_est(i), u_standard);
            fprintf('           %.8f ± %.8f (k=2, 95%%置信)\n', x_est(i), u_expanded);
        end
    end
    
    % 2. 单位权中误差
    fprintf('\n单位权中误差:\n');
    sigma0_rtls = sqrt(sigma0_sq);
    fprintf('  σ₀ = %.6f\n', sigma0_rtls);
    fprintf('  σ₀² = %.6f\n', sigma0_sq);
    
    % 3. 算法运行时间
    fprintf('\n算法运行时间:\n');
    fprintf('  总耗时: %.3f 秒\n', time_rtls);
    fprintf('  迭代次数: %d\n', iter_info.outer_iter);
    
    % 4. 粗差检测摘要
    if isfield(iter_info, 'rejected_idx') && isfield(iter_info, 'downweighted_idx')
        rejected_idx = iter_info.rejected_idx;
        downweighted_idx = iter_info.downweighted_idx;
        fprintf('\n粗差检测摘要:\n');
        fprintf('  总观测数: %d\n', m);
        fprintf('  剔除点数: %d (%.1f%%)\n', length(rejected_idx), length(rejected_idx)/m*100);
        fprintf('  降权点数: %d (%.1f%%)\n', length(downweighted_idx), length(downweighted_idx)/m*100);
    end
    
    fprintf('\n========================================\n');
    
catch ME
    fprintf('\n✗ 计算失败: %s\n', ME.message);
    rethrow(ME);
end

%% 可视化结果
try
    % 提取自变量
    if size(A, 2) == 2 && all(A(:, 2) == 1)
        x_plot = A(:, 1);
    else
        x_plot = (1:m)';
    end
    
    % 创建图形
    figure('Position', [100, 100, 1400, 600]);
    set(gcf, 'DefaultTextFontName', 'Microsoft YaHei');
    set(gcf, 'DefaultAxesFontName', 'Microsoft YaHei');
    
    % 子图1: 拟合结果
    subplot(1, 2, 1);
    hold on;
    
    % 绘制所有观测点
    scatter(x_plot, y, 30, [0.6, 0.6, 0.6], 'filled', 'DisplayName', '观测数据');
    
    % 绘制拟合线
    y_fitted = A * x_est;
    [x_sorted, idx] = sort(x_plot);
    plot(x_sorted, y_fitted(idx), 'b-', 'LineWidth', 2.5, 'DisplayName', 'RTLS拟合');
    
    % 标记剔除点（红色叉号）
    if isfield(iter_info, 'rejected_idx') && ~isempty(iter_info.rejected_idx)
        rejected_idx = iter_info.rejected_idx;
        scatter(x_plot(rejected_idx), y(rejected_idx), 100, [0.8, 0, 0], ...
            'x', 'LineWidth', 2.5, 'DisplayName', sprintf('剔除点 (%d)', length(rejected_idx)));
    end
    
    % 标记降权点（橙色圆圈）
    if isfield(iter_info, 'downweighted_idx') && ~isempty(iter_info.downweighted_idx)
        downweighted_idx = iter_info.downweighted_idx;
        scatter(x_plot(downweighted_idx), y(downweighted_idx), 80, [1, 0.5, 0], ...
            'o', 'LineWidth', 2, 'DisplayName', sprintf('降权点 (%d)', length(downweighted_idx)));
    end
    
    xlabel('自变量 x', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('因变量 y', 'FontSize', 12, 'FontWeight', 'bold');
    title('RTLS 拟合结果', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    box on;
    hold off;
    
    % 子图2: 标准化残差图
    subplot(1, 2, 2);
    hold on;
    
    % 计算标准化残差
    residuals = y - y_fitted;
    if isfield(iter_info, 'v_t_final')
        std_residuals = iter_info.v_t_final;
    else
        std_residuals = residuals / sqrt(sigma0_sq);
    end
    
    % 绘制所有标准化残差
    scatter(x_plot, std_residuals, 30, [0.6, 0.6, 0.6], 'filled', 'DisplayName', '标准化残差');
    
    % 绘制阈值线
    plot([min(x_plot), max(x_plot)], [options.k0, options.k0], 'g--', 'LineWidth', 2, 'DisplayName', sprintf('保守区阈值 k₀=%.1f', options.k0));
    plot([min(x_plot), max(x_plot)], [-options.k0, -options.k0], 'g--', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([min(x_plot), max(x_plot)], [options.k1, options.k1], 'r--', 'LineWidth', 2, 'DisplayName', sprintf('排除区阈值 k₁=%.1f', options.k1));
    plot([min(x_plot), max(x_plot)], [-options.k1, -options.k1], 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot([min(x_plot), max(x_plot)], [0, 0], 'k-', 'LineWidth', 1);
    
    % 标记剔除点
    if isfield(iter_info, 'rejected_idx') && ~isempty(iter_info.rejected_idx)
        rejected_idx = iter_info.rejected_idx;
        scatter(x_plot(rejected_idx), std_residuals(rejected_idx), 100, [0.8, 0, 0], ...
            'x', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    
    % 标记降权点
    if isfield(iter_info, 'downweighted_idx') && ~isempty(iter_info.downweighted_idx)
        downweighted_idx = iter_info.downweighted_idx;
        scatter(x_plot(downweighted_idx), std_residuals(downweighted_idx), 80, [1, 0.5, 0], ...
            'o', 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    
    xlabel('自变量 x', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('标准化残差', 'FontSize', 12, 'FontWeight', 'bold');
    title('标准化残差分布图', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    box on;
    hold off;
    
    % 添加总标题（兼容旧版本MATLAB）
    annotation('textbox', [0.35, 0.95, 0.3, 0.04], ...
        'String', 'RTLS 抗差估计结果', ...
        'FontSize', 16, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
        'FontName', 'Microsoft YaHei');
    
catch ME
    fprintf('✗ 可视化失败: %s\n', ME.message);
end

%% 保存结果到MAT文件
% 提取剔除点和降权点索引
if isfield(iter_info, 'rejected_idx')
    rtls_rejected_idx = iter_info.rejected_idx;
else
    rtls_rejected_idx = [];
end

if isfield(iter_info, 'downweighted_idx')
    rtls_downweighted_idx = iter_info.downweighted_idx;
else
    rtls_downweighted_idx = [];
end

% 保存RTLS结果到MAT文件
rtls_results = struct();
rtls_results.X = x_est;
rtls_results.sigma0_sq = sigma0_sq;
rtls_results.Q_x = Q_x;
rtls_results.rejected_idx = rtls_rejected_idx;
rtls_results.downweighted_idx = rtls_downweighted_idx;
rtls_results.iter_info = iter_info;
rtls_results.time = time_rtls;
rtls_results.method = 'RTLS (Robust TLS)';
rtls_results.A = A;
rtls_results.L = y;

% 保存到桌面
save_path = 'C:\Users\24159\Desktop\RTLS.mat';
save(save_path, 'rtls_results', '-v7.3');


