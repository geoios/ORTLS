clear all; clc; close all;

%% ========== 从Excel文件读取数据 ==========
file_path = 'C:\Users\24159\Desktop\synthetic_data.xlsx';
fprintf('正在读取Excel文件: %s\n', file_path);

% 尝试读取Excel文件
try
    % 读取数据，假设第一列是x观测值，第二列是y观测值
    data = xlsread(file_path);
    
    % 检查数据格式
    if size(data, 2) < 2
        error('Excel文件需要至少包含两列数据：x观测值和y观测值');
    end
    
    % 提取x和y观测值
    x_obs = data(:, 1);
    y_obs = data(:, 2);
    
    n = length(x_obs);  % 数据点数量
    fprintf('成功读取 %d 个数据点\n', n);
    
    % 检查数据中是否有NaN值
    nan_indices = find(isnan(x_obs) | isnan(y_obs));
    if ~isempty(nan_indices)
        fprintf('警告: 数据中包含 %d 个NaN值，将自动剔除\n', length(nan_indices));
        
        % 剔除NaN值
        valid_indices = ~(isnan(x_obs) | isnan(y_obs));
        x_obs = x_obs(valid_indices);
        y_obs = y_obs(valid_indices);
        n = length(x_obs);
        fprintf('剔除NaN后剩余 %d 个数据点\n', n);
    end
    
catch ME
    fprintf('读取Excel文件失败: %s\n', ME.message);
    fprintf('将使用默认的仿真数据\n');
    
    % 如果读取失败，使用默认仿真数据
    n = 20;  % 数据点数量
    a_true = 2;  % 真实斜率
    b_true = -1;  % 真实截距
    sigma_x = 0.1;  % x方向噪声标准差
    sigma_y = 0.1;  % y方向噪声标准差
    
    % 生成仿真数据
    rng(57);
    x_true = rand(n, 1) * 10;  % U(0, 10)
    x_obs = x_true + randn(n, 1) * sigma_x;
    y_calculated = x_obs * a_true + b_true;
    y_obs = y_calculated + randn(n, 1) * sigma_y;
    
    fprintf('使用仿真数据: n = %d, a_true = %.1f, b_true = %.1f\n', n, a_true, b_true);
end

%% ========== 参数设置 ==========
m = 2;  % 参数个数（a和b）
alpha = 0.02;  % 显著性水平

% 噪声标准差（需要根据实际情况调整）
sigma_x = 0.1;  % x方向噪声标准差
sigma_y = 0.1;  % y方向噪声标准差

% 真实参数（如果使用仿真数据则有真实值，否则未知）
if exist('a_true', 'var') && exist('b_true', 'var')
    fprintf('真实参数: a = %.1f, b = %.1f\n', a_true, b_true);
else
    fprintf('真实参数未知，将仅进行粗差探测\n');
    a_true = NaN;
    b_true = NaN;
end

fprintf('噪声标准差: σ_x = %.3f, σ_y = %.3f\n', sigma_x, sigma_y);
fprintf('数据点数量: %d\n', n);
fprintf('显著性水平: α = %.2f\n\n', alpha);

%% ========== 使用分量压缩法进行粗差探测（迭代版本）==========
fprintf('========== 使用分量压缩法进行粗差探测（迭代版本）==========\n');

% ========== 开始计时整个算法流程 ==========
fprintf('开始计时整个算法流程（参数设置 + 粗差探测 + 剔除后重新估计）...\n');
tic;

% 构建设计矩阵
A_obs = [x_obs, ones(n, 1)];  % 设计矩阵：[x, 1]
L_obs = y_obs;  % 观测值：y

% 设置权阵
py = 1/sigma_y^2 * ones(n, 1);  % L的权
px1 = 1/sigma_x^2 * ones(n, 1);  % A1的权
px2 = 1e14 * ones(n, 1);  % A2的权（常数项，近似无噪声）
P = [py, px1, px2]';  % 3 x n

% 计算临界值（F检验，自由度为1和n-m）
df1 = 1;
df2 = n - m;
F_critical = sqrt(finv(1 - alpha, df1, df2));

% 调用迭代粗差探测函数
fprintf('开始迭代粗差探测...\n');
[detected_component, w_tests_component, v_component, x_hat_component, iter_info] = detect_outlier_v_iterative(A_obs, L_obs, P, alpha);

% 显示粗差探测结果
detected_positions = find(detected_component == 1);
num_detected = length(detected_positions);

fprintf('\n========== 粗差探测结果（迭代式数据探测）==========\n');
fprintf('临界值 (F检验, α=%.2f): %.4f\n', alpha, F_critical);
fprintf('迭代次数: %d\n', iter_info.iter);
fprintf('检测到的粗差点数量: %d (%.1f%%)\n', num_detected, num_detected/n*100);
fprintf('剩余有效点数量: %d (%.1f%%)\n', iter_info.n_final, iter_info.n_final/n*100);

% 【诊断信息】检查是否过度剔除
if num_detected > n * 0.3  % 如果剔除超过30%的点
    fprintf('\n⚠️ 警告: 剔除了%.1f%%的点，可能存在过度剔除！\n', num_detected/n*100);
    fprintf('诊断和建议:\n');
    fprintf('  1. 当前参数设置:\n');
    fprintf('     sigma_y = %.4f, sigma_x = %.4f\n', sigma_y, sigma_x);
    fprintf('     alpha = %.2f, F_critical = %.4f\n', alpha, F_critical);
    fprintf('  2. 可能的原因:\n');
    fprintf('     - sigma值设置过小，导致阈值过于严格\n');
    fprintf('     - alpha值设置过大（当前%.2f），建议降低到0.001-0.01\n', alpha);
    fprintf('     - 数据本身的噪声水平高于设定的sigma\n');
    fprintf('  3. 建议的解决方案:\n');
    fprintf('     a) 增大sigma_x和sigma_y的值（如从0.1改为0.5或1.0）\n');
    fprintf('     b) 减小alpha的值（如从%.2f改为0.01或0.001）\n', alpha);
    fprintf('     c) 从数据估计实际噪声水平: 使用MAD估计\n');
    fprintf('  4. 快速测试:\n');
    fprintf('     将第60-64行改为: sigma_x = 0.5; sigma_y = 0.5; 或\n');
    fprintf('     将第60行改为: alpha = 0.01;\n\n');
end

if num_detected > 0
    fprintf('\n--- 粗差点详细信息 ---\n');
    fprintf('粗差点位置: ');
    for i = 1:min(20, num_detected)
        fprintf('%d ', detected_positions(i));
    end
    if num_detected > 20
        fprintf('... (还有 %d 个点)', num_detected - 20);
    end
    fprintf('\n');
    
    % 显示前10个粗差点的详细信息
    if num_detected > 0
        fprintf('\n前10个粗差点的详细信息:\n');
        fprintf('%-6s %-12s %-12s %-12s %-12s\n', '序号', 'x观测值', 'y观测值', '残差', '|w|值');
        fprintf('%s\n', repmat('-', 1, 60));
        for i = 1:min(10, num_detected)
            idx = detected_positions(i);
            fprintf('%-6d %-12.6f %-12.6f %-12.6f %-12.4f\n', ...
                idx, x_obs(idx), y_obs(idx), v_component(idx), abs(w_tests_component(idx)));
        end
        if num_detected > 10
            fprintf('... (还有 %d 个粗差点未显示)\n', num_detected - 10);
        end
    end
    
    % 统计粗差的分布情况
    fprintf('\n--- 粗差统计信息 ---\n');
    w_outliers = abs(w_tests_component(detected_positions));
    fprintf('|w|统计量范围: [%.4f, %.4f]\n', min(w_outliers), max(w_outliers));
    fprintf('|w|统计量平均值: %.4f\n', mean(w_outliers));
    fprintf('|w|统计量中位数: %.4f\n', median(w_outliers));
    
    residuals_outliers = abs(v_component(detected_positions));
    fprintf('粗差点残差范围: [%.6f, %.6f]\n', min(residuals_outliers), max(residuals_outliers));
    fprintf('粗差点残差平均值: %.6f\n', mean(residuals_outliers));
else
    fprintf('未检测到粗差点\n');
end

% 显示所有点的统计信息（对比）
fprintf('\n--- 所有数据点的统计信息 ---\n');
fprintf('总数据点数: %d\n', n);
fprintf('正常点数量: %d (%.1f%%)\n', n - num_detected, (n - num_detected)/n*100);
fprintf('粗差点数量: %d (%.1f%%)\n', num_detected, num_detected/n*100);

fprintf('\n|w|统计量总体分布:\n');
fprintf('  最小值: %.4f\n', min(abs(w_tests_component)));
fprintf('  最大值: %.4f\n', max(abs(w_tests_component)));
fprintf('  平均值: %.4f\n', mean(abs(w_tests_component)));
fprintf('  中位数: %.4f\n', median(abs(w_tests_component)));
fprintf('  标准差: %.4f\n', std(abs(w_tests_component)));

fprintf('\n残差总体分布:\n');
fprintf('  最小值: %.6f\n', min(abs(v_component)));
fprintf('  最大值: %.6f\n', max(abs(v_component)));
fprintf('  平均值: %.6f\n', mean(abs(v_component)));
fprintf('  中位数: %.6f\n', median(abs(v_component)));
fprintf('  标准差: %.6f\n', std(abs(v_component)));

% 计算正常点的统计信息
if num_detected > 0 && num_detected < n
    normal_idx = ~detected_component;
    fprintf('\n正常点的|w|统计量分布:\n');
    fprintf('  最小值: %.4f\n', min(abs(w_tests_component(normal_idx))));
    fprintf('  最大值: %.4f\n', max(abs(w_tests_component(normal_idx))));
    fprintf('  平均值: %.4f\n', mean(abs(w_tests_component(normal_idx))));
    fprintf('  中位数: %.4f\n', median(abs(w_tests_component(normal_idx))));
    
    fprintf('\n正常点的残差分布:\n');
    fprintf('  最小值: %.6f\n', min(abs(v_component(normal_idx))));
    fprintf('  最大值: %.6f\n', max(abs(v_component(normal_idx))));
    fprintf('  平均值: %.6f\n', mean(abs(v_component(normal_idx))));
    fprintf('  中位数: %.6f\n', median(abs(v_component(normal_idx))));
end

%% ========== 参数估计结果 ==========
fprintf('\n========== 参数估计结果 ==========\n');

% 计算单位权中误差 sigma_0
% 使用WTLS方法计算残差和单位权中误差
% 构建协方差矩阵
Q_y = diag(1 ./ py);
Q_A = zeros(n*2, n*2);
Q_A(1:n, 1:n) = diag(1 ./ px1);
Q_A(n+1:end, n+1:end) = diag(1 ./ px2);

% 计算Q_y_tilde
x_kron_T = kron(x_hat_component', eye(n));
x_kron = kron(x_hat_component, eye(n));
Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
Q_y_tilde_inv = inv(Q_y_tilde);

% 计算残差
e_hat = L_obs - A_obs * x_hat_component;

% 单位权方差
sigma_0_squared = (e_hat' * Q_y_tilde_inv * e_hat) / (n - m);
sigma_0 = sqrt(sigma_0_squared);

% 计算参数的协方差矩阵
vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
E_A_hat = reshape(vec_E_A, n, 2);
A_tilde = A_obs - E_A_hat;
Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);

% 参数的标准差（不确定度）
sigma_a = sqrt(sigma_0_squared * Q_x(1, 1));
sigma_b = sqrt(sigma_0_squared * Q_x(2, 2));

fprintf('估计参数: a = %.6f ± %.6f, b = %.6f ± %.6f\n', ...
    x_hat_component(1), sigma_a, x_hat_component(2), sigma_b);
fprintf('单位权中误差 σ_0 = %.6f\n', sigma_0);
fprintf('参数不确定度: σ_a = %.6f, σ_b = %.6f\n', sigma_a, sigma_b);
fprintf('参数相关系数: ρ_ab = %.6f\n', Q_x(1, 2) / sqrt(Q_x(1, 1) * Q_x(2, 2)));

if ~isnan(a_true) && ~isnan(b_true)
    fprintf('真实参数: a = %.6f, b = %.6f\n', a_true, b_true);
    fprintf('参数误差: Δa = %.6f, Δb = %.6f\n', ...
        x_hat_component(1) - a_true, x_hat_component(2) - b_true);
end

%% ========== 抗差剔除后重新估计参数 ==========
if num_detected > 0
    fprintf('\n========== 抗差剔除后重新估计参数 ==========\n');
    
    % 剔除检测到的粗差点
    valid_idx = find(detected_component == 0);
        n_valid = length(valid_idx);
        
    if n_valid >= m
            % 提取有效数据
        x_obs_clean = x_obs(valid_idx);
        y_obs_clean = y_obs(valid_idx);
        A_obs_clean = [x_obs_clean, ones(n_valid, 1)];
        L_obs_clean = y_obs_clean;
        
        % 重新构建权阵
        py_clean = 1/sigma_y^2 * ones(n_valid, 1);
        px1_clean = 1/sigma_x^2 * ones(n_valid, 1);
        px2_clean = 1e14 * ones(n_valid, 1);
        P_clean = [py_clean, px1_clean, px2_clean]';
        
        % 使用WTLS方法进行参数估计
        % 构建协方差矩阵
        Q_y_clean = diag(1 ./ py_clean);  % L的协方差矩阵
        Q_A_clean = zeros(n_valid*2, n_valid*2);
        Q_A_clean(1:n_valid, 1:n_valid) = diag(1 ./ px1_clean);  % x列的协方差
        Q_A_clean(n_valid+1:end, n_valid+1:end) = diag(1 ./ px2_clean);  % 常数列的协方差
        
        [x_hat_clean_component, sigma_0_clean, Q_x_clean] = WTLS_solver_with_stats(A_obs_clean, L_obs_clean, Q_y_clean, Q_A_clean);
        
        % 计算参数的标准差（不确定度）
        sigma_a_clean = sqrt(sigma_0_clean^2 * Q_x_clean(1, 1));
        sigma_b_clean = sqrt(sigma_0_clean^2 * Q_x_clean(2, 2));
        
        fprintf('剔除 %d 个粗差点后剩余 %d 个数据点\n', num_detected, n_valid);
        fprintf('重新估计参数: a = %.6f ± %.6f, b = %.6f ± %.6f\n', ...
            x_hat_clean_component(1), sigma_a_clean, x_hat_clean_component(2), sigma_b_clean);
        fprintf('单位权中误差 σ_0 = %.6f\n', sigma_0_clean);
        fprintf('参数不确定度: σ_a = %.6f, σ_b = %.6f\n', sigma_a_clean, sigma_b_clean);
        
        if ~isnan(a_true) && ~isnan(b_true)
            fprintf('参数误差: Δa = %.6f, Δb = %.6f\n', ...
                x_hat_clean_component(1) - a_true, x_hat_clean_component(2) - b_true);
        end
    else
        fprintf('警告: 剔除粗差点后剩余数据点不足，无法重新估计参数\n');
    end
end

% ========== 结束计时整个算法流程 ==========
time_total = toc;
fprintf('\n完整算法流程总运行时间: %.4f秒\n', time_total);

%% ========== 绘制图形结果 ==========
fprintf('\n========== 绘制图形结果 ==========\n');

% 设置图形字体
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    
% 图1：原始数据散点图和拟合线
figure('Position', [100, 100, 1200, 500]);
    
% 子图1：原始数据
subplot(1, 2, 1);
    hold on;
    grid on;

% 绘制所有数据点
scatter(x_obs, y_obs, 60, 'b', 'filled', 'DisplayName', '正常点');
if num_detected > 0
    scatter(x_obs(detected_positions), y_obs(detected_positions), 80, 'r', 'filled', 'DisplayName', '粗差点');
end

% 绘制拟合线
x_range = [min(x_obs), max(x_obs)];
y_fit = x_range * x_hat_component(1) + x_hat_component(2);
plot(x_range, y_fit, 'k-', 'LineWidth', 2, 'DisplayName', '拟合线 (所有数据)');

% 如果进行了抗差剔除，绘制剔除后的拟合线
if num_detected > 0 && n_valid >= m
    y_fit_clean = x_range * x_hat_clean_component(1) + x_hat_clean_component(2);
    plot(x_range, y_fit_clean, 'g--', 'LineWidth', 2, 'DisplayName', '拟合线 (剔除粗差后)');
end

xlabel('x观测值', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('y观测值', 'FontSize', 14, 'FontWeight', 'bold');
title('数据散点图及拟合线', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

% 子图2：残差图
subplot(1, 2, 2);
hold on;
grid on;

% 计算残差
residuals = y_obs - (x_obs * x_hat_component(1) + x_hat_component(2));

% 绘制残差
plot(1:n, residuals, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '残差');
plot([1, n], [0, 0], 'k-', 'LineWidth', 1, 'DisplayName', '零线');

% 标记粗差点
if num_detected > 0
    plot(detected_positions, residuals(detected_positions), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '粗差点');
end

xlabel('数据点序号', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('残差', 'FontSize', 14, 'FontWeight', 'bold');
title('残差图', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

% 图2：w检验统计量
figure('Position', [100, 100, 1000, 600]);
hold on;
grid on;

% 绘制w检验值
bar(1:n, abs(w_tests_component), 'FaceColor', [0.7, 0.7, 0.9], 'EdgeColor', 'none', 'DisplayName', '|w|检验值');

% 标记超过临界值的点
if num_detected > 0
    bar(detected_positions, abs(w_tests_component(detected_positions)), ...
        'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', '粗差点 (|w|>临界值)');
end

% 添加临界值线
plot([1, n], [F_critical, F_critical], 'r--', 'LineWidth', 2.5, 'DisplayName', sprintf('临界值=%.4f', F_critical));

xlabel('数据点序号', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('|w|检验统计量', 'FontSize', 14, 'FontWeight', 'bold');
title('w检验统计量', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
xlim([0, n+1]);
set(gca, 'FontSize', 12);
hold off;

%% ========== 保存结果到文件 ==========
fprintf('\n========== 保存结果到文件 ==========\n');

% 准备保存结果
results_table = table((1:n)', x_obs, y_obs, residuals, w_tests_component, ...
    detected_component, ...
    'VariableNames', {'序号', 'x观测值', 'y观测值', '残差', 'w检验值', '是否粗差'});

% 添加参数估计结果
params_table = table({'a'; 'b'}, [x_hat_component(1); x_hat_component(2)], ...
    'VariableNames', {'参数', '估计值'});

if num_detected > 0 && n_valid >= m
    params_table = [x_hat_clean_component(1); x_hat_clean_component(2)];
end

% 保存Jazaeri方法结果到MAT文件（供点云文件读取）
X = x_hat_component;  % 参数估计
detected = detected_component;  % 粗差检测结果
w_tests = w_tests_component;  % w检验统计量
residuals_v = v_component;  % 残差
iter_info_save = iter_info;  % 迭代信息
F_critical_save = F_critical;  % 临界值
time_save = time_total;  % 完整算法流程的运行时间（包括参数设置+粗差探测+重新估计）

% 保存到桌面的jaz.mat文件
jaz_save_path = 'C:\Users\24159\Desktop\jaz.mat';
save(jaz_save_path, 'X', 'detected', 'w_tests', 'residuals_v', 'iter_info_save', 'F_critical_save', ...
    'x_obs', 'y_obs', 'n', 'num_detected', 'alpha', 'time_save');
fprintf('✓ Jazaeri方法结果已保存到: %s\n', jaz_save_path);
fprintf('  包含变量: X, detected, w_tests, residuals_v, iter_info_save, F_critical_save, time_save\n\n');

% 保存到Excel文件
output_file = 'C:\Users\24159\Desktop\outlier_detection_results.xlsx';
try
    writetable(results_table, output_file, 'Sheet', '粗差检测结果');
    writetable(params_table, output_file, 'Sheet', '参数估计结果');
    fprintf('✓ 详细结果已保存到Excel: %s\n', output_file);
catch
    fprintf('⚠ 无法保存Excel文件，但MAT文件已保存\n');
end

%% ========== 显示总结信息 ==========
fprintf('\n========== 分析总结 ==========\n');
fprintf('1. 数据点总数: %d\n', n);
fprintf('2. 检测到粗差点: %d (%.1f%%)\n', num_detected, num_detected/n*100);
fprintf('3. 参数估计结果:\n');
fprintf('   使用所有数据: a = %.6f ± %.6f, b = %.6f ± %.6f\n', ...
    x_hat_component(1), sigma_a, x_hat_component(2), sigma_b);
fprintf('   单位权中误差: σ_0 = %.6f\n', sigma_0);

if num_detected > 0 && n_valid >= m
    fprintf('   剔除粗差后: a = %.6f ± %.6f, b = %.6f ± %.6f\n', ...
        x_hat_clean_component(1), sigma_a_clean, x_hat_clean_component(2), sigma_b_clean);
    fprintf('   单位权中误差: σ_0 = %.6f\n', sigma_0_clean);
end

fprintf('4. 临界值 (α=%.2f): %.4f\n', alpha, F_critical);
fprintf('5. 最大 |w| 值: %.4f\n', max(abs(w_tests_component)));
fprintf('6. 完整算法流程总运行时间: %.4f秒\n', time_total);
fprintf('   (包括: 参数设置 + 粗差探测 + 剔除后重新估计)\n');

fprintf('\n分析完成！\n');

%% ========== WTLS求解器函数（带统计信息）==========
function [x_hat, sigma_0, Q_x] = WTLS_solver_with_stats(A, y, Q_y, Q_A)
    % WTLS求解器 - 基于Jazaeri论文的加权总体最小二乘法（带统计信息）
    % 输入:
    %   A - 设计矩阵 [m x n]
    %   y - 观测向量 [m x 1]
    %   Q_y - y的协方差矩阵 [m x m]
    %   Q_A - vec(A)的协方差矩阵 [mn x mn]
    % 输出:
    %   x_hat - 参数估计 [n x 1]
    %   sigma_0 - 单位权中误差
    %   Q_x - 参数的协因数矩阵
    
    [m, n] = size(A);
    
    % 初始化：使用OLS作为起点
    x_hat = (A' * A) \ (A' * y);
    
    epsilon = 1e-10;
    max_iter = 50;
    
    fprintf('  WTLS迭代开始...\n');
    
    for iter = 1:max_iter
        % 估计残差
        e_hat = y - A * x_hat;
        
        % 计算 Q_y_tilde = Q_y + (x^T ⊗ I_m) Q_A (x ⊗ I_m)
        x_kron_T = kron(x_hat', eye(m));  % m × mn
        x_kron = kron(x_hat, eye(m));      % mn × m
        Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
        Q_y_tilde_inv = inv(Q_y_tilde);
        
        % 计算 E_A_hat = -vec^{-1}(Q_A (x ⊗ I_m) Q_y_tilde^{-1} e_hat)
        vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
        E_A_hat = reshape(vec_E_A, m, n);
        
        % 预测值
        A_tilde = A - E_A_hat;
        y_tilde = y - E_A_hat * x_hat;
        
        % 更新参数
        x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * y_tilde);
        
        % 收敛检查
        delta = norm(x_hat_new - x_hat);
        
        if mod(iter, 5) == 0 || iter == 1
            fprintf('    迭代 %d: a = %.6f, b = %.6f, delta = %.2e\n', ...
                iter, x_hat_new(1), x_hat_new(2), delta);
        end
        
        x_hat = x_hat_new;
        
        if delta < epsilon
            fprintf('    WTLS收敛 (迭代%d次)!\n', iter);
            break;
        end
    end
    
    if delta >= epsilon
        fprintf('    警告: WTLS未完全收敛，达到最大迭代次数\n');
    end
    
    % 计算统计信息
    e_hat = y - A * x_hat;
    x_kron_T = kron(x_hat', eye(m));
    x_kron = kron(x_hat, eye(m));
    Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
    Q_y_tilde_inv = inv(Q_y_tilde);
    
    vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
    E_A_hat = reshape(vec_E_A, m, n);
    A_tilde = A - E_A_hat;
    
    % 单位权方差
    sigma_0_squared = (e_hat' * Q_y_tilde_inv * e_hat) / (m - n);
    sigma_0 = sqrt(sigma_0_squared);
    
    % 参数的协因数矩阵
    Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);
end

%% ========== WTLS求解器函数 ==========
function [x_hat, iter_count] = WTLS_solver(A, y, Q_y, Q_A)
    % WTLS求解器 - 基于Jazaeri论文的加权总体最小二乘法
    % 输入:
    %   A - 设计矩阵 [m x n]
    %   y - 观测向量 [m x 1]
    %   Q_y - y的协方差矩阵 [m x m]
    %   Q_A - vec(A)的协方差矩阵 [mn x mn]
    % 输出:
    %   x_hat - 参数估计 [n x 1]
    %   iter_count - 迭代次数
    
    [m, n] = size(A);
    
    % 初始化：使用OLS作为起点
    x_hat = (A' * A) \ (A' * y);
    
    epsilon = 1e-10;
    max_iter = 50;
    
    fprintf('  WTLS迭代开始...\n');
    
    iter_count = 0;
    for iter = 1:max_iter
        iter_count = iter;
        % 估计残差
        e_hat = y - A * x_hat;
        
        % 计算 Q_y_tilde = Q_y + (x^T ⊗ I_m) Q_A (x ⊗ I_m)
        x_kron_T = kron(x_hat', eye(m));  % m × mn
        x_kron = kron(x_hat, eye(m));      % mn × m
        Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
        Q_y_tilde_inv = inv(Q_y_tilde);
        
        % 计算 E_A_hat = -vec^{-1}(Q_A (x ⊗ I_m) Q_y_tilde^{-1} e_hat)
        vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
        % vec按列堆叠，所以前m个元素是第1列，后m个是第2列，等等
        E_A_hat = reshape(vec_E_A, m, n);
        
        % 预测值
        A_tilde = A - E_A_hat;
        y_tilde = y - E_A_hat * x_hat;
        
        % 更新参数
        x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * y_tilde);
        
        % 收敛检查
        delta = norm(x_hat_new - x_hat);
        
        if mod(iter, 5) == 0 || iter == 1
            fprintf('    迭代 %d: a = %.6f, b = %.6f, delta = %.2e\n', ...
                iter, x_hat_new(1), x_hat_new(2), delta);
        end
        
        x_hat = x_hat_new;
        
        if delta < epsilon
            fprintf('    WTLS收敛 (迭代%d次)!\n', iter);
            break;
        end
    end
    
    if delta >= epsilon
        fprintf('    警告: WTLS未完全收敛，达到最大迭代次数\n');
    end
end

%% ========== 迭代式粗差探测函数 ==========
function [detected, w_tests, v, x_hat, results] = detect_outlier_v_iterative(A_obs, L_obs, P, alpha)
% DETECT_OUTLIER_V_ITERATIVE 迭代式粗差探测函数（分量压缩法）
% 每次迭代剔除所有检测到的粗差，直到没有粗差为止
%
% 输入参数：
%   A_obs - 设计矩阵 [n x m]
%   L_obs - 观测值向量 [n x 1]
%   P     - 权阵 [3 x n]，格式为 [py; px1; px2]'
%   alpha - 显著性水平（可选，默认0.05）
%
% 输出参数：
%   detected   - 粗差检测结果 [n x 1]，相对于原始数据
%   w_tests    - w检验统计量 [n x 1]，最终迭代的结果
%   v          - 残差向量 [n x 1]，最终迭代的结果
%   x_hat      - 参数估计值 [m x 1]，最终迭代的结果
%   results    - 结构体，包含迭代信息

% 参数检查
if nargin < 4
    alpha = 0.02;
end

% 获取维度信息
[n_total, m] = size(A_obs);

% 初始化
detected_total = false(n_total, 1);  % 相对于原始数据的检测结果
valid_idx = (1:n_total)';  % 当前有效数据的索引
iter_count = 0;
max_iter = 50;  % 最大迭代次数
inner_iter_total = 0;  % 累积内层TLS迭代次数

fprintf('  [迭代检测] 开始迭代粗差检测...\n');

% 迭代检测
while iter_count < max_iter
    iter_count = iter_count + 1;
    
    % 当前有效数据
    A_current = A_obs(valid_idx, :);
    L_current = L_obs(valid_idx);
    n_current = length(valid_idx);
    
    % 检查是否还有足够的数据点
    if n_current < m + 3
        fprintf('  [迭代检测] 第%d次: 有效点数不足，停止检测\n', iter_count);
        break;
    end
    
    % 提取当前权阵
    P_current = P(:, valid_idx);
    
    % ========== 步骤1: 参数估计（使用WTLS） ==========
    % 构建协方差矩阵
    py_current = P_current(1, :)';
    px1_current = P_current(2, :)';
    px2_current = P_current(3, :)';
    
    Q_y_current = diag(1 ./ py_current);
    Q_A_current = zeros(n_current*2, n_current*2);
    Q_A_current(1:n_current, 1:n_current) = diag(1 ./ px1_current);
    Q_A_current(n_current+1:end, n_current+1:end) = diag(1 ./ px2_current);
    
    % 调用WTLS求解器
    [x_hat, iter_tls] = WTLS_solver(A_current, L_current, Q_y_current, Q_A_current);
    inner_iter_total = inner_iter_total + iter_tls;  % 累积内层迭代次数
    
    % ========== 步骤2: 计算残差和w检验 ==========
    e_hat = L_current - A_current * x_hat;
    
    % 计算Q_y_tilde
    x_kron_T = kron(x_hat', eye(n_current));
    x_kron = kron(x_hat, eye(n_current));
    Q_y_tilde = Q_y_current + x_kron_T * Q_A_current * x_kron;
    Q_y_tilde_inv = inv(Q_y_tilde);
    
    % 计算设计矩阵的改正
    vec_E_A = -Q_A_current * x_kron * Q_y_tilde_inv * e_hat;
    E_A_hat = reshape(vec_E_A, n_current, 2);
    A_tilde = A_current - E_A_hat;
    
    % 单位权方差
    sigma_0_squared = (e_hat' * Q_y_tilde_inv * e_hat) / (n_current - m);
    
    % 参数的协因数矩阵
    Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);
    
    % 计算残差的协方差矩阵
    Q_v = Q_y_tilde - A_tilde * Q_x * A_tilde';
    
    % w检验统计量
    w_tests_current = e_hat ./ sqrt(sigma_0_squared * diag(Q_v));
    
    % ========== 步骤3: 计算临界值 ==========
    df1 = 1;
    df2 = n_current - m;
    F_critical = sqrt(finv(1 - alpha, df1, df2));
    
    % ========== 步骤4: 找出所有超过阈值的粗差点 ==========
    detected_current = abs(w_tests_current) > F_critical;
    n_outliers_current = sum(detected_current);
    
    % 判断是否检测到粗差
    if n_outliers_current > 0
        % 找到粗差点在原始数据中的索引
        outlier_idx_local = find(detected_current);
        outlier_idx_global = valid_idx(outlier_idx_local);
        
        % 获取这些粗差点的w值
        w_outliers = w_tests_current(outlier_idx_local);
        
        fprintf('  [迭代检测] 第%d次: 检测到%d个粗差点, 剔除它们\n', ...
            iter_count, n_outliers_current);
        fprintf('    粗差点索引: %s\n', mat2str(outlier_idx_global'));
        fprintf('    对应|w|值: [%.3f-%.3f], 阈值=%.3f\n', ...
            min(abs(w_outliers)), max(abs(w_outliers)), F_critical);
        
        % 标记为粗差
        detected_total(outlier_idx_global) = true;
        
        % 从有效数据中移除所有检测到的粗差点
        valid_idx(outlier_idx_local) = [];
    else
        % 没有检测到更多粗差，停止迭代
        fprintf('  [迭代检测] 第%d次: 未检测到粗差 (max|w|=%.4f < %.4f), 停止检测\n', ...
            iter_count, max(abs(w_tests_current)), F_critical);
        break;
    end
end

% 最终结果
fprintf('  [迭代检测] 迭代结束: 共进行%d次迭代，检测到%d个粗差点\n', ...
    iter_count, sum(detected_total));

% 使用最终有效数据重新估计参数
A_final = A_obs(~detected_total, :);
L_final = L_obs(~detected_total);
P_final = P(:, ~detected_total);

py_final = P_final(1, :)';
px1_final = P_final(2, :)';
px2_final = P_final(3, :)';

Q_y_final = diag(1 ./ py_final);
Q_A_final = zeros(length(py_final)*2, length(py_final)*2);
Q_A_final(1:length(py_final), 1:length(py_final)) = diag(1 ./ px1_final);
Q_A_final(length(py_final)+1:end, length(py_final)+1:end) = diag(1 ./ px2_final);

[x_hat_final, iter_tls_final] = WTLS_solver(A_final, L_final, Q_y_final, Q_A_final);
inner_iter_total = inner_iter_total + iter_tls_final;  % 累积最终估计的内层迭代次数

% 计算所有点（包括粗差点）相对于最终模型的残差和w统计量
v_all = L_obs - A_obs * x_hat_final;
w_tests_all = zeros(n_total, 1);

% 对于非粗差点，计算准确的w统计量
valid_final = ~detected_total;
if sum(valid_final) >= m
    % 重新计算最终模型的统计信息
    e_hat_final = L_final - A_final * x_hat_final;
    x_kron_T_final = kron(x_hat_final', eye(sum(valid_final)));
    x_kron_final = kron(x_hat_final, eye(sum(valid_final)));
    Q_y_tilde_final = Q_y_final + x_kron_T_final * Q_A_final * x_kron_final;
    Q_y_tilde_inv_final = inv(Q_y_tilde_final);
    
    vec_E_A_final = -Q_A_final * x_kron_final * Q_y_tilde_inv_final * e_hat_final;
    E_A_hat_final = reshape(vec_E_A_final, sum(valid_final), 2);
    A_tilde_final = A_final - E_A_hat_final;
    
    sigma_0_squared_final = (e_hat_final' * Q_y_tilde_inv_final * e_hat_final) / (sum(valid_final) - m);
    Q_x_final = inv(A_tilde_final' * Q_y_tilde_inv_final * A_tilde_final);
    Q_v_final = Q_y_tilde_final - A_tilde_final * Q_x_final * A_tilde_final';
    
    w_tests_all(valid_final) = v_all(valid_final) ./ sqrt(sigma_0_squared_final * diag(Q_v_final));
    
    % 对于粗差点，给一个简化的w统计量估计
    sigma_robust = sqrt(sigma_0_squared_final) * median(sqrt(diag(Q_v_final)));
    w_tests_all(detected_total) = v_all(detected_total) / sigma_robust;
end

% 组织输出结果
detected = detected_total;
w_tests = w_tests_all;
v = v_all;
x_hat = x_hat_final;

if nargout > 4
    results = struct();
    results.iter = iter_count;
    results.inner_iter_total = inner_iter_total;  % 添加内层迭代总数
    results.n_outliers = sum(detected_total);
    results.outlier_indices = find(detected_total);
    results.n_final = sum(~detected_total);
end

end