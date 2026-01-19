%% 随机粗差收敛性诊断实验目前来说最中版本10.11最终版本论文实验
% 目的：测试随机粗差设计下，分方向方法和总体方法的收敛性
% 诊断哪个方法会出现收敛问题

clear; clc; close all;
rng(48);%39  47
fprintf('========== 随机粗差收敛性诊断 ==========\n');
fprintf('测试随机粗差下两种方法的收敛性能\n\n');

%% 实验参数设置
n = 20;  % 样本量
true_params = [2; -1];  % 真实参数 y = 2x - 3
outlier_ratio_3sigma = 0.10;  % 10%的粗差为3倍标准差
outlier_ratio_5sigma = 0.05;  % 5%的粗差为5倍标准差
outlier_ratio = outlier_ratio_3sigma + outlier_ratio_5sigma;  % 总粗差比例15%
outlier_magnitude_3sigma = 3;  % 3倍粗差
outlier_magnitude_5sigma = 5;  % 5倍粗差
max_experiments = 3000;  % 测试次数（少量测试便于观察）

% 存储收敛信息
convergence_results = struct();
convergence_results.directional_success = [];
convergence_results.overall_success = [];
convergence_results.standard_success = [];
convergence_results.directional_iters = [];
convergence_results.overall_iters = [];
convergence_results.standard_iters = [];
convergence_results.directional_params = [];
convergence_results.overall_params = [];
convergence_results.standard_params = [];
% 新增：WTLS数据探测法结果存储
convergence_results.wtls_success = [];
convergence_results.wtls_iters = [];
convergence_results.wtls_params = [];
convergence_results.wtls_errors = {};
convergence_results.directional_errors = {};
convergence_results.overall_errors = {};
convergence_results.standard_errors = {};
% 新增：RTLS方法结果存储
convergence_results.rtls_success = [];
convergence_results.rtls_iters = [];
convergence_results.rtls_params = [];
convergence_results.rtls_errors = {};
% 运行时间记录
convergence_results.directional_times = [];
convergence_results.overall_times = [];
convergence_results.wtls_times = [];
convergence_results.rtls_times = [];

fprintf('实验设置：\n');
fprintf('- 样本量: %d\n', n);
fprintf('- 真实参数: y = %.1fx + %.1f\n', true_params(1), true_params(2));
fprintf('- 粗差比例: %.0f%% (%.0f%%为3倍标准差, %.0f%%为5倍标准差)\n', ...
    outlier_ratio*100, outlier_ratio_3sigma*100, outlier_ratio_5sigma*100);
fprintf('- 测试次数: %d\n\n', max_experiments);

%% 开始诊断实验
for exp_idx = 1:max_experiments
    fprintf('========== 实验 %d/%d ==========\n', exp_idx, max_experiments);
    
    %% 数据生成（与原实验一致）
    % 步骤1：生成基础数据
    x = rand(n,1) * 10;  % U(0, 10) - 均匀分布，范围0到10
    
    % 步骤2：x真值加正常噪声（方差0.01，标准差0.1）
    x_noisy = x + randn(n, 1) * sqrt(0.0001);  % x + N(0, 0.01)

    % 步骤3：用加噪声的x计算y（y = 2x - 3）
    y_calculated = x_noisy * true_params(1) + true_params(2);

    % 步骤4：y加正常噪声（方差0.01，标准差0.1）
    y_noisy = y_calculated + randn(n, 1) * sqrt(0.0001);  % y + N(0, 0.01)

    % 步骤5：随机选择点加入异常值（10%为3倍，5%为5倍）
    % 正常噪声标准差：0.01
    sigma_x_noise = sqrt(0.0001);  % x方向噪声标准差
    sigma_y_noise =sqrt(0.0001);  % y方向噪声标准差
    
    n_outliers_3sigma = round(n * outlier_ratio_3sigma);  % 3倍粗差的点数
    n_outliers_5sigma = round(n * outlier_ratio_5sigma);  % 5倍粗差的点数
    n_outliers_total = n_outliers_3sigma + n_outliers_5sigma;
    
    if n_outliers_total > 0
        % 随机选择所有粗差点
        all_outlier_indices = randperm(n, n_outliers_total);
        % 前n_outliers_3sigma个点为3倍粗差
        outlier_indices_3sigma = all_outlier_indices(1:n_outliers_3sigma);
        % 剩余的点为5倍粗差
        if n_outliers_5sigma > 0
            outlier_indices_5sigma = all_outlier_indices(n_outliers_3sigma+1:end);
        else
            outlier_indices_5sigma = [];
        end
    else
        outlier_indices_3sigma = [];
        outlier_indices_5sigma = [];
        all_outlier_indices = [];
    end
    
    % 步骤6：添加固定大小粗差（10%为3倍，5%为5倍）
    fprintf('添加混合粗差（%.0f%%为3倍，%.0f%%为5倍标准差）...\n', ...
        outlier_ratio_3sigma*100, outlier_ratio_5sigma*100);
    
    % 计算粗差大小
    outlier_size_x_3sigma = outlier_magnitude_3sigma * sigma_x_noise;  % X方向3倍粗差
    outlier_size_y_3sigma = outlier_magnitude_3sigma * sigma_y_noise;  % Y方向3倍粗差
    outlier_size_x_5sigma = outlier_magnitude_5sigma * sigma_x_noise;  % X方向5倍粗差
    outlier_size_y_5sigma = outlier_magnitude_5sigma * sigma_y_noise;  % Y方向5倍粗差
    
    % 添加3倍粗差
    for i = 1:length(outlier_indices_3sigma)
        idx = outlier_indices_3sigma(i);
        % 随机方向（+1或-1），x和y方向相反确保点偏离拟合线
        sign_outlier = sign(randn(1));  % 随机方向：+1或-1
        if sign_outlier == 0  % 防止randn(1)恰好为0
            sign_outlier = 1;
        end
        x_noisy(idx) = x_noisy(idx) + sign_outlier * outlier_size_x_3sigma;
        y_noisy(idx) = y_noisy(idx) - sign_outlier * outlier_size_y_3sigma;  % 方向相反
    end
    
    % 添加5倍粗差
    for i = 1:length(outlier_indices_5sigma)
        idx = outlier_indices_5sigma(i);
        % 随机方向（+1或-1），x和y方向相反确保点偏离拟合线
        sign_outlier = sign(randn(1));  % 随机方向：+1或-1
        if sign_outlier == 0  % 防止randn(1)恰好为0
            sign_outlier = 1;
        end
        x_noisy(idx) = x_noisy(idx) + sign_outlier * outlier_size_x_5sigma;
        y_noisy(idx) = y_noisy(idx) - sign_outlier * outlier_size_y_5sigma;  % 方向相反
    end
    
    % 构建观测方程
    A = [x_noisy, ones(n,1)];
    L = y_noisy;
    P_initial = ones(3, n);
    P_initial(3, :) = 1e6;  % 常数项权重极大值
    
    %% 测试分方向残差法
    fprintf('\n--- 测试分方向残差法 ---\n');
    try
        tic;
        [X_dir, ~, iter_info_dir] = iterative_weight_optimization_with_timeout(A, L, P_initial, 5.0);
        time_dir = toc;
        
        if time_dir > 5.0
            % 超时情况
            convergence_results.directional_success(end+1) = 0;
            convergence_results.directional_iters(end+1) = NaN;
            convergence_results.directional_params(end+1, :) = [NaN, NaN];
            convergence_results.directional_errors{end+1} = '超时(5秒)';
            convergence_results.directional_times(end+1) = time_dir;
            
            fprintf('✗ 分方向方法超时(>5秒)\n');
            fprintf('  计算时间: %.3f秒\n', time_dir);
        else
            % 正常收敛
            convergence_results.directional_success(end+1) = 1;
            convergence_results.directional_iters(end+1) = iter_info_dir.total_iterations;
            convergence_results.directional_params(end+1, :) = X_dir';
            convergence_results.directional_errors{end+1} = '';
            convergence_results.directional_times(end+1) = time_dir;
            
            fprintf('✓ 分方向方法收敛成功\n');
            fprintf('  迭代次数: %d\n', iter_info_dir.total_iterations);
            fprintf('  计算时间: %.3f秒\n', time_dir);
            fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_dir(1), X_dir(2));
            fprintf('  参数误差: 斜率误差=%.6f, 截距误差=%.6f\n', X_dir(1)-true_params(1), X_dir(2)-true_params(2));
        end
        
    catch ME
        convergence_results.directional_success(end+1) = 0;
        convergence_results.directional_iters(end+1) = NaN;
        convergence_results.directional_params(end+1, :) = [NaN, NaN];
        convergence_results.directional_errors{end+1} = ME.message;
        convergence_results.directional_times(end+1) = NaN;
        
        fprintf('✗ 分方向方法收敛失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
    %% 测试总体残差法
    fprintf('\n--- 测试总体残差法 ---\n');
    try
        tic;
        [X_overall, ~, iter_info_overall] = overall_residual_weight_optimization_with_timeout(A, L, P_initial, 5.0);
        time_overall = toc;
        
        if time_overall > 5.0
            % 超时情况
            convergence_results.overall_success(end+1) = 0;
            convergence_results.overall_iters(end+1) = NaN;
            convergence_results.overall_params(end+1, :) = [NaN, NaN];
            convergence_results.overall_errors{end+1} = '超时(5秒)';
            convergence_results.overall_times(end+1) = time_overall;
            
            fprintf('✗ 总体方法超时(>5秒)\n');
            fprintf('  计算时间: %.3f秒\n', time_overall);
        else
            % 正常收敛
            convergence_results.overall_success(end+1) = 1;
            convergence_results.overall_iters(end+1) = iter_info_overall.total_iterations;
            convergence_results.overall_params(end+1, :) = X_overall';
            convergence_results.overall_errors{end+1} = '';
            convergence_results.overall_times(end+1) = time_overall;
            
            fprintf('✓ 总体方法收敛成功\n');
            fprintf('  迭代次数: %d\n', iter_info_overall.total_iterations);
            fprintf('  计算时间: %.3f秒\n', time_overall);
            fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_overall(1), X_overall(2));
            fprintf('  参数误差: 斜率误差=%.6f, 截距误差=%.6f\n', X_overall(1)-true_params(1), X_overall(2)-true_params(2));
        end
        
    catch ME
        convergence_results.overall_success(end+1) = 0;
        convergence_results.overall_iters(end+1) = NaN;
        convergence_results.overall_params(end+1, :) = [NaN, NaN];
        convergence_results.overall_errors{end+1} = ME.message;
        convergence_results.overall_times(end+1) = NaN;
        
        fprintf('✗ 总体方法收敛失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
    %% 测试标准TLS方法
    fprintf('\n--- 测试标准TLS方法 ---\n');
    try
        tic;
        X_standard = standard_TLS(A, L, P_initial);
        time_standard = toc;
        
        convergence_results.standard_success(end+1) = 1;
        convergence_results.standard_iters(end+1) = 1;  % 标准方法不需要迭代
        convergence_results.standard_params(end+1, :) = X_standard';
        convergence_results.standard_errors{end+1} = '';
        
        fprintf('✓ 标准TLS方法计算成功\n');
        fprintf('  计算时间: %.3f秒\n', time_standard);
        fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_standard(1), X_standard(2));
        fprintf('  参数误差: 斜率误差=%.6f, 截距误差=%.6f\n', X_standard(1)-true_params(1), X_standard(2)-true_params(2));
        
    catch ME
        convergence_results.standard_success(end+1) = 0;
        convergence_results.standard_iters(end+1) = NaN;
        convergence_results.standard_params(end+1, :) = [NaN, NaN];
        convergence_results.standard_errors{end+1} = ME.message;
        
        fprintf('✗ 标准TLS方法计算失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
    %% 测试WTLS数据探测法（Amiri-Simkooei & Jazaeri, 2013）
    fprintf('\n--- 测试WTLS数据探测法（含数据探测）---\n');
    try
        tic;
        [X_wtls, wtls_info] = wtls_data_snooping_estimate(A, L);
        time_wtls = toc;
        
        if wtls_info.converged
            convergence_results.wtls_success(end+1) = 1;
            convergence_results.wtls_iters(end+1) = wtls_info.total_iterations;
            convergence_results.wtls_params(end+1, :) = X_wtls';
            convergence_results.wtls_errors{end+1} = '';
            convergence_results.wtls_times(end+1) = time_wtls;
            
            fprintf('✓ WTLS方法收敛成功\n');
            fprintf('  迭代次数: %d（外层数据探测%d次，已移除粗差%d个）\n', ...
                    wtls_info.inner_iterations, wtls_info.snooping_iterations, wtls_info.outliers_removed);
            fprintf('  计算时间: %.3f秒\n', time_wtls);
            fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_wtls(1), X_wtls(2));
        else
            convergence_results.wtls_success(end+1) = 0;
            convergence_results.wtls_iters(end+1) = wtls_info.total_iterations;
            convergence_results.wtls_params(end+1, :) = [NaN, NaN];
            convergence_results.wtls_errors{end+1} = '未收敛';
            convergence_results.wtls_times(end+1) = time_wtls;
            
            fprintf('✗ WTLS方法未收敛\n');
        end
    catch ME
        convergence_results.wtls_success(end+1) = 0;
        convergence_results.wtls_iters(end+1) = NaN;
        convergence_results.wtls_params(end+1, :) = [NaN, NaN];
        convergence_results.wtls_errors{end+1} = ME.message;
        convergence_results.wtls_times(end+1) = NaN;
        
        fprintf('✗ WTLS方法计算失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
    %% 测试RTLS方法（Lv等人，抗差总体最小二乘）
    fprintf('\n--- 测试RTLS方法（Lv Com-Com）---\n');
    try
        tic;
        % 准备RTLS输入：协因数矩阵（使用单位矩阵）
        Q_c_rtls = eye(n);
        % RTLS参数设置
        options_rtls = struct();
        options_rtls.k0 = 1.5;
        options_rtls.k1 = 2.5;
        options_rtls.max_iter = 20;
        options_rtls.tol = 1e-2;
        options_rtls.max_inner_iter = 5;
        
        % 调用RTLS估计器
        [X_rtls, ~, ~, rtls_info] = RTLS_Eqn(A, L, Q_c_rtls, [], options_rtls);
        time_rtls = toc;
        
        if rtls_info.converged
            convergence_results.rtls_success(end+1) = 1;
            convergence_results.rtls_iters(end+1) = rtls_info.outer_iter;
            convergence_results.rtls_params(end+1, :) = X_rtls';
            convergence_results.rtls_errors{end+1} = '';
            convergence_results.rtls_times(end+1) = time_rtls;
            
            fprintf('✓ RTLS方法收敛成功\n');
            fprintf('  迭代次数: %d（外层迭代，已剔除粗差%d个）\n', ...
                    rtls_info.outer_iter, length(rtls_info.rejected_idx));
            fprintf('  计算时间: %.3f秒\n', time_rtls);
            fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_rtls(1), X_rtls(2));
        else
            convergence_results.rtls_success(end+1) = 0;
            convergence_results.rtls_iters(end+1) = rtls_info.outer_iter;
            convergence_results.rtls_params(end+1, :) = [NaN, NaN];
            convergence_results.rtls_errors{end+1} = '未收敛';
            convergence_results.rtls_times(end+1) = time_rtls;
            
            fprintf('✗ RTLS方法未收敛\n');
        end
    catch ME
        convergence_results.rtls_success(end+1) = 0;
        convergence_results.rtls_iters(end+1) = NaN;
        convergence_results.rtls_params(end+1, :) = [NaN, NaN];
        convergence_results.rtls_errors{end+1} = ME.message;
        convergence_results.rtls_times(end+1) = NaN;
        
        fprintf('✗ RTLS方法计算失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
    fprintf('\n');
end

%% 收敛性统计分析
fprintf('========== 收敛性统计分析 ==========\n\n');

% 成功率统计
dir_success_rate = sum(convergence_results.directional_success) / length(convergence_results.directional_success) * 100;
overall_success_rate = sum(convergence_results.overall_success) / length(convergence_results.overall_success) * 100;
standard_success_rate = sum(convergence_results.standard_success) / length(convergence_results.standard_success) * 100;
wtls_success_rate = sum(convergence_results.wtls_success) / length(convergence_results.wtls_success) * 100;
rtls_success_rate = sum(convergence_results.rtls_success) / length(convergence_results.rtls_success) * 100;

fprintf('【收敛成功率】\n');
fprintf('分方向残差法: %.1f%% (%d/%d)\n', dir_success_rate, sum(convergence_results.directional_success), length(convergence_results.directional_success));
fprintf('总体残差法:   %.1f%% (%d/%d)\n', overall_success_rate, sum(convergence_results.overall_success), length(convergence_results.overall_success));
fprintf('标准TLS法:    %.1f%% (%d/%d)\n', standard_success_rate, sum(convergence_results.standard_success), length(convergence_results.standard_success));
fprintf('WTLS数据探测: %.1f%% (%d/%d)\n', wtls_success_rate, sum(convergence_results.wtls_success), length(convergence_results.wtls_success));
fprintf('RTLS方法:     %.1f%% (%d/%d)\n', rtls_success_rate, sum(convergence_results.rtls_success), length(convergence_results.rtls_success));

% 迭代次数统计
valid_dir_iters = convergence_results.directional_iters(~isnan(convergence_results.directional_iters));
valid_overall_iters = convergence_results.overall_iters(~isnan(convergence_results.overall_iters));

if ~isempty(valid_dir_iters)
    fprintf('\n【迭代次数统计】\n');
    fprintf('分方向残差法: %.1f ± %.1f (范围: %d-%d)\n', mean(valid_dir_iters), std(valid_dir_iters), min(valid_dir_iters), max(valid_dir_iters));
end

if ~isempty(valid_overall_iters)
    fprintf('总体残差法:   %.1f ± %.1f (范围: %d-%d)\n', mean(valid_overall_iters), std(valid_overall_iters), min(valid_overall_iters), max(valid_overall_iters));
end

% 平均运行时间统计（过滤异常值）
valid_dir_times = convergence_results.directional_times(~isnan(convergence_results.directional_times));
valid_overall_times = convergence_results.overall_times(~isnan(convergence_results.overall_times));
valid_wtls_times = convergence_results.wtls_times(~isnan(convergence_results.wtls_times));
valid_rtls_times = convergence_results.rtls_times(~isnan(convergence_results.rtls_times));

% 过滤掉异常长的时间（使用中位数±3倍MAD方法，更稳健）
fprintf('\n【平均运行时间统计 (已过滤异常值)】\n');

if ~isempty(valid_dir_times)
    median_dir = median(valid_dir_times);
    mad_dir = median(abs(valid_dir_times - median_dir));
    % 过滤：保留中位数±3*1.4826*MAD范围内的值（1.4826是MAD到标准差的转换因子）
    filtered_dir_times = valid_dir_times(abs(valid_dir_times - median_dir) <= 3 * 1.4826 * mad_dir);
    fprintf('分方向残差法: %.4f ± %.4f秒 (范围: %.4f-%.4f秒) [使用%d/%d个样本]\n', ...
        mean(filtered_dir_times), std(filtered_dir_times), min(filtered_dir_times), max(filtered_dir_times), ...
        length(filtered_dir_times), length(valid_dir_times));
end

if ~isempty(valid_overall_times)
    median_overall = median(valid_overall_times);
    mad_overall = median(abs(valid_overall_times - median_overall));
    filtered_overall_times = valid_overall_times(abs(valid_overall_times - median_overall) <= 3 * 1.4826 * mad_overall);
    fprintf('总体残差法:   %.4f ± %.4f秒 (范围: %.4f-%.4f秒) [使用%d/%d个样本]\n', ...
        mean(filtered_overall_times), std(filtered_overall_times), min(filtered_overall_times), max(filtered_overall_times), ...
        length(filtered_overall_times), length(valid_overall_times));
end

if ~isempty(valid_wtls_times)
    median_wtls = median(valid_wtls_times);
    mad_wtls = median(abs(valid_wtls_times - median_wtls));
    filtered_wtls_times = valid_wtls_times(abs(valid_wtls_times - median_wtls) <= 3 * 1.4826 * mad_wtls);
    fprintf('WTLS数据探测: %.4f ± %.4f秒 (范围: %.4f-%.4f秒) [使用%d/%d个样本]\n', ...
        mean(filtered_wtls_times), std(filtered_wtls_times), min(filtered_wtls_times), max(filtered_wtls_times), ...
        length(filtered_wtls_times), length(valid_wtls_times));
end

if ~isempty(valid_rtls_times)
    median_rtls = median(valid_rtls_times);
    mad_rtls = median(abs(valid_rtls_times - median_rtls));
    filtered_rtls_times = valid_rtls_times(abs(valid_rtls_times - median_rtls) <= 3 * 1.4826 * mad_rtls);
    fprintf('RTLS方法:     %.4f ± %.4f秒 (范围: %.4f-%.4f秒) [使用%d/%d个样本]\n', ...
        mean(filtered_rtls_times), std(filtered_rtls_times), min(filtered_rtls_times), max(filtered_rtls_times), ...
        length(filtered_rtls_times), length(valid_rtls_times));
end

% 计算加速比（基于过滤后的数据）
if ~isempty(filtered_dir_times) && ~isempty(filtered_overall_times)
    speedup_overall_vs_dir = mean(filtered_dir_times) / mean(filtered_overall_times);
    fprintf('\n【加速比 (基于过滤后数据)】\n');
    fprintf('总体方法 vs 分方向方法: %.2fx\n', speedup_overall_vs_dir);
end
if ~isempty(filtered_wtls_times) && ~isempty(filtered_dir_times)
    speedup_wtls_vs_dir = mean(filtered_dir_times) / mean(filtered_wtls_times);
    fprintf('WTLS vs 分方向方法: %.2fx\n', speedup_wtls_vs_dir);
end
if ~isempty(filtered_overall_times) && ~isempty(filtered_wtls_times)
    speedup_wtls_vs_overall = mean(filtered_overall_times) / mean(filtered_wtls_times);
    fprintf('WTLS vs 总体方法: %.2fx\n', speedup_wtls_vs_overall);
end

% 参数估计精度统计
valid_dir_params = convergence_results.directional_params(convergence_results.directional_success == 1, :);
valid_overall_params = convergence_results.overall_params(convergence_results.overall_success == 1, :);

if ~isempty(valid_dir_params)
    dir_slope_rmse = sqrt(mean((valid_dir_params(:,1) - true_params(1)).^2));
    dir_intercept_rmse = sqrt(mean((valid_dir_params(:,2) - true_params(2)).^2));
    fprintf('\n【参数估计精度 - 分方向方法】\n');
    fprintf('斜率RMSE: %.6f, 截距RMSE: %.6f\n', dir_slope_rmse, dir_intercept_rmse);
end

if ~isempty(valid_overall_params)
    overall_slope_rmse = sqrt(mean((valid_overall_params(:,1) - true_params(1)).^2));
    overall_intercept_rmse = sqrt(mean((valid_overall_params(:,2) - true_params(2)).^2));
    fprintf('\n【参数估计精度 - 总体方法】\n');
    fprintf('斜率RMSE: %.6f, 截距RMSE: %.6f\n', overall_slope_rmse, overall_intercept_rmse);
end

% 失败原因分析
fprintf('\n【失败原因分析】\n');
dir_failures = find(convergence_results.directional_success == 0);
overall_failures = find(convergence_results.overall_success == 0);

if ~isempty(dir_failures)
    fprintf('分方向方法失败实验: %s\n', mat2str(dir_failures));
    for i = 1:length(dir_failures)
        fprintf('  实验%d: %s\n', dir_failures(i), convergence_results.directional_errors{dir_failures(i)});
    end
end

if ~isempty(overall_failures)
    fprintf('总体方法失败实验: %s\n', mat2str(overall_failures));
    for i = 1:length(overall_failures)
        fprintf('  实验%d: %s\n', overall_failures(i), convergence_results.overall_errors{overall_failures(i)});
    end
end

%% 参数分布可视化
if ~isempty(valid_dir_params) || ~isempty(valid_overall_params)
    fprintf('\n========== 参数分布可视化 ==========\n');
    
    % 设置全局字体为Times New Roman
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    
    % 创建图形窗口
    figure('Position', [100, 100, 1200, 800]);
    
    % 斜率分布对比
    subplot(2, 2, 1);
    hold on;
    
    % 计算统一的x范围（与抗差方法对比图保持一致）
    slope_min_kde = min([min(valid_dir_params(:,1)), min(valid_overall_params(:,1))]);
    slope_max_kde = max([max(valid_dir_params(:,1)), max(valid_overall_params(:,1))]);
    slope_range_kde = linspace(slope_min_kde, slope_max_kde, 200);
    
    if ~isempty(valid_dir_params)
        % 分方向方法斜率核密度估计（红色，与抗差方法对比图一致）
        [f_dir_slope, xi_dir_slope] = ksdensity(valid_dir_params(:,1), slope_range_kde);
        plot(xi_dir_slope, f_dir_slope, 'r-', 'LineWidth', 2, 'DisplayName', 'Full-Component');
    end
    
    if ~isempty(valid_overall_params)
        % 总体方法斜率核密度估计（蓝色，与抗差方法对比图一致）
        [f_overall_slope, xi_overall_slope] = ksdensity(valid_overall_params(:,1), slope_range_kde);
        plot(xi_overall_slope, f_overall_slope, 'b-', 'LineWidth', 2, 'DisplayName', 'Component-Compressed');
    end
% 新增：WTLS核密度曲线
successful_wtls = convergence_results.wtls_success == 1;
valid_wtls_params = convergence_results.wtls_params(successful_wtls, :);
if ~isempty(valid_wtls_params)
    [f_wtls_slope, xi_wtls_slope] = ksdensity(valid_wtls_params(:,1));
    plot(xi_wtls_slope, f_wtls_slope, 'm-', 'LineWidth', 2, 'DisplayName', 'WTLS (Data-snooping)');
end
    
    % 真实值标记
    xline(true_params(1), 'k--', 'LineWidth', 2, 'DisplayName', 'True Value');
    
    xlabel('Slope Estimate', 'FontSize', 15, 'FontName', 'Times New Roman');
    ylabel('Probability Density', 'FontSize', 15, 'FontName', 'Times New Roman');
    title('Slope Parameter Distribution (KDE)', 'FontSize', 17, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    legend('Location', 'eastoutside', 'FontSize', 15, 'FontName', 'Times New Roman');
    grid on;
    hold off;
    
    % 截距分布对比
    subplot(2, 2, 2);
    hold on;
    
    % 计算统一的x范围（与抗差方法对比图保持一致）
    intercept_min_kde = min([min(valid_dir_params(:,2)), min(valid_overall_params(:,2))]);
    intercept_max_kde = max([max(valid_dir_params(:,2)), max(valid_overall_params(:,2))]);
    intercept_range_kde = linspace(intercept_min_kde, intercept_max_kde, 200);
    
    if ~isempty(valid_dir_params)
        % 分方向方法截距核密度估计（红色，与抗差方法对比图一致）
        [f_dir_intercept, xi_dir_intercept] = ksdensity(valid_dir_params(:,2), intercept_range_kde);
        plot(xi_dir_intercept, f_dir_intercept, 'r-', 'LineWidth', 2, 'DisplayName', 'Full-Component');
    end
    
    if ~isempty(valid_overall_params)
        % 总体方法截距核密度估计（蓝色，与抗差方法对比图一致）
        [f_overall_intercept, xi_overall_intercept] = ksdensity(valid_overall_params(:,2), intercept_range_kde);
        plot(xi_overall_intercept, f_overall_intercept, 'b-', 'LineWidth', 2, 'DisplayName', 'Component-Compressed');
    end
% 新增：WTLS核密度曲线
if ~isempty(valid_wtls_params)
    [f_wtls_intercept, xi_wtls_intercept] = ksdensity(valid_wtls_params(:,2));
    plot(xi_wtls_intercept, f_wtls_intercept, 'm-', 'LineWidth', 2, 'DisplayName', 'WTLS (Data-snooping)');
end
    
    % 真实值标记
    xline(true_params(2), 'k--', 'LineWidth', 2, 'DisplayName', 'True Value');
    
    xlabel('Intercept Estimate', 'FontSize', 15, 'FontName', 'Times New Roman');
    ylabel('Probability Density', 'FontSize', 15, 'FontName', 'Times New Roman');
    title('Intercept Parameter Distribution (KDE)', 'FontSize', 17, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    legend('Location', 'best', 'FontSize', 15, 'FontName', 'Times New Roman');
    grid on;
    hold off;
    
    % 斜率直方图对比
    subplot(2, 2, 3);
    hold on;
    
    if ~isempty(valid_dir_params)
        histogram(valid_dir_params(:,1), 'Normalization', 'probability', ...
                 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                 'DisplayName', 'Full-Component');
    end
    
    if ~isempty(valid_overall_params)
        histogram(valid_overall_params(:,1), 'Normalization', 'probability', ...
                 'FaceColor', [0, 0, 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                 'DisplayName', 'Component-Compressed');
    end
% 新增：WTLS直方图
if ~isempty(valid_wtls_params)
    histogram(valid_wtls_params(:,1), 'Normalization', 'probability', ...
             'FaceColor', [0.9, 0.7, 0.9], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
             'DisplayName', 'WTLS (Data-snooping)');
end
    
    % 真实值标记
    xline(true_params(1), 'k--', 'LineWidth', 2, 'DisplayName', 'True Value');
    
    xlabel('Slope Estimate', 'FontSize', 15, 'FontName', 'Times New Roman');
    ylabel('Probability', 'FontSize', 15, 'FontName', 'Times New Roman');
    title('Slope Parameter Distribution (Histogram)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    legend('Location', 'best', 'FontSize', 15, 'FontName', 'Times New Roman');
    grid on;
    hold off;
    
    % 截距直方图对比
    subplot(2, 2, 4);
    hold on;
    
    if ~isempty(valid_dir_params)
        histogram(valid_dir_params(:,2), 'Normalization', 'probability', ...
                 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                 'DisplayName', 'Full-Component');
    end
    
    if ~isempty(valid_overall_params)
        histogram(valid_overall_params(:,2), 'Normalization', 'probability', ...
                 'FaceColor', [0, 0, 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                 'DisplayName', 'Component-Compressed');
    end
% 新增：WTLS直方图
if ~isempty(valid_wtls_params)
    histogram(valid_wtls_params(:,2), 'Normalization', 'probability', ...
             'FaceColor', [0.9, 0.7, 0.9], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
             'DisplayName', 'WTLS (Data-snooping)');
end
    
    % 真实值标记
    xline(true_params(2), 'k--', 'LineWidth', 2, 'DisplayName', 'True Value');
    
    xlabel('Intercept Estimate', 'FontSize', 15, 'FontName', 'Times New Roman');
    ylabel('Probability', 'FontSize', 15, 'FontName', 'Times New Roman');
    title('Intercept Parameter Distribution (Histogram)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    legend('Location', 'best', 'FontSize', 15, 'FontName', 'Times New Roman');
    grid on;
    hold off;
    
    % 保存图片（分辨率300 DPI）
    print(gcf, '随机粗差参数分布对比.png', '-dpng', '-r300');
    fprintf('参数分布图已保存为: 随机粗差参数分布对比.png (300 DPI)\n');
end

%% 结论
fprintf('\n========== 诊断结论 ==========\n');
if dir_success_rate < overall_success_rate
    fprintf('🔍 分方向方法在随机粗差下更容易收敛失败\n');
    fprintf('   可能原因：分方向残差分解在随机粗差情况下数值不稳定\n');
elseif overall_success_rate < dir_success_rate
    fprintf('🔍 总体方法在随机粗差下更容易收敛失败\n');
    fprintf('   可能原因：总体残差对随机粗差敏感，权重更新机制过于激进\n');
else
    fprintf('🔍 两种方法的收敛性能相近\n');
    fprintf('   说明随机粗差对两种方法的影响相似\n');
end

if dir_success_rate < 80 || overall_success_rate < 80
    fprintf('\n⚠️  建议：随机粗差设计存在数值稳定性问题\n');
    fprintf('   考虑调整抗差函数参数或降权策略\n');
end

%% ====== 绘制参数分布对比图 ======
% 使用已有的收敛性诊断实验结果
fprintf('\n========== 绘制参数分布对比图 ==========\n');

    % 提取成功收敛的结果
    successful_directional = convergence_results.directional_success == 1;
    successful_overall = convergence_results.overall_success == 1;
    successful_standard = convergence_results.standard_success == 1;
    successful_wtls = convergence_results.wtls_success == 1;
    successful_rtls = convergence_results.rtls_success == 1;

    if sum(successful_directional) > 0 && sum(successful_overall) > 0 && sum(successful_standard) > 0
        % 获取成功收敛的参数结果
        results_sep = convergence_results.directional_params(successful_directional, :);
        results_overall = convergence_results.overall_params(successful_overall, :);
        results_standard = convergence_results.standard_params(successful_standard, :);
        if sum(successful_wtls) > 0
            results_wtls = convergence_results.wtls_params(successful_wtls, :);
        else
            results_wtls = [];
        end
        if sum(successful_rtls) > 0
            results_rtls = convergence_results.rtls_params(successful_rtls, :);
        else
            results_rtls = [];
        end
        
        % 计算统计结果
        mean_sep = mean(results_sep);
        mean_overall = mean(results_overall);
        mean_standard = mean(results_standard);
        if ~isempty(results_wtls)
            mean_wtls = mean(results_wtls);
        else
            mean_wtls = [NaN, NaN];
        end
        if ~isempty(results_rtls)
            mean_rtls = mean(results_rtls);
        else
            mean_rtls = [NaN, NaN];
        end
        
        fprintf('使用已有实验结果：\n');
        fprintf('- 分方向方法成功样本: %d/%d (成功率: %.1f%%)\n', sum(successful_directional), max_experiments, sum(successful_directional)/max_experiments*100);
        fprintf('- 总体方法成功样本: %d/%d (成功率: %.1f%%)\n', sum(successful_overall), max_experiments, sum(successful_overall)/max_experiments*100);
        fprintf('- 标准TLS方法成功样本: %d/%d (成功率: %.1f%%)\n', sum(successful_standard), max_experiments, sum(successful_standard)/max_experiments*100);
        if ~isempty(results_wtls)
            fprintf('- WTLS数据探测法成功样本: %d/%d (成功率: %.1f%%)\n', sum(successful_wtls), max_experiments, sum(successful_wtls)/max_experiments*100);
        end
        if ~isempty(results_rtls)
            fprintf('- RTLS方法成功样本: %d/%d (成功率: %.1f%%)\n', sum(successful_rtls), max_experiments, sum(successful_rtls)/max_experiments*100);
        end
        fprintf('- 实际用于统计的样本数: 分方向=%d, 总体=%d, 标准=%d', size(results_sep,1), size(results_overall,1), size(results_standard,1));
        if ~isempty(results_wtls)
            fprintf(', WTLS=%d', size(results_wtls,1));
        end
        if ~isempty(results_rtls)
            fprintf(', RTLS=%d', size(results_rtls,1));
        end
        fprintf('\n');
    
    %% ====== 偏差和方差统计 ======
    fprintf('\n========== 偏差和方差统计 ==========\n');
    
    % 真实参数
    true_slope = true_params(1);    % 2
    true_intercept = true_params(2); % -3
    
    % 偏差计算
    bias_slope_sep = mean_sep(1) - true_slope;
    bias_slope_overall = mean_overall(1) - true_slope;
    bias_slope_standard = mean_standard(1) - true_slope;
    if ~isempty(results_wtls)
        bias_slope_wtls = mean_wtls(1) - true_slope;
    else
        bias_slope_wtls = NaN;
    end
    bias_intercept_sep = mean_sep(2) - true_intercept;
    bias_intercept_overall = mean_overall(2) - true_intercept;
    bias_intercept_standard = mean_standard(2) - true_intercept;
    if ~isempty(results_wtls)
        bias_intercept_wtls = mean_wtls(2) - true_intercept;
    else
        bias_intercept_wtls = NaN;
    end
    
    % 方差计算
    var_slope_sep = var(results_sep(:,1));
    var_slope_overall = var(results_overall(:,1));
    var_slope_standard = var(results_standard(:,1));
    if ~isempty(results_wtls)
        var_slope_wtls = var(results_wtls(:,1));
    else
        var_slope_wtls = NaN;
    end
    var_intercept_sep = var(results_sep(:,2));
    var_intercept_overall = var(results_overall(:,2));
    var_intercept_standard = var(results_standard(:,2));
    if ~isempty(results_wtls)
        var_intercept_wtls = var(results_wtls(:,2));
    else
        var_intercept_wtls = NaN;
    end
    
    % 标准差计算
    std_slope_sep = std(results_sep(:,1));
    std_slope_overall = std(results_overall(:,1));
    std_slope_standard = std(results_standard(:,1));
    if ~isempty(results_wtls)
        std_slope_wtls = std(results_wtls(:,1));
    else
        std_slope_wtls = NaN;
    end
    std_intercept_sep = std(results_sep(:,2));
    std_intercept_overall = std(results_overall(:,2));
    std_intercept_standard = std(results_standard(:,2));
    if ~isempty(results_wtls)
        std_intercept_wtls = std(results_wtls(:,2));
    else
        std_intercept_wtls = NaN;
    end
    
    % MSE计算
    mse_slope_sep = bias_slope_sep^2 + var_slope_sep;
    mse_slope_overall = bias_slope_overall^2 + var_slope_overall;
    mse_slope_standard = bias_slope_standard^2 + var_slope_standard;
    if ~isempty(results_wtls)
        mse_slope_wtls = bias_slope_wtls^2 + var_slope_wtls;
    else
        mse_slope_wtls = NaN;
    end
    mse_intercept_sep = bias_intercept_sep^2 + var_intercept_sep;
    mse_intercept_overall = bias_intercept_overall^2 + var_intercept_overall;
    mse_intercept_standard = bias_intercept_standard^2 + var_intercept_standard;
    if ~isempty(results_wtls)
        mse_intercept_wtls = bias_intercept_wtls^2 + var_intercept_wtls;
    else
        mse_intercept_wtls = NaN;
    end
    
    % 相对偏差和相对方差
    rel_bias_slope_sep = abs(bias_slope_sep) / abs(true_slope) * 100;
    rel_bias_slope_overall = abs(bias_slope_overall) / abs(true_slope) * 100;
    rel_bias_slope_standard = abs(bias_slope_standard) / abs(true_slope) * 100;
    if ~isempty(results_wtls)
        rel_bias_slope_wtls = abs(bias_slope_wtls) / abs(true_slope) * 100;
    else
        rel_bias_slope_wtls = NaN;
    end
    rel_bias_intercept_sep = abs(bias_intercept_sep) / abs(true_intercept) * 100;
    rel_bias_intercept_overall = abs(bias_intercept_overall) / abs(true_intercept) * 100;
    rel_bias_intercept_standard = abs(bias_intercept_standard) / abs(true_intercept) * 100;
    if ~isempty(results_wtls)
        rel_bias_intercept_wtls = abs(bias_intercept_wtls) / abs(true_intercept) * 100;
    else
        rel_bias_intercept_wtls = NaN;
    end
    
    rel_var_slope_sep = var_slope_sep / (true_slope^2) * 100;
    rel_var_slope_overall = var_slope_overall / (true_slope^2) * 100;
    rel_var_slope_standard = var_slope_standard / (true_slope^2) * 100;
    if ~isempty(results_wtls)
        rel_var_slope_wtls = var_slope_wtls / (true_slope^2) * 100;
    else
        rel_var_slope_wtls = NaN;
    end
    rel_var_intercept_sep = var_intercept_sep / (true_intercept^2) * 100;
    rel_var_intercept_overall = var_intercept_overall / (true_intercept^2) * 100;
    rel_var_intercept_standard = var_intercept_standard / (true_intercept^2) * 100;
    if ~isempty(results_wtls)
        rel_var_intercept_wtls = var_intercept_wtls / (true_intercept^2) * 100;
    else
        rel_var_intercept_wtls = NaN;
    end
    
    % 输出结果 - 使用科学计数法显示小数值
    fprintf('\n【斜率参数统计】\n');
    fprintf('真实值: %.6f\n', true_slope);
    fprintf('样本数: 分方向=%d, 总体=%d, 标准=%d', size(results_sep,1), size(results_overall,1), size(results_standard,1));
    if ~isempty(results_wtls)
        fprintf(', WTLS=%d\n', size(results_wtls,1));
    else
        fprintf('\n');
    end
    fprintf('分方向方法: 均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
        mean_sep(1), bias_slope_sep, var_slope_sep, std_slope_sep, mse_slope_sep);
    fprintf('总体方法:   均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
        mean_overall(1), bias_slope_overall, var_slope_overall, std_slope_overall, mse_slope_overall);
    fprintf('标准TLS法:  均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
        mean_standard(1), bias_slope_standard, var_slope_standard, std_slope_standard, mse_slope_standard);
    if ~isempty(results_wtls)
        fprintf('WTLS数据探测: 均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
            mean_wtls(1), bias_slope_wtls, var_slope_wtls, std_slope_wtls, mse_slope_wtls);
    end
    fprintf('相对偏差: 分方向=%.6f%%, 总体=%.6f%%, 标准=%.6f%%', rel_bias_slope_sep, rel_bias_slope_overall, rel_bias_slope_standard);
    if ~isempty(results_wtls)
        fprintf(', WTLS=%.6f%%\n', rel_bias_slope_wtls);
    else
        fprintf('\n');
    end
    fprintf('相对方差: 分方向=%.6f%%, 总体=%.6f%%, 标准=%.6f%%', rel_var_slope_sep, rel_var_slope_overall, rel_var_slope_standard);
    if ~isempty(results_wtls)
        fprintf(', WTLS=%.6f%%\n', rel_var_slope_wtls);
    else
        fprintf('\n');
    end
    
    fprintf('\n【截距参数统计】\n');
    fprintf('真实值: %.6f\n', true_intercept);
    fprintf('样本数: 分方向=%d, 总体=%d, 标准=%d', size(results_sep,1), size(results_overall,1), size(results_standard,1));
    if ~isempty(results_wtls)
        fprintf(', WTLS=%d\n', size(results_wtls,1));
    else
        fprintf('\n');
    end
    fprintf('分方向方法: 均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
        mean_sep(2), bias_intercept_sep, var_intercept_sep, std_intercept_sep, mse_intercept_sep);
    fprintf('总体方法:   均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
        mean_overall(2), bias_intercept_overall, var_intercept_overall, std_intercept_overall, mse_intercept_overall);
    fprintf('标准TLS法:  均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
        mean_standard(2), bias_intercept_standard, var_intercept_standard, std_intercept_standard, mse_intercept_standard);
    if ~isempty(results_wtls)
        fprintf('WTLS数据探测: 均值=%.6f, 偏差=%.2e, 方差=%.2e, 标准差=%.6f, MSE=%.2e\n', ...
            mean_wtls(2), bias_intercept_wtls, var_intercept_wtls, std_intercept_wtls, mse_intercept_wtls);
    end
    fprintf('相对偏差: 分方向=%.6f%%, 总体=%.6f%%, 标准=%.6f%%', rel_bias_intercept_sep, rel_bias_intercept_overall, rel_bias_intercept_standard);
    if ~isempty(results_wtls)
        fprintf(', WTLS=%.6f%%\n', rel_bias_intercept_wtls);
    else
        fprintf('\n');
    end
    fprintf('相对方差: 分方向=%.6f%%, 总体=%.6f%%, 标准=%.6f%%', rel_var_intercept_sep, rel_var_intercept_overall, rel_var_intercept_standard);
    if ~isempty(results_wtls)
        fprintf(', WTLS=%.6f%%\n', rel_var_intercept_wtls);
    else
        fprintf('\n');
    end
    
    % 添加详细的计算过程说明
    fprintf('\n【计算过程说明】\n');
    fprintf('偏差 = 均值 - 真实值\n');
    fprintf('方差 = E[(X - E[X])²]\n');
    fprintf('标准差 = sqrt(方差)\n');
    fprintf('MSE = 偏差² + 方差\n');
    fprintf('相对偏差 = |偏差| / |真实值| × 100%%\n');
    fprintf('相对方差 = 方差 / 真实值² × 100%%\n');
    
    % 显示实际的数据范围
    fprintf('\n【数据范围】\n');
    fprintf('斜率 - 分方向方法: [%.6f, %.6f], 范围=%.2e\n', min(results_sep(:,1)), max(results_sep(:,1)), max(results_sep(:,1))-min(results_sep(:,1)));
    fprintf('斜率 - 总体方法:   [%.6f, %.6f], 范围=%.2e\n', min(results_overall(:,1)), max(results_overall(:,1)), max(results_overall(:,1))-min(results_overall(:,1)));
    fprintf('斜率 - 标准TLS法:  [%.6f, %.6f], 范围=%.2e\n', min(results_standard(:,1)), max(results_standard(:,1)), max(results_standard(:,1))-min(results_standard(:,1)));
    if ~isempty(results_wtls)
        fprintf('斜率 - WTLS数据探测: [%.6f, %.6f], 范围=%.2e\n', min(results_wtls(:,1)), max(results_wtls(:,1)), max(results_wtls(:,1))-min(results_wtls(:,1)));
    end
    fprintf('截距 - 分方向方法: [%.6f, %.6f], 范围=%.2e\n', min(results_sep(:,2)), max(results_sep(:,2)), max(results_sep(:,2))-min(results_sep(:,2)));
    fprintf('截距 - 总体方法:   [%.6f, %.6f], 范围=%.2e\n', min(results_overall(:,2)), max(results_overall(:,2)), max(results_overall(:,2))-min(results_overall(:,2)));
    fprintf('截距 - 标准TLS法:  [%.6f, %.6f], 范围=%.2e\n', min(results_standard(:,2)), max(results_standard(:,2)), max(results_standard(:,2))-min(results_standard(:,2)));
    if ~isempty(results_wtls)
        fprintf('截距 - WTLS数据探测: [%.6f, %.6f], 范围=%.2e\n', min(results_wtls(:,2)), max(results_wtls(:,2)), max(results_wtls(:,2))-min(results_wtls(:,2)));
    end
    
    % 性能比较
    fprintf('\n【性能比较 (基于MSE)】\n');
    % 斜率参数比较
    mse_slope_list = [mse_slope_sep, mse_slope_overall, mse_slope_standard];
    slope_methods = {'分方向方法', '总体方法', '标准TLS法'};
    if ~isempty(results_wtls)
        mse_slope_list = [mse_slope_list, mse_slope_wtls];
        slope_methods = [slope_methods, {'WTLS数据探测'}];
    end
    [~, best_slope_idx] = min(mse_slope_list);
    fprintf('斜率参数: %s最优 (MSE: %.2e)\n', slope_methods{best_slope_idx}, min(mse_slope_list));
    
    % 截距参数比较
    mse_intercept_list = [mse_intercept_sep, mse_intercept_overall, mse_intercept_standard];
    intercept_methods = {'分方向方法', '总体方法', '标准TLS法'};
    if ~isempty(results_wtls)
        mse_intercept_list = [mse_intercept_list, mse_intercept_wtls];
        intercept_methods = [intercept_methods, {'WTLS数据探测'}];
    end
    [~, best_intercept_idx] = min(mse_intercept_list);
    fprintf('截距参数: %s最优 (MSE: %.2e)\n', intercept_methods{best_intercept_idx}, min(mse_intercept_list));

% ====== 抗差方法对比图（包含WTLS数据探测法和RTLS方法）======
figure('Position', [100, 100, 1200, 800]);

% 提取WTLS和RTLS成功结果（RTLS用于计算bin边界，但不绘制）
successful_wtls = convergence_results.wtls_success == 1;
valid_wtls_params = convergence_results.wtls_params(successful_wtls, :);
successful_rtls = convergence_results.rtls_success == 1;
valid_rtls_params = convergence_results.rtls_params(successful_rtls, :);

% 斜率参数分布对比 - 抗差方法（包含WTLS和RTLS）
subplot(2,1,1);
hold on;

% 计算斜率的统一bin边界 - 包含所有四种抗差方法（包含RTLS用于确定范围，但不绘制）
slope_min = min([min(results_sep(:,1)), min(results_overall(:,1))]);
slope_max = max([max(results_sep(:,1)), max(results_overall(:,1))]);
if ~isempty(valid_wtls_params)
    slope_min = min([slope_min, min(valid_wtls_params(:,1))]);
    slope_max = max([slope_max, max(valid_wtls_params(:,1))]);
end
if ~isempty(valid_rtls_params)
    slope_min = min([slope_min, min(valid_rtls_params(:,1))]);
    slope_max = max([slope_max, max(valid_rtls_params(:,1))]);
end
slope_bins = linspace(slope_min, slope_max, 41);  % 40个区间，41个边界点

% 绘制抗差方法的直方图
% Newton Full-Com (WTLS): 红色, Newton Com-Com (Component-Compressed): 蓝色, Mahboub Full-Com (Full-Component): 绿色
if ~isempty(valid_wtls_params)
    h1 = histogram(valid_wtls_params(:,1), slope_bins, 'Normalization', 'probability', 'FaceColor', [1, 0, 0], 'EdgeColor', [0.8, 0, 0], 'LineWidth', 0.5, 'FaceAlpha', 0.4, 'DisplayName', 'Newton Full-Com (Histogram)');
end
h2 = histogram(results_overall(:,1), slope_bins, 'Normalization', 'probability', 'FaceColor', [0, 0, 1], 'EdgeColor', [0, 0, 0.8], 'LineWidth', 0.5, 'FaceAlpha', 0.4, 'DisplayName', 'Newton Com-Com (Histogram)');
h3 = histogram(results_sep(:,1), slope_bins, 'Normalization', 'probability', 'FaceColor', [0, 0.8, 0], 'EdgeColor', [0, 0.6, 0], 'LineWidth', 0.5, 'FaceAlpha', 0.4, 'DisplayName', 'Mahboub Full-Com (Histogram)');

% 绘制核密度曲线
slope_range = linspace(slope_min, slope_max, 200);
[f_sep_slope, xi_sep_slope] = ksdensity(results_sep(:,1), slope_range);
[f_overall_slope, xi_overall_slope] = ksdensity(results_overall(:,1), slope_range);
if ~isempty(valid_wtls_params)
    [f_wtls_slope, xi_wtls_slope] = ksdensity(valid_wtls_params(:,1), slope_range);
end

% 缩放密度曲线以匹配直方图的高度
max_sep_hist = max(h3.Values);
max_overall_hist = max(h2.Values);
f_sep_scaled = f_sep_slope * max_sep_hist / max(f_sep_slope);
f_overall_scaled = f_overall_slope * max_overall_hist / max(f_overall_slope);
if ~isempty(valid_wtls_params)
    max_wtls_hist = max(h1.Values);
    f_wtls_scaled = f_wtls_slope * max_wtls_hist / max(f_wtls_slope);
end

% 绘制密度曲线
if ~isempty(valid_wtls_params)
    plot(xi_wtls_slope, f_wtls_scaled, 'Color', [1, 0, 0], 'LineWidth', 2, 'DisplayName', 'Newton Full-Com (KDE)');  % 红色
end
plot(xi_overall_slope, f_overall_scaled, 'Color', [0, 0, 1], 'LineWidth', 2, 'DisplayName', 'Newton Com-Com (KDE)');  % 蓝色
plot(xi_sep_slope, f_sep_scaled, 'Color', [0, 0.8, 0], 'LineWidth', 2, 'DisplayName', 'Mahboub Full-Com (KDE)');  % 绿色

% 添加参考线
xline(true_params(1), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');

title('Slope Parameter Distribution', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Slope Value', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Frequency/Density', 'FontSize', 15, 'FontName', 'Times New Roman');
legend('Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;

% 截距参数分布对比 - 抗差方法（包含WTLS）
subplot(2,1,2);
hold on;

% 计算截距的统一bin边界 - 包含所有四种抗差方法（包含RTLS用于确定范围，但不绘制）
intercept_min = min([min(results_sep(:,2)), min(results_overall(:,2))]);
intercept_max = max([max(results_sep(:,2)), max(results_overall(:,2))]);
if ~isempty(valid_wtls_params)
    intercept_min = min([intercept_min, min(valid_wtls_params(:,2))]);
    intercept_max = max([intercept_max, max(valid_wtls_params(:,2))]);
end
if ~isempty(valid_rtls_params)
    intercept_min = min([intercept_min, min(valid_rtls_params(:,2))]);
    intercept_max = max([intercept_max, max(valid_rtls_params(:,2))]);
end
intercept_bins = linspace(intercept_min, intercept_max, 41);  % 40个区间，41个边界点

% 绘制抗差方法的直方图
% Newton Full-Com (WTLS): 红色, Newton Com-Com (Component-Compressed): 蓝色, Mahboub Full-Com (Full-Component): 绿色
if ~isempty(valid_wtls_params)
    h5 = histogram(valid_wtls_params(:,2), intercept_bins, 'Normalization', 'probability', 'FaceColor', [1, 0, 0], 'EdgeColor', [0.8, 0, 0], 'LineWidth', 0.5, 'FaceAlpha', 0.4, 'DisplayName', 'Newton Full-Com (Histogram)');
end
h6 = histogram(results_overall(:,2), intercept_bins, 'Normalization', 'probability', 'FaceColor', [0, 0, 1], 'EdgeColor', [0, 0, 0.8], 'LineWidth', 0.5, 'FaceAlpha', 0.4, 'DisplayName', 'Newton Com-Com (Histogram)');
h7 = histogram(results_sep(:,2), intercept_bins, 'Normalization', 'probability', 'FaceColor', [0, 0.8, 0], 'EdgeColor', [0, 0.6, 0], 'LineWidth', 0.5, 'FaceAlpha', 0.4, 'DisplayName', 'Mahboub Full-Com (Histogram)');

% 绘制核密度曲线
intercept_range = linspace(intercept_min, intercept_max, 200);
[f_sep_intercept, xi_sep_intercept] = ksdensity(results_sep(:,2), intercept_range);
[f_overall_intercept, xi_overall_intercept] = ksdensity(results_overall(:,2), intercept_range);
if ~isempty(valid_wtls_params)
    [f_wtls_intercept, xi_wtls_intercept] = ksdensity(valid_wtls_params(:,2), intercept_range);
end

% 缩放密度曲线以匹配直方图的高度
max_sep_hist_intercept = max(h7.Values);
max_overall_hist_intercept = max(h6.Values);
f_sep_intercept_scaled = f_sep_intercept * max_sep_hist_intercept / max(f_sep_intercept);
f_overall_intercept_scaled = f_overall_intercept * max_overall_hist_intercept / max(f_overall_intercept);
if ~isempty(valid_wtls_params)
    max_wtls_hist_intercept = max(h5.Values);
    f_wtls_intercept_scaled = f_wtls_intercept * max_wtls_hist_intercept / max(f_wtls_intercept);
end

% 绘制密度曲线
if ~isempty(valid_wtls_params)
    plot(xi_wtls_intercept, f_wtls_intercept_scaled, 'Color', [1, 0, 0], 'LineWidth', 2, 'DisplayName', 'Newton Full-Com (KDE)');  % 红色
end
plot(xi_overall_intercept, f_overall_intercept_scaled, 'Color', [0, 0, 1], 'LineWidth', 2, 'DisplayName', 'Newton Com-Com (KDE)');  % 蓝色
plot(xi_sep_intercept, f_sep_intercept_scaled, 'Color', [0, 0.8, 0], 'LineWidth', 2, 'DisplayName', 'Mahboub Full-Com (KDE)');  % 绿色

% 添加参考线
xline(true_params(2), 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');

title('Intercept Parameter Distribution', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Intercept Value', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Frequency/Density', 'FontSize', 15, 'FontName', 'Times New Roman');
legend('Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;

% 保存抗差方法对比图（分辨率300 DPI）
print(gcf, 'Robust_Methods_Comparison.png', '-dpng', '-r300');
fprintf('抗差方法对比图已保存为: Robust_Methods_Comparison.png (300 DPI)\n');

% ====== 标准TLS方法单独展示 ======
figure('Position', [1300, 100, 800, 800]);

% 斜率参数分布 - 标准TLS方法
subplot(2,1,1);
hold on;

% 计算斜率的bin边界 - 只包含标准方法
slope_min_std = min(results_standard(:,1));
slope_max_std = max(results_standard(:,1));
slope_bins_std = linspace(slope_min_std, slope_max_std, 41);

% 绘制标准方法的直方图
h3 = histogram(results_standard(:,1), slope_bins_std, 'Normalization', 'probability', 'FaceColor', [0.7, 0.7, 0.9], 'EdgeColor', [0.5, 0.5, 0.8], 'LineWidth', 0.5, 'FaceAlpha', 0.4);

% 绘制标准方法的核密度曲线
slope_range_std = linspace(slope_min_std, slope_max_std, 200);
[f_standard_slope, xi_standard_slope] = ksdensity(results_standard(:,1), slope_range_std);
f_standard_scaled = f_standard_slope * max(h3.Values) / max(f_standard_slope);
plot(xi_standard_slope, f_standard_scaled, 'Color', [0.0, 0.0, 0.8], 'LineWidth', 2);

% 添加参考线
xline(true_params(1), 'k--', 'LineWidth', 3);

title('Slope Parameter Distribution - Ordinary TLS', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Slope Value', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Probability Density', 'FontSize', 15, 'FontName', 'Times New Roman');
legend('Ordinary TLS (Histogram)', 'Ordinary TLS (KDE)', 'True Value', 'Location', 'best', 'FontSize', 15, 'FontName', 'Times New Roman');
grid on;

% 截距参数分布 - 标准TLS方法
subplot(2,1,2);
hold on;

% 计算截距的bin边界 - 只包含标准方法
intercept_min_std = min(results_standard(:,2));
intercept_max_std = max(results_standard(:,2));
intercept_bins_std = linspace(intercept_min_std, intercept_max_std, 41);

% 绘制标准方法的直方图
h6 = histogram(results_standard(:,2), intercept_bins_std, 'Normalization', 'probability', 'FaceColor', [0.7, 0.7, 0.9], 'EdgeColor', [0.5, 0.5, 0.8], 'LineWidth', 0.5, 'FaceAlpha', 0.4);

% 绘制标准方法的核密度曲线
intercept_range_std = linspace(intercept_min_std, intercept_max_std, 200);
[f_standard_intercept, xi_standard_intercept] = ksdensity(results_standard(:,2), intercept_range_std);
f_standard_intercept_scaled = f_standard_intercept * max(h6.Values) / max(f_standard_intercept);
plot(xi_standard_intercept, f_standard_intercept_scaled, 'Color', [0.0, 0.0, 0.8], 'LineWidth', 2);

% 添加参考线
xline(true_params(2), 'k--', 'LineWidth', 3);
xline(mean_standard(2), 'b--', 'LineWidth', 2.5);

title('Intercept Parameter Distribution - Ordinary TLS', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Intercept Value', 'FontSize', 15, 'FontName', 'Times New Roman');
ylabel('Probability Density', 'FontSize', 15, 'FontName', 'Times New Roman');
legend('Ordinary TLS (Histogram)', 'Ordinary TLS (KDE)', 'True Value', 'Ordinary TLS Mean', 'Location', 'best', 'FontSize', 15, 'FontName', 'Times New Roman');
grid on;

% 保存标准TLS方法图（分辨率300 DPI）
print(gcf, 'Ordinary_TLS_Distribution.png', '-dpng', '-r300');
fprintf('标准TLS方法图已保存为: Ordinary_TLS_Distribution.png (300 DPI)\n');

    fprintf('\n========== 绘图完成 ==========\n');
else
    fprintf('\n========== 绘图失败 ==========\n');
    fprintf('没有足够的成功收敛样本进行绘图\n');
end

%% 算法函数定义（与原实验完全一致）

function [X_hat, info] = wtls_data_snooping_estimate(A, L)
% 最小实现：WTLS + 数据探测（Amiri-Simkooei & Jazaeri, 2013）
% 仅用于对比绘图：输入A=[x,1], L=y

    y = L;
    x_obs = A(:,1);
    m = size(A, 1);
    n = size(A, 2); % 2

    % 初始协方差（单位权）
    Q_y = eye(m);
    Q_A = zeros(m*n, m*n);
    % vec([x,1]) = [x; 1] → 前m是x列，有误差；后m是常数列，无误差
    Q_A(1:m, 1:m) = eye(m);
    Q_A(m+1:end, m+1:end) = 0;

    % 外层数据探测
    % ============================================================
    % 数据探测法的停止条件：
    % 1. 主要停止条件：假设检验找不到粗差（wmax <= Fcrit）→ 自然停止
    % 2. 安全上限：max_snoop 防止以下异常情况：
    %    a) 数值误差导致w统计量在临界值附近波动，形成循环
    %    b) 数据质量极差，可能误删正常点，导致算法无法收敛
    %    c) 计算效率：防止计算时间过长
    % 3. 最小观测数限制：m_current > n+2（至少需要n+2个观测才能进行假设检验）
    % ============================================================
    % 理论上，当假设检验找不到粗差时，算法会通过 else break 自然停止
    % max_snoop 只是一个安全上限，应该设置得足够大（如数据量的50%）
    % 但也不能太大，避免误删过多正常点
    max_snoop = 15;%min(ceil(m * 0.5), 50);  % 最多探测数据量的50%，但不超过50个
    snoop_iters = 0;
    removed = 0;
    converged = true;

    point_no = (1:m)';

    while snoop_iters < max_snoop
        snoop_iters = snoop_iters + 1;
        m_current = size(y, 1);
        
        % 检查最小观测数要求
        if m_current <= n + 2
            % 观测数不足，无法进行假设检验，停止算法
            break;
        end
        
        A_cur = [x_obs, ones(m_current,1)];

        % 初值：OLS
        x_hat = (A_cur' * A_cur) \ (A_cur' * y);

        % 内层WTLS
        eps_tol = 1e-6;
        max_inner = 100;
        inner_done = 0;

        for it = 1:max_inner
            inner_done = it;
            e_hat = y - A_cur * x_hat;

            x_kron_T = kron(x_hat', eye(m_current)); % m × mn
            x_kron   = kron(x_hat , eye(m_current)); % mn × m
            Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
            Q_y_tilde_inv = inv(Q_y_tilde);

            vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
            E_A_hat = [vec_E_A(1:m_current), vec_E_A(m_current+1:end)];

            A_tilde = A_cur - E_A_hat;
            y_tilde = y - E_A_hat * x_hat;

            x_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * y_tilde);
            if norm(x_new - x_hat) < eps_tol, break; end
            x_hat = x_new;
        end

        % w检验（假设检验）
        e_total = y - A_cur * x_hat;
        Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);
        Q_e_norm = Q_y_tilde - A_tilde * Q_x * A_tilde';

        sigma0_sq = (e_total' * Q_y_tilde_inv * e_total) / (m_current - n);
        sigma0 = sqrt(sigma0_sq);

        w = zeros(m_current,1);
        for i = 1:m_current
            e_i = zeros(m_current,1); e_i(i) = 1;
            num = e_i' * Q_y_tilde_inv * e_total;
            den = sigma0 * sqrt(e_i' * Q_y_tilde_inv * Q_e_norm * Q_y_tilde_inv * e_i);
            w(i) = num / den;
        end

        [wmax, idx] = max(abs(w));
        Fcrit = sqrt(finv(0.95, 1, m_current - n));
        
        % ============================================================
        % 主要停止条件：假设检验找不到粗差（wmax <= Fcrit）
        % 这是算法的自然停止条件，理论上应该通过这个条件停止
        % ============================================================
        if wmax > Fcrit && m_current > n+2
            % 找到粗差，删除并继续探测
            keep = true(m_current,1); keep(idx) = false;
            x_obs = x_obs(keep);
            y = y(keep);
            point_no = point_no(keep);
            Q_y = Q_y(keep, keep);
            keep_exp = [keep; keep];
            Q_A = Q_A(keep_exp, keep_exp);
            removed = removed + 1;
            continue;  % 继续下一次探测
        else
            % ============================================================
            % 自然停止：假设检验找不到粗差（wmax <= Fcrit）
            % 或者观测数不足（m_current <= n+2）
            % 这是算法的正常停止条件
            % ============================================================
            break;  % 自然停止，不再探测
        end
    end
    
    % 检查是否因为达到max_snoop而停止（可能是异常情况）
    if snoop_iters >= max_snoop && removed > 0
        % 达到上限但仍有粗差被删除，可能是数据中粗差较多
        % 这是正常情况，不需要警告
    end

    X_hat = x_hat;
    info = struct('converged', converged, ...
                  'inner_iterations', inner_done, ...
                  'snooping_iterations', snoop_iters, ...
                  'outliers_removed', removed, ...
                  'total_iterations', inner_done);
end



function X = TLS_newton_2(A, L, P, X0_init)
[m, n] = size(P);
PP = P(2:m, :);
for i = 1:n
    Pi{i} = diag(PP(:,i));
end

if nargin < 4 || isempty(X0_init)
    p0 = P(1,:);
    P0 = diag(p0);
    X0 = pinv(A' * P0 * A) * A' * P0 * L;
else
    X0 = X0_init;
end

cita = 1;
iter_count = 0;
while cita > 1e-10
    iter_count = iter_count + 1;
    v = L - A * X0;
    H1 = 0;
    w = zeros(n,1);
    for i = 1:n
        p_i = P(:,i);
        zero_pos = find(p_i == 0, 1);
        if isempty(zero_pos)
            w(i) = p_i(1) / (1 + p_i(1) * X0' * pinv(Pi{i}) * X0);
            e{i} = w(i) * v(i) * pinv(Pi{i}) * X0;
            E(i,:) = e{i}';
        else
            w(i) = 0;
            if zero_pos == 1
                v(i) = A(i,:) * X0 - L(i);
                e{i} = zeros(size(Pi{i},1),1);
            else
                k = zero_pos - 1;
                v(i) = 0;
                e_vector = zeros(size(Pi{i},1),1);
                e_vector(k) = (L(i) - A(i,:)*X0) / X0(k);
                e{i} = e_vector;
            end
            E(i,:) = e{i}';
        end
    end
    
    W = diag(w);
    H2 = -4 * A' * W * E;
    H3 = -A' * W * A;
    H4 = -4 * E' * W * E;
    
    for i = 1:n
        if w(i) > 0
            H1 = H1 + w(i)^2 * v(i)^2 * pinv(Pi{i});
        end
    end
    
    F = (A + E)' * W * v;
    H = H1 + H2 + H3 + H4;
    dX = pinv(H) * F;
    X = X0 - dX;
    cita = norm(dX);
    X0 = X;
end
end

function X = standard_TLS(A, L, P_initial)
% 标准TLS方法，不进行迭代权重调整
% 输入：
%   A - 设计矩阵
%   L - 观测向量
%   P_initial - 初始权重矩阵
% 输出：
%   X - 参数估计结果
% 直接使用TLS_newton_2函数求解，不进行权重迭代
X = TLS_newton_2(A, L, P_initial);
end

function [X, residuals, iter_info] = iterative_weight_optimization_with_timeout(A, L, P_initial, timeout_seconds)
% 带超时机制的分方向残差权重优化
% 输入：
%   A - 设计矩阵
%   L - 观测向量
%   P_initial - 初始权重矩阵
%   timeout_seconds - 超时时间（秒）
% 输出：
%   X - 参数估计结果
%   residuals - 残差信息
%   iter_info - 迭代信息

start_time = tic;
[m, n] = size(P_initial);
% 使用绝对误差，考虑参数理论精度
% 基于观测精度（标准差0.1）和样本量，理论精度约为 0.1/sqrt(n)
param_tol_abs = 1e-3;  % 绝对误差阈值，考虑参数理论精度
max_iterations = 100;  % 添加最大迭代次数限制

% 初始化
iter_info = struct();
iter_info.total_iterations = 0;
P = P_initial;
param_diff = 1;
iter_count = 0;

% 初始最小二乘解
X0 = TLS_newton_2(A, L, P);

while param_diff > param_tol_abs && iter_count < max_iterations

    iter_count = iter_count + 1;
    if iter_count > 1
        X_prev = X0;
    end

    % 计算分方向残差
    v = L - A * X0;
    e_y = zeros(n,1);
    e_x1 = zeros(n,1);
    
    for i = 1:n
        p_i = P(:,i);
        zero_pos = find(p_i(1:2) == 0, 1);

        if isempty(zero_pos)
            % 标准情况
            p_simplified = p_i(1:2);
            B_i_simplified = [1, X0(1)];
            Pi_inv_simplified = diag(1./p_simplified);

            BPiB_simplified = B_i_simplified * Pi_inv_simplified * B_i_simplified';
            ei_simplified = Pi_inv_simplified * B_i_simplified' * (1/BPiB_simplified) * v(i);

            e_y(i) = ei_simplified(1);
            e_x1(i) = ei_simplified(2);
        else
            % 零权处理
            if zero_pos == 1
                e_y(i) = v(i);
                e_x1(i) = 0;
            else
                e_y(i) = 0;
                e_x1(i) = v(i);
            end
        end
    end

    % 计算统一的单位权中误差
    r = n - size(A,2);

    % 计算每个观测的整体权重
    w = zeros(n,1);
    PP = P(2:m, :);
    for i = 1:n
        Pi_diag = diag(PP(:,i));
        p_i = P(:,i);
        if all(p_i ~= 0)
            w(i) = p_i(1) / (1 + p_i(1) * X0' * pinv(Pi_diag) * X0);
        else
            w(i) = 0;
        end
    end

    rho = v' * diag(w) * v;
    sigma0 = sqrt(rho / r);

    % 更新权重矩阵
    k0 = 1.5;
    k1 = 2.5;
    min_weight = 0.0;   

    for i = 1:n
        e_bar_y = abs(e_y(i)) / sigma0;
        if e_bar_y <= k0
            q_y = 1.0;
        elseif e_bar_y <= k1
            q_y = (k0/e_bar_y) * ((k1 - e_bar_y)/(k1 - k0))^2;
        else
            q_y = min_weight;
        end
        P(1,i) = P_initial(1,i) * q_y;

        e_bar_x1 = abs(e_x1(i)) / sigma0;
        if e_bar_x1 <= k0
            q_x1 = 1.0;
        elseif e_bar_x1 <= k1
            q_x1 = (k0/e_bar_x1) * ((k1 - e_bar_x1)/(k1 - k0))^2;
        else
            q_x1 = min_weight;
        end
        P(2,i) = P_initial(2,i) * q_x1;
    end

    % 更新参数估计
    X0 = TLS_newton_2(A, L, P, X0);

    % 使用绝对误差判定收敛
    if iter_count > 1
        param_diff = norm(X0 - X_prev);  % 绝对误差
    end
    if iter_count < 3
        param_diff = 1;
    end
    
    % 添加调试信息（每10次迭代输出一次）
    if mod(iter_count, 10) == 0
        fprintf('  分方向方法迭代 %d: 参数绝对变化 = %.6e\n', iter_count, param_diff);
    end
end

% 检查是否达到最大迭代次数
if iter_count >= max_iterations
    warning('分方向方法达到最大迭代次数 %d，可能未完全收敛', max_iterations);
end

% 计算最终结果
iter_info.total_iterations = iter_count;
X = X0;
residuals = struct('e_y', [], 'e_x1', []);
end

function [X, residuals, iter_info] = overall_residual_weight_optimization_with_timeout(A, L, P_initial, timeout_seconds)
% 带超时机制的总体残差权重优化
% 输入：
%   A - 设计矩阵
%   L - 观测向量
%   P_initial - 初始权重矩阵
%   timeout_seconds - 超时时间（秒）
% 输出：
%   X - 参数估计结果
%   residuals - 残差信息
%   iter_info - 迭代信息

start_time = tic;
[m, n] = size(P_initial);
% 使用绝对误差，考虑参数理论精度
% 基于观测精度（标准差0.1）和样本量，理论精度约为 0.1/sqrt(n)
param_tol_abs = 1e-2;  % 绝对误差阈值，考虑参数理论精度
max_iterations = 100;  % 添加最大迭代次数限制

% 初始化
iter_info = struct();
iter_info.total_iterations = 0;
P = P_initial;
param_diff = 1;
iter_count = 0;

% 初始最小二乘解
X0 = TLS_newton_2(A, L, P);

while param_diff > param_tol_abs && iter_count < max_iterations
    % 检查超时
    if toc(start_time) > timeout_seconds
        warning('总体方法超时(%.1f秒)，提前终止', timeout_seconds);
        break;
    end
    
    iter_count = iter_count + 1;
    if iter_count > 1
        X_prev = X0;
    end

    % 计算总体残差
    v = L - A * X0;
    
    % 检查数值稳定性
    if any(isnan(X0)) || any(isinf(X0))
        warning('总体方法参数出现NaN或Inf，算法异常终止');
        X = X_prev;  % 返回上一步的结果
        break;
    end
    
    % 检查参数发散
    if norm(X0) > 1e6
        warning('总体方法参数发散(norm>1e6)，算法异常终止');
        X = X_prev;  % 返回上一步的结果
        break;
    end

    % 计算总体残差的单位权中误差
    r = n - size(A,2);

    % 计算每个观测的整体权重
    w = zeros(n,1);
    PP = P(2:m, :);
    for i = 1:n
        Pi_diag = diag(PP(:,i));
        p_i = P(:,i);
        if all(p_i ~= 0)
            w(i) = p_i(1) / (1 + p_i(1) * X0' * pinv(Pi_diag) * X0);
        else
            w(i) = 0;
        end
    end

    rho = v' * diag(w) * v;
    sigma0 = sqrt(rho / r);
    
    % 检查单位权中误差的有效性
    if isnan(sigma0) || isinf(sigma0) || sigma0 < 1e-10
        warning('总体方法单位权中误差异常(%.2e)，算法终止', sigma0);
        X = X_prev;
        break;
    end

    % 基于总体残差更新所有方向的权重
    for i = 1:n
        e_bar = abs(v(i)) / sigma0;

        k0 = 1.5;
        k1 = 2.5;
        if e_bar <= k0
            q = 1.0;
        elseif e_bar <= k1
            q = (k0/e_bar) * ((k1 - e_bar)/(k1 - k0))^2;
        else
            q = 0;
        end

        P(1,i) = P_initial(1,i) * q;
        P(2,i) = P_initial(2,i) * q;
    end
    
    % 检查权重矩阵是否全部为零
    if sum(P(1,:)) == 0 || sum(P(2,:)) == 0
        warning('总体方法权重全部为零，算法终止');
        X = X_prev;
        break;
    end
    
    % 检查权重有效性（至少保留15%的观测）
    effective_obs = sum(P(1,:) > 0.01 * P_initial(1,1));
    min_required_obs = max(ceil(n * 0.15), 5);  % 至少15%或5个观测，取较大者
    
    if effective_obs < min_required_obs
        warning('总体方法有效观测过少(%d/%d)，低于最低要求(%d)，算法终止', ...
                effective_obs, n, min_required_obs);
        X = X_prev;
        break;
    end
    
    % 使用try-catch保护TLS求解
    try
        X0_new = TLS_newton_2(A, L, P, X0);
        
        % 检查TLS求解后的结果
        if any(isnan(X0_new)) || any(isinf(X0_new))
            warning('总体方法TLS求解返回NaN或Inf，算法终止');
            X = X_prev;
            break;
        end
        
        % 检查参数变化是否异常巨大（可能数值不稳定）
        if iter_count > 1 && norm(X0_new - X0) > 100 * norm(X0)
            warning('总体方法参数变化异常巨大，算法可能不稳定，终止');
            X = X_prev;
            break;
        end
        
        X0 = X0_new;
        
    catch ME
        warning('总体方法TLS求解异常: %s，算法终止', ME.message);
        X = X_prev;
        break;
    end

    % 使用绝对误差判定收敛
    if iter_count > 1
        param_diff = norm(X0 - X_prev);  % 绝对误差
    end
%     if iter_count < 3
%         param_diff = 1;
%     end
    % 添加调试信息（每10次迭代输出一次）
    if mod(iter_count, 10) == 0
        fprintf('  总体方法迭代 %d: 参数绝对变化 = %.6e\n', iter_count, param_diff);
    end
end

% 检查是否达到最大迭代次数
if iter_count >= max_iterations
    warning('总体方法达到最大迭代次数 %d，可能未完全收敛', max_iterations);
end

iter_info.total_iterations = iter_count;
X = X0;
v_final = L - A * X;
residuals = struct('v', v_final);
end

