%% 全分量抗差方法收敛性诊断实验
% 目的：测试随机粗差设计下，全分量抗差方法的收敛性
% 诊断方法在随机粗差情况下的表现

clear; clc; close all;
rng(42);
fprintf('========== 全分量抗差方法收敛性诊断 ==========\n');
fprintf('测试随机粗差下全分量抗差方法的收敛性能\n\n');

%% 实验参数设置
n = 20;  % 样本量
true_params = [2; -3];  % 真实参数 y = 2x - 3
outlier_ratio_3sigma = 0.10;  % 10%的粗差为3倍标准差
outlier_ratio_5sigma = 0.05;  % 5%的粗差为5倍标准差
outlier_ratio = outlier_ratio_3sigma + outlier_ratio_5sigma;  % 总粗差比例15%
outlier_magnitude_3sigma = 3;  % 3倍粗差
outlier_magnitude_5sigma = 5;  % 5倍粗差
max_experiments = 1000;  % 测试次数

% 存储收敛信息
convergence_results = struct();
convergence_results.robust_success = [];
convergence_results.robust_iters = [];
convergence_results.robust_params = [];
convergence_results.robust_errors = {};
convergence_results.robust_times = [];

% 存储全分量数据探测方法的结果
convergence_results.full_detection_success = [];
convergence_results.full_detection_iters = [];
convergence_results.full_detection_params = [];
convergence_results.full_detection_errors = {};
convergence_results.full_detection_times = [];

fprintf('实验设置：\n');
fprintf('- 样本量: %d\n', n);
fprintf('- 真实参数: y = %.1fx + %.1f\n', true_params(1), true_params(2));
fprintf('- 粗差比例: %.0f%% (%.0f%%为3倍标准差, %.0f%%为5倍标准差)\n', ...
    outlier_ratio*100, outlier_ratio_3sigma*100, outlier_ratio_5sigma*100);
fprintf('- 测试次数: %d\n', max_experiments);
fprintf('- 方法: 全分量抗差法（迭代权重优化） + 全分量数据探测法（迭代检测）\n\n');

%% 开始诊断实验
for exp_idx = 1:max_experiments
    fprintf('========== 实验 %d/%d ==========\n', exp_idx, max_experiments);
    
    %% 数据生成（与原实验一致）
    % 步骤1：生成基础数据
    x = rand(n,1) * 10;  % U(0, 10) - 均匀分布，范围0到10
    
    % 步骤2：x真值加正常噪声（方差0.01，标准差0.1）
    x_noisy = x + randn(n, 1) * sqrt(0.001);  % x + N(0, 0.01)

    % 步骤3：用加噪声的x计算y（y = 2x - 3）
    y_calculated = x_noisy * true_params(1) + true_params(2);

    % 步骤4：y加正常噪声（方差0.01，标准差0.1）
    y_noisy = y_calculated + randn(n, 1) *sqrt(0.001);  % y + N(0, 0.01)

    % 步骤5：随机选择点加入异常值（10%为3倍，5%为5倍）
    % 正常噪声标准差：sqrt(0.001)
    sigma_x_noise = sqrt(0.001);  % x方向噪声标准差
    sigma_y_noise = sqrt(0.001);  % y方向噪声标准差
    
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
    P_initial(3, :) = 1e12;  % 常数项权重极大值
    
    %% 测试全分量抗差方法
    fprintf('\n--- 测试全分量抗差方法 ---\n');
    try
        tic;
        [X_robust, ~, iter_info_robust] = iterative_weight_optimization_with_timeout(A, L, P_initial, 5.0);
        time_robust = toc;
        
        if time_robust > 5.0
            % 超时情况
            convergence_results.robust_success(end+1) = 0;
            convergence_results.robust_iters(end+1) = NaN;
            convergence_results.robust_params(end+1, :) = [NaN, NaN];
            convergence_results.robust_errors{end+1} = '超时(5秒)';
            convergence_results.robust_times(end+1) = time_robust;
            
            fprintf('✗ 全分量抗差方法超时(>5秒)\n');
            fprintf('  计算时间: %.3f秒\n', time_robust);
        else
            % 正常收敛
            convergence_results.robust_success(end+1) = 1;
            convergence_results.robust_iters(end+1) = iter_info_robust.total_iterations;
            convergence_results.robust_params(end+1, :) = X_robust';
            convergence_results.robust_errors{end+1} = '';
            convergence_results.robust_times(end+1) = time_robust;
            
            fprintf('✓ 全分量抗差方法收敛成功\n');
            fprintf('  迭代次数: %d\n', iter_info_robust.total_iterations);
            fprintf('  计算时间: %.3f秒\n', time_robust);
            fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_robust(1), X_robust(2));
            fprintf('  参数误差: 斜率误差=%.6f, 截距误差=%.6f\n', X_robust(1)-true_params(1), X_robust(2)-true_params(2));
        end
        
    catch ME
        convergence_results.robust_success(end+1) = 0;
        convergence_results.robust_iters(end+1) = NaN;
        convergence_results.robust_params(end+1, :) = [NaN, NaN];
        convergence_results.robust_errors{end+1} = ME.message;
        convergence_results.robust_times(end+1) = NaN;
        
        fprintf('✗ 全分量抗差方法收敛失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
    %% 测试全分量数据探测方法（迭代版本）
    fprintf('\n--- 测试全分量数据探测方法（迭代版本） ---\n');
    try
        tic;
        % 使用与sjtcduibi相同的方法
        % 构建权重矩阵P（3n×3n块对角矩阵格式）
        sigma_x_det = sqrt(0.001);
        sigma_y_det = sqrt(0.001);
        
        % ========== 迭代检测循环 ==========
        m_det = 2;  % 参数个数
        alpha_det = 0.1;  % 显著性水平
        max_iter = 100;  % 最大迭代次数
        iter_count = 0;
        detected_full = false(n, 1);  % 初始化检测结果
        
        % 当前使用的数据
        A_current = A;
        L_current = L;
        n_current = n;
        valid_idx_global = (1:n)';  % 全局有效点索引
        
        while iter_count < max_iter
            iter_count = iter_count + 1;
            
            % 构建当前数据的权重矩阵
            P_newton_full = zeros(3*n_current, 3*n_current);
            Sigma_global_newton = zeros(3*n_current, 3*n_current);
            for i = 1:n_current
                block_start = (i-1)*3 + 1;
                block_end = i*3;
                P_newton_full(block_start:block_end, block_start:block_end) = ...
                    diag([1/sigma_y_det^2, 1/sigma_x_det^2, 1e12]);
                Sigma_global_newton(block_start:block_end, block_start:block_end) = ...
                    diag([sigma_y_det^2, sigma_x_det^2, 1e-12]);
            end
            
            % 调用牛顿法估计参数
            PP_newton = inv(Sigma_global_newton);
            x_hat_newton = TLS_XG_newton3(A_current, L_current, PP_newton);
            
            % 计算临界值
            df1 = 1;
            df2 = n_current - m_det;
            F_critical = sqrt(finv(1 - alpha_det, df1, df2));
            
            % 计算w检验统计量（按照sjtcduibi的方法）
            v_newton = L_current - A_current * x_hat_newton;
            Q_e_newton = Sigma_global_newton;
            k_det = m_det + 1;
            [H_newton, e_A_newton, B_newton, e_newton] = Hessian(A_current, L_current, PP_newton, x_hat_newton);
            
            % 计算残差协方差
            Pv_inv_newton = B_newton * Q_e_newton * B_newton';
            if rcond(Pv_inv_newton) < eps
                P_v_newton = pinv(Pv_inv_newton);
            else
                P_v_newton = inv(Pv_inv_newton);
            end
            
            C_newton = B_newton' * P_v_newton;
            
            % 计算J矩阵
            Gamma_newton = (A_current + e_A_newton)' * P_v_newton;
            if rcond(H_newton) < eps
                dx_dL = -pinv(H_newton) * Gamma_newton;
            else
                dx_dL = -H_newton \ Gamma_newton;
            end
            
            dC_dx = zeros(n_current*k_det*n_current, m_det);
            for i = 1:m_det
                Tk_i = zeros(n_current, n_current * (m_det + 1));
                for j = 1:n_current
                    col_index = (j - 1) * (m_det + 1) + i + 1;
                    Tk_i(j, col_index) = 1;
                end
                dC_dxk_i = Tk_i' * P_v_newton - 2 * B_newton' * P_v_newton * Tk_i * Q_e_newton * B_newton' * P_v_newton;
                dC_dx(:,i) = dC_dxk_i(:);
            end
            dC_dL_vec = dC_dx * dx_dL;
            
            de_dL = zeros(n_current*k_det, n_current);
            for j = 1:n_current
                dC_dL_j = dC_dL_vec(:,j);
                dC_dL_j_reshaped = reshape(dC_dL_j, n_current*k_det, n_current);
                dl = zeros(n_current,1);
                dl(j) = 1;
                de_dL(:, j) = Q_e_newton * (dC_dL_j_reshaped * L_current + C_newton * dl - dC_dL_j_reshaped * A_current * x_hat_newton - C_newton * A_current * dx_dL(:,j));
            end
            J0 = de_dL;
            
            % 计算 Ji = ∂e/∂ai
            Ja = cell(1, m_det);
            dx_da_all = cell(1, m_det);
            
            for param_idx = 1:m_det
                de_da_i = zeros(m_det, n_current);
                
                for obs_idx = 1:n_current
                    R1 = zeros(n_current, m_det);
                    R1(obs_idx, param_idx) = 1;
                    R2 = zeros(n_current, 1);
                    R2(obs_idx) = -x_hat_newton(param_idx);
                    
                    de_da = Q_e_newton * B_newton' * P_v_newton * R2;
                    dET_da = zeros(m_det, n_current);
                    for j = 1:n_current
                        dET_da(:, j) = de_da((m_det+1)*(j-1)+2:(m_det+1)*(j-1)+(m_det+1));
                    end
                    dF_da_single = (R1' + dET_da) * P_v_newton * v_newton + (A_current + e_A_newton)' * P_v_newton * R2;
                    de_da_i(:, obs_idx) = dF_da_single;
                end
                
                if rcond(H_newton) < eps
                    dx_da_i = -pinv(H_newton) * de_da_i;
                else
                    dx_da_i = -H_newton \ de_da_i;
                end
                dx_da_all{param_idx} = dx_da_i;
                
                dC_dx = zeros(n_current*k_det*n_current, m_det);
                for i = 1:m_det
                    Tk_i = zeros(n_current, n_current * (m_det + 1));
                    for j = 1:n_current
                        col_index = (j - 1) * (m_det + 1) + i + 1;
                        Tk_i(j, col_index) = 1;
                    end
                    dC_dxk_i = Tk_i' * P_v_newton - 2 * B_newton' * P_v_newton * Tk_i * Q_e_newton * B_newton' * P_v_newton;
                    dC_dx(:,i) = dC_dxk_i(:);
                end
                dC_dx_vec = dC_dx * dx_da_all{param_idx};
                
                term1 = zeros(k_det*n_current, n_current);
                for j = 1:n_current
                    dC_da_j = dC_dx_vec(:,j);
                    dC_da_j_reshaped = reshape(dC_da_j, n_current*k_det, n_current);
                    da = zeros(n_current,m_det);
                    da(j,param_idx) = 1;
                    term1(:,j) = dC_da_j_reshaped * v_newton - C_newton * da * x_hat_newton;
                end
                de_da = Q_e_newton * (term1 - C_newton * A_current * dx_da_all{param_idx});
                Ja{param_idx} = de_da;
            end
            J_total = [J0, Ja{1}, Ja{2}];
            
            % 构建Q_total
            Q_L_newton = sigma_y_det^2 * eye(n_current);
            QA1_newton = sigma_x_det^2 * eye(n_current);
            QA2_newton = 1e-12 * eye(n_current);
            Q_total = blkdiag(Q_L_newton, QA1_newton, QA2_newton);
            
            % 单位权方差
            sit0_1_newton = (v_newton' * P_v_newton * v_newton) / (n_current - m_det);
            
            % 计算Sigma_e
            Sigma_e = sit0_1_newton * J_total * Q_total * J_total';
            
            % 提取e_L, e_A1, e_A2
            e_hat_reshaped_newton = reshape(e_newton, k_det, n_current)';
            e_L = e_hat_reshaped_newton(:, 1);
            e_A1 = e_hat_reshaped_newton(:, 2);
            
            % 计算w检验统计量
            diag_Sigma_e = diag(Sigma_e);
            row_indices = 1:k_det:k_det*n_current;
            row_indices_A1 = 2:k_det:k_det*n_current;
            
            var_eL = diag_Sigma_e(row_indices) / sit0_1_newton;
            sigma_eL = sqrt(var_eL);
            var_eA1 = diag_Sigma_e(row_indices_A1) / sit0_1_newton;
            sigma_eA1 = sqrt(var_eA1);
            
            w_tests_newton_L = e_L ./ (sqrt(sit0_1_newton) * sigma_eL);
            w_tests_newton_A1 = e_A1 ./ (sqrt(sit0_1_newton) * sigma_eA1);
            
            detected_newton_L = abs(w_tests_newton_L) > F_critical;
            detected_newton_A1 = abs(w_tests_newton_A1) > F_critical;
            
            % 综合判断当前迭代的检测结果
            detected_current_iter = detected_newton_L | detected_newton_A1;
            num_detected_current = sum(detected_current_iter);
            
            fprintf('  迭代 %d: 检测到 %d 个粗差点', iter_count, num_detected_current);
            if num_detected_current > 0
                fprintf(' (全局索引: %s)\n', mat2str(valid_idx_global(detected_current_iter)'));
            else
                fprintf('\n');
            end
            
            % 如果本次迭代没有检测到粗差，退出循环
            if num_detected_current == 0
                X_full_detection = x_hat_newton;
                break;
            end
            
            % 更新全局检测结果
            detected_full(valid_idx_global(detected_current_iter)) = true;
            
            % 剔除检测到的粗差点，准备下一次迭代
            valid_idx_local = find(~detected_current_iter);
            if length(valid_idx_local) < m_det
                % 如果剩余点数不足，使用当前结果并退出
                X_full_detection = x_hat_newton;
                fprintf('  警告：剩余点数不足，停止迭代\n');
                break;
            end
            
            % 更新当前数据
            A_current = A_current(valid_idx_local, :);
            L_current = L_current(valid_idx_local);
            valid_idx_global = valid_idx_global(valid_idx_local);
            n_current = length(valid_idx_local);
        end
        
        % 迭代结束，记录最终结果
        if iter_count >= max_iter
            fprintf('  警告：达到最大迭代次数 %d\n', max_iter);
            X_full_detection = x_hat_newton;
        end
        
        time_full_detection = toc;
        
        % 记录结果
        convergence_results.full_detection_success(end+1) = 1;
        convergence_results.full_detection_iters(end+1) = iter_count;  % 记录实际迭代次数
        convergence_results.full_detection_params(end+1, :) = X_full_detection';
        convergence_results.full_detection_errors{end+1} = '';
        convergence_results.full_detection_times(end+1) = time_full_detection;
        
        fprintf('✓ 全分量数据探测方法（迭代）收敛成功\n');
        fprintf('  迭代次数: %d\n', iter_count);
        fprintf('  计算时间: %.3f秒\n', time_full_detection);
        fprintf('  估计参数: 斜率=%.6f, 截距=%.6f\n', X_full_detection(1), X_full_detection(2));
        fprintf('  参数误差: 斜率误差=%.6f, 截距误差=%.6f\n', X_full_detection(1)-true_params(1), X_full_detection(2)-true_params(2));
        fprintf('  检测到粗差点数: %d (全局)\n', sum(detected_full));
        
    catch ME
        convergence_results.full_detection_success(end+1) = 0;
        convergence_results.full_detection_iters(end+1) = NaN;
        convergence_results.full_detection_params(end+1, :) = [NaN, NaN];
        convergence_results.full_detection_errors{end+1} = ME.message;
        convergence_results.full_detection_times(end+1) = NaN;
        
        fprintf('✗ 全分量数据探测方法失败\n');
        fprintf('  错误信息: %s\n', ME.message);
    end
    
                fprintf('\n');
end

%% 收敛性统计分析
fprintf('========== 收敛性统计分析 ==========\n\n');

% 成功率统计
robust_success_rate = sum(convergence_results.robust_success) / length(convergence_results.robust_success) * 100;
full_detection_success_rate = sum(convergence_results.full_detection_success) / length(convergence_results.full_detection_success) * 100;

fprintf('【收敛成功率】\n');
fprintf('全分量抗差法: %.1f%% (%d/%d)\n', robust_success_rate, sum(convergence_results.robust_success), length(convergence_results.robust_success));
fprintf('全分量数据探测法: %.1f%% (%d/%d)\n', full_detection_success_rate, sum(convergence_results.full_detection_success), length(convergence_results.full_detection_success));

% 迭代次数统计
valid_robust_iters = convergence_results.robust_iters(~isnan(convergence_results.robust_iters));
valid_full_detection_iters = convergence_results.full_detection_iters(~isnan(convergence_results.full_detection_iters));

fprintf('\n【迭代次数统计】\n');
if ~isempty(valid_robust_iters)
    fprintf('全分量抗差法: %.1f ± %.1f (范围: %d-%d)\n', mean(valid_robust_iters), std(valid_robust_iters), min(valid_robust_iters), max(valid_robust_iters));
end
if ~isempty(valid_full_detection_iters)
    fprintf('全分量数据探测法: %.1f ± %.1f (范围: %d-%d)\n', mean(valid_full_detection_iters), std(valid_full_detection_iters), min(valid_full_detection_iters), max(valid_full_detection_iters));
end

% 平均运行时间统计（过滤异常值）
valid_robust_times = convergence_results.robust_times(~isnan(convergence_results.robust_times));
valid_full_detection_times = convergence_results.full_detection_times(~isnan(convergence_results.full_detection_times));

% 过滤掉异常长的时间（使用中位数±3倍MAD方法，更稳健）
fprintf('\n【平均运行时间统计 (已过滤异常值)】\n');

if ~isempty(valid_robust_times)
    median_robust = median(valid_robust_times);
    mad_robust = median(abs(valid_robust_times - median_robust));
    % 过滤：保留中位数±3*1.4826*MAD范围内的值
    filtered_robust_times = valid_robust_times(abs(valid_robust_times - median_robust) <= 3 * 1.4826 * mad_robust);
    fprintf('全分量抗差法: %.4f ± %.4f秒 (范围: %.4f-%.4f秒) [使用%d/%d个样本]\n', ...
        mean(filtered_robust_times), std(filtered_robust_times), min(filtered_robust_times), max(filtered_robust_times), ...
        length(filtered_robust_times), length(valid_robust_times));
end

if ~isempty(valid_full_detection_times)
    median_full_detection = median(valid_full_detection_times);
    mad_full_detection = median(abs(valid_full_detection_times - median_full_detection));
    filtered_full_detection_times = valid_full_detection_times(abs(valid_full_detection_times - median_full_detection) <= 3 * 1.4826 * mad_full_detection);
    fprintf('全分量数据探测法: %.4f ± %.4f秒 (范围: %.4f-%.4f秒) [使用%d/%d个样本]\n', ...
        mean(filtered_full_detection_times), std(filtered_full_detection_times), min(filtered_full_detection_times), max(filtered_full_detection_times), ...
        length(filtered_full_detection_times), length(valid_full_detection_times));
end

% 参数估计精度统计
valid_robust_params = convergence_results.robust_params(convergence_results.robust_success == 1, :);
valid_full_detection_params = convergence_results.full_detection_params(convergence_results.full_detection_success == 1, :);

if ~isempty(valid_robust_params)
    robust_slope_rmse = sqrt(mean((valid_robust_params(:,1) - true_params(1)).^2));
    robust_intercept_rmse = sqrt(mean((valid_robust_params(:,2) - true_params(2)).^2));
    fprintf('\n【参数估计精度 - 全分量抗差法】\n');
    fprintf('斜率RMSE: %.6f, 截距RMSE: %.6f\n', robust_slope_rmse, robust_intercept_rmse);
end

if ~isempty(valid_full_detection_params)
    full_detection_slope_rmse = sqrt(mean((valid_full_detection_params(:,1) - true_params(1)).^2));
    full_detection_intercept_rmse = sqrt(mean((valid_full_detection_params(:,2) - true_params(2)).^2));
    fprintf('\n【参数估计精度 - 全分量数据探测法】\n');
    fprintf('斜率RMSE: %.6f, 截距RMSE: %.6f\n', full_detection_slope_rmse, full_detection_intercept_rmse);
end

% 失败原因分析
fprintf('\n【失败原因分析】\n');
robust_failures = find(convergence_results.robust_success == 0);
full_detection_failures = find(convergence_results.full_detection_success == 0);

fprintf('=== 全分量抗差法 ===\n');
if ~isempty(robust_failures)
    fprintf('失败实验: %s\n', mat2str(robust_failures));
    for i = 1:min(length(robust_failures), 5)  % 只显示前5个
        fprintf('  实验%d: %s\n', robust_failures(i), convergence_results.robust_errors{robust_failures(i)});
    end
else
    fprintf('所有实验均成功收敛\n');
end

fprintf('\n=== 全分量数据探测法 ===\n');
if ~isempty(full_detection_failures)
    fprintf('失败实验: %s\n', mat2str(full_detection_failures));
    for i = 1:min(length(full_detection_failures), 5)
        fprintf('  实验%d: %s\n', full_detection_failures(i), convergence_results.full_detection_errors{full_detection_failures(i)});
    end
else
    fprintf('所有实验均成功收敛\n');
end

%% 参数分布可视化
if ~isempty(valid_robust_params) || ~isempty(valid_full_detection_params)
    fprintf('\n========== 参数分布可视化 ==========\n');
    
    % 设置全局字体为Times New Roman
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    
    % 创建一个图框，包含两个子图
    figure('Position', [100, 100, 1600, 600]);
    
    % === 子图1：斜率参数分布（直方图+密度图） ===
    subplot(1, 2, 1);
    hold on;
    
    % 设置双y轴
    yyaxis left;
    
    % 绘制直方图
    if ~isempty(valid_robust_params)
        histogram(valid_robust_params(:,1), 30, 'Normalization', 'pdf', ...
                 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                 'DisplayName', 'Robust (Histogram)');
    end
    
    if ~isempty(valid_full_detection_params)
        histogram(valid_full_detection_params(:,1), 30, 'Normalization', 'pdf', ...
                 'FaceColor', [0, 0, 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                 'DisplayName', 'Detection (Histogram)');
    end
    
    ylabel('Probability Density', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ax = gca;
    ax.YColor = 'k';
    
    % 切换到右侧y轴绘制密度曲线
    yyaxis right;
    
    if ~isempty(valid_robust_params)
        % 全分量抗差法斜率核密度估计（红色粗线）
        [f_robust_slope, xi_robust_slope] = ksdensity(valid_robust_params(:,1));
        plot(xi_robust_slope, f_robust_slope, 'r-', 'LineWidth', 3.5, 'DisplayName', 'Robust (KDE)');
    end
    
    if ~isempty(valid_full_detection_params)
        % 全分量数据探测法斜率核密度估计（蓝色细线）
        [f_full_det_slope, xi_full_det_slope] = ksdensity(valid_full_detection_params(:,1));
        plot(xi_full_det_slope, f_full_det_slope, 'b-', 'LineWidth', 2, 'DisplayName', 'Detection (KDE)');
    end
    
    ylabel('Kernel Density', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ax.YColor = 'k';
    
    % 真实值标记
    yyaxis left;
    xline(true_params(1), 'k--', 'LineWidth', 2.5, 'DisplayName', 'True Value');
    
    xlabel('Slope Estimate', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    title('Slope Parameter Distribution', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    legend('Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');
    grid on;
    set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
    hold off;
    
    % === 子图2：截距参数分布（直方图+密度图） ===
    subplot(1, 2, 2);
    hold on;
    
    % 设置双y轴
    yyaxis left;
    
    % 绘制直方图
    if ~isempty(valid_robust_params)
        histogram(valid_robust_params(:,2), 30, 'Normalization', 'pdf', ...
                 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                 'DisplayName', 'Robust (Histogram)');
    end
    
    if ~isempty(valid_full_detection_params)
        histogram(valid_full_detection_params(:,2), 30, 'Normalization', 'pdf', ...
                 'FaceColor', [0, 0, 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                 'DisplayName', 'Detection (Histogram)');
    end
    
    ylabel('Probability Density', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ax = gca;
    ax.YColor = 'k';
    
    % 切换到右侧y轴绘制密度曲线
    yyaxis right;
    
    if ~isempty(valid_robust_params)
        % 全分量抗差法截距核密度估计（红色粗线）
        [f_robust_intercept, xi_robust_intercept] = ksdensity(valid_robust_params(:,2));
        plot(xi_robust_intercept, f_robust_intercept, 'r-', 'LineWidth', 3.5, 'DisplayName', 'Robust (KDE)');
    end
    
    if ~isempty(valid_full_detection_params)
        % 全分量数据探测法截距核密度估计（蓝色细线）
        [f_full_det_intercept, xi_full_det_intercept] = ksdensity(valid_full_detection_params(:,2));
        plot(xi_full_det_intercept, f_full_det_intercept, 'b-', 'LineWidth', 2, 'DisplayName', 'Detection (KDE)');
    end
    
    ylabel('Kernel Density', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ax.YColor = 'k';
    
    % 真实值标记
    yyaxis left;
    xline(true_params(2), 'k--', 'LineWidth', 2.5, 'DisplayName', 'True Value');
    
    xlabel('Intercept Estimate', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    title('Intercept Parameter Distribution', 'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    legend('Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');
    grid on;
    set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
    hold off;
    
    % 保存图片（分辨率300 DPI）
    print(gcf, '两种方法参数分布对比.png', '-dpng', '-r300');
    fprintf('参数分布对比图已保存为: 两种方法参数分布对比.png (300 DPI)\n');
end

%% 结论
fprintf('\n========== 诊断结论 ==========\n');

fprintf('\n=== 全分量抗差方法（迭代权重优化） ===\n');
if robust_success_rate >= 80
    fprintf('✓ 全分量抗差方法在随机粗差下表现良好\n');
    fprintf('  成功率达到 %.1f%%，数值稳定性较好\n', robust_success_rate);
else
    fprintf('⚠️  全分量抗差方法在随机粗差下存在收敛问题\n');
    fprintf('   成功率仅为 %.1f%%，建议调整抗差函数参数或优化算法\n', robust_success_rate);
end

fprintf('\n=== 全分量数据探测方法（迭代检测） ===\n');
if full_detection_success_rate >= 80
    fprintf('✓ 全分量数据探测方法在随机粗差下表现良好\n');
    fprintf('  成功率达到 %.1f%%，数值稳定性较好\n', full_detection_success_rate);
else
    fprintf('⚠️  全分量数据探测方法在随机粗差下存在问题\n');
    fprintf('   成功率仅为 %.1f%%\n', full_detection_success_rate);
end

fprintf('\n=== 方法对比总结 ===\n');
if ~isempty(valid_robust_times) && ~isempty(valid_full_detection_times)
    fprintf('平均运行时间对比：\n');
    fprintf('  - 全分量抗差法: %.4f秒\n', mean(filtered_robust_times));
    fprintf('  - 全分量数据探测法: %.4f秒\n', mean(filtered_full_detection_times));
    speedup = mean(filtered_robust_times) / mean(filtered_full_detection_times);
    if speedup > 1
        fprintf('  - 数据探测法快 %.2f倍\n', speedup);
    else
        fprintf('  - 抗差法快 %.2f倍\n', 1/speedup);
    end
end

if ~isempty(valid_robust_params) && ~isempty(valid_full_detection_params)
    fprintf('\n参数估计精度对比：\n');
    fprintf('  - 全分量抗差法: 斜率RMSE=%.6f, 截距RMSE=%.6f\n', robust_slope_rmse, robust_intercept_rmse);
    fprintf('  - 全分量数据探测法: 斜率RMSE=%.6f, 截距RMSE=%.6f\n', full_detection_slope_rmse, full_detection_intercept_rmse);
end


%% 算法函数定义（与原实验完全一致）

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
param_tol_abs = 5e-3;  % 绝对误差阈值，考虑参数理论精度
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

function [x_tls, e_hat, iter] = TLS_XG_newton3(A_obs, L_obs, P)
% TLS_XG_newton3 使用牛顿法求解加权总体最小二乘问题

    tol = 1e-10;
    % 获取问题维度
    [n, m] = size(A_obs);
    k = m + 1;
    
    % 只计算一次Q_e和I_n，不预提取块
    Q_e = inv(P);
    I_n = eye(n);
    
    % 提取L对应的权矩阵用于初始OLS估计
    P_L = zeros(n, n);
    for i = 1:n
        P_L(i, i) = P((i-1)*k+1, (i-1)*k+1);
    end
    
    % 初始OLS估计
    x0 = (A_obs' * P_L * A_obs) \ (A_obs' * P_L * L_obs);
    
    % TLS牛顿迭代法
    x = x0;
    iter = 0;
    dx_norm = 1;

    while dx_norm > tol
        iter = iter + 1;
        
        % 1. 计算残差向量
        v = L_obs - A_obs * x;
        
        % 2. 构建b向量和B矩阵
        b = [1, x'];
        B = kron(I_n, b);
        
        % 3. 计算 P_v = (B * Q_e * B')^{-1}
        Pv_inv = B * Q_e * B';
        if rcond(Pv_inv) < eps
            P_v = pinv(Pv_inv);
        else
            P_v = inv(Pv_inv);
        end
        
        % 预计算常用项
        vp = P_v * v;
        B_T_vp = B' * vp;
        
        % 4. 计算误差估计 ê
        e_hat = Q_e * B_T_vp;
        
        % 5. 提取系数矩阵的误差
        e_hat_reshaped = reshape(e_hat, k, n)';
        e_A = e_hat_reshaped(:, 2:end);
        
        % 6. 计算修正后的系数矩阵
        A_corr = A_obs + e_A;
        A_corr2 = A_obs + 2*e_A;
        % 7. 计算梯度 F
        F = A_corr' * vp;
        
        % 8. 预计算Hessian矩阵的公共项
        P_v_e_A = P_v * e_A;
        P_v_A_obs = P_v * A_obs;
        
        H1 = - A_corr' * P_v *A_corr2;

        % H5矩阵 - 优化计算（移除冗余循环）
        H5 = zeros(m, m);
        
        for i = 1:m
            % 计算 v_tk
            v_tk = zeros(k*n, 1);
            v_tk(i+1:k:(k*n)) = vp;

            % 预计算kron项
            E_bk = kron(P_v_e_A(:, i)', b);
            A_bk = kron(P_v_A_obs(:, i)', b);
            
            % 计算de/dxi
            de_dxi = Q_e * (v_tk - 2*E_bk' - A_bk');

            % 提取dE^T/dxi（向量化）
            de_dxi_reshaped = reshape(de_dxi, k, n);
            dET_dxi = de_dxi_reshaped(2:k, :);
            
            % 计算H5的第i列
            H5(:, i) = dET_dxi * vp;
        end

        H = H1 + H5;
        
        % 9. 牛顿迭代更新
        if rcond(H) < eps
            dx = pinv(H) * F;
        else
            dx = H \ F;
        end
        
        x = x - dx;
        dx_norm = norm(dx);
    end
    
    x_tls = x; 
    
end

function [H, e_A, B, e_hat] = Hessian(A, L, P, x)
% Hessian 计算Hessian矩阵及相关误差项

    % 获取问题维度
    [n, m] = size(A);
    k = m + 1;
    
    % 只计算一次Q_e和I_n，不预提取块
    Q_e = pinv(P);
    I_n = eye(n);

    % 1. 计算残差向量
    v = L - A * x;
        
    % 2. 构建b向量和B矩阵
    b = [1, x'];
    B = kron(I_n, b);
        
    % 3. 计算 P_v = (B * Q_e * B')^{-1}
    Pv_inv = B * Q_e * B';
    if rcond(Pv_inv) < eps
        P_v = pinv(Pv_inv);
    else
        P_v = inv(Pv_inv);
    end
        
    % 预计算常用项
    vp = P_v * v;
    B_T_vp = B' * vp;
        
    % 4. 计算误差估计 ê
    e_hat = Q_e * B_T_vp;
        
    % 5. 提取系数矩阵的误差
    e_hat_reshaped = reshape(e_hat, k, n)';
    e_A = e_hat_reshaped(:, 2:end);
        
    % 6. 计算修正后的系数矩阵
    A_corr = A + e_A;
    A_corr2 = A + 2*e_A;
        
    % 8. 预计算Hessian矩阵的公共项
    P_v_e_A = P_v * e_A;
    P_v_A_obs = P_v * A;
        
    H1 = - A_corr' * P_v *A_corr2;

    % H5矩阵 - 优化计算（移除冗余循环）
    H5 = zeros(m, m);
        
    for i = 1:m
        % 计算 v_tk
        v_tk = zeros(k*n, 1);
        v_tk(i+1:k:(k*n)) = vp;
        % 预计算kron项
        E_bk = kron(P_v_e_A(:, i)', b);
        A_bk = kron(P_v_A_obs(:, i)', b);
            
        % 计算de/dxi
        de_dxi = Q_e * (v_tk - 2*E_bk' - A_bk');
        % 提取dE^T/dxi（向量化）
        de_dxi_reshaped = reshape(de_dxi, k, n);
        dET_dxi = de_dxi_reshaped(2:k, :);
           
        % 计算H5的第i列
        H5(:, i) = dET_dxi * vp;
    end

    H = H1 + H5;

end
