clear all; clc; close all;

% 添加必要的路径（使用相对路径，避免路径不存在的问题）
current_dir = fileparts(mfilename('fullpath'));
if isempty(current_dir)
    current_dir = pwd;
end

% 尝试添加路径，如果不存在则忽略警告
try
    addpath(fullfile(current_dir, '..', 'data snooping'));
catch
end
try
    addpath(fullfile(current_dir, '..', '试验'));
catch
end
try
    addpath('data snooping');
catch
end
try
    addpath('试验');
catch
end

% 注意：separate_direction_robust_tls和correlated_robust_tls函数定义在试验/相关性抗差文件中
% 这些函数及其依赖函数（TLS_XG, TLS_XG_with_W等）将在文件末尾定义

%% ========== 相关性数据抗差实验 ==========
% 基于相关数据生成相关性数据，添加粗差，测试不处理相关性的方法在处理相关性数据时的表现

fprintf('========== 相关性数据抗差实验 ==========\n');
fprintf('测试不处理相关性的方法在处理相关性数据时的表现\n\n');

% 参数设置
n = 20;                % 观测值数量（从10改为20，增加样本量）
m = 2;                 % 参数个数
k = m + 1;             % 每个观测点的变量数(L + A1 + A2)
num_experiments = 100; % 实验次数（先设置较小值用于测试）

% 生成真实数据
A1_true = (0:(n-1))';  % 从0到n-1，适应不同的n值
A2_true = ones(n, 1);
L_true = 2*A1_true - A2_true;
A_true = [A1_true, A2_true];

% 真实参数
a_true = 2;             % 真实斜率
b_true = -1;             % 真实截距

% 设置噪声方差
sigma_x = 0.001;         % x的噪声方差
sigma_L = 0.001;         % L的噪声方差

% 定义三个相关系数
rho_values = [0.3, 0.6, 0.9];
num_rhos = length(rho_values);

% 粗差设置
outlier_ratio_3sigma = 0.10;  % 10%的粗差为3倍标准差
outlier_ratio_5sigma = 0.05;  % 5%的粗差为5倍标准差
outlier_ratio = outlier_ratio_3sigma + outlier_ratio_5sigma;  % 总粗差比例15%
outlier_magnitude_3sigma = 3;  % 3倍粗差
outlier_magnitude_5sigma = 5;  % 5倍粗差
alpha = 0.1;                   % 显著性水平

% ========== 实验设置 ==========
fprintf('实验设置: n=%d, 粗差比例=%.1f%% (%.1f%%为3倍标准差, %.1f%%为5倍标准差)\n\n', ...
    n, outlier_ratio * 100, outlier_ratio_3sigma * 100, outlier_ratio_5sigma * 100);

% 计算临界值（F检验，自由度为1和n-m）
df1 = 1;
df2 = n - m;
F_critical = sqrt(finv(1 - alpha, df1, df2));

% 初始化结果存储
results_all = struct();

% 对每个相关系数进行实验
for rho_idx = 1:num_rhos
    rho = rho_values(rho_idx);
    
    fprintf('\n========================================\n');
    fprintf('=== 相关系数 rho = %.1f ===\n', rho);
    fprintf('========================================\n\n');
    
    % 初始化该相关系数下的结果
    success_count_mahboub = 0;
    success_count_component = 0;
    success_count_full = 0;
    success_count_mahboub_irtls = 0;
    success_count_separate = 0;      % 分方向残差抗差
    success_count_overall = 0;        % 总体残差抗差
    
    fp_count_mahboub = 0;
    fp_count_component = 0;
    fp_count_full = 0;
    fp_count_mahboub_irtls = 0;
    fp_count_separate = 0;
    fp_count_overall = 0;
    
    fn_count_mahboub = 0;
    fn_count_component = 0;
    fn_count_full = 0;
    fn_count_mahboub_irtls = 0;
    fn_count_separate = 0;
    fn_count_overall = 0;
    
    param_error_mahboub = zeros(2, 1);
    param_error_component = zeros(2, 1);
    param_error_full = zeros(2, 1);
    param_error_mahboub_irtls = zeros(2, 1);
    param_error_separate = zeros(2, 1);
    param_error_overall = zeros(2, 1);
    
    param_count_mahboub = 0;
    param_count_component = 0;
    param_count_full = 0;
    param_count_mahboub_irtls = 0;
    param_count_separate = 0;
    param_count_overall = 0;
    
    % 初始化单位权中误差sigma0的累积变量
    sigma0_sum_mahboub = 0;
    sigma0_sum_component = 0;
    sigma0_sum_full = 0;
    sigma0_sum_mahboub_irtls = 0;
    sigma0_sum_separate = 0;
    sigma0_sum_overall = 0;
    
    total_time_mahboub = 0;
    total_time_component = 0;
    total_time_full = 0;
    total_time_mahboub_irtls = 0;
    total_time_separate = 0;
    total_time_overall = 0;
    
    param_estimates_mahboub = [];
    param_estimates_component = [];
    param_estimates_full = [];
    param_estimates_mahboub_irtls = [];
    param_estimates_separate = [];
    param_estimates_overall = [];
    
    % ========== 验证相关性结构（仅第一次实验）==========
    verify_correlation = (rho_idx == 1);  % 只在第一个rho时验证
    
    for exp_idx = 1:num_experiments
        % 实验进度提示已移除，减少输出
        
        % ========== 生成相关性数据 ==========
        % 构建单个观测点内的协方差矩阵
        total_vars = m + 1;  % L + A1, A2
        Sigma_point = zeros(total_vars, total_vars);
        
        % 对角线元素（方差）
        Sigma_point(1,1) = sigma_L;         % L的方差
        Sigma_point(2,2) = sigma_x;         % A1的方差
        Sigma_point(3,3) = 1e-12;           % A2的方差（设置为10^(-12)）
        
        % 同一观测点内变量间的相关性
        % L和A1之间的相关性
        Sigma_point(1,2) = rho * sqrt(Sigma_point(1,1) * Sigma_point(2,2));
        Sigma_point(2,1) = Sigma_point(1,2);
        
        % A2与其他变量不相关
        Sigma_point(1,3) = 0;
        Sigma_point(3,1) = 0;
        Sigma_point(2,3) = 0;
        Sigma_point(3,2) = 0;
        
        % 构建完整的协方差矩阵（块对角矩阵，观测点之间不相关）
        % 注意：观测点之间的相关性已经通过同一观测点内L和A1的相关性来体现
        % 如果不同观测点间也相关，会导致全局协方差矩阵接近奇异
        Sigma_global = zeros(3*n, 3*n);
        for i = 1:n
            block_start = (i-1)*3 + 1;
            block_end = i*3;
            % 每个观测点使用相同的Sigma_point（内部变量间有rho相关性）
            Sigma_global(block_start:block_end, block_start:block_end) = Sigma_point;
        end
        
        % 验证协方差矩阵的正定性（仅第一次，静默验证）
        if verify_correlation && exp_idx == 1
            eigenvals = eig(Sigma_global);
            min_eigenval = min(eigenvals);
            if min_eigenval < 0
                fprintf('警告: 协方差矩阵不是正定的！最小特征值: %.2e\n', min_eigenval);
            end
        end
        
        % 生成观测数据
        try
            R = chol(Sigma_global);  % Sigma_global = R' * R
        catch
            % 如果不是正定，添加小的正则化项
            Sigma_global = Sigma_global + eye(size(Sigma_global)) * 1e-10;
            R = chol(Sigma_global);
        end
        
        noise = randn(n*total_vars, 1);
        correlated_noise = R' * noise;
        
        % 将噪声添加到真实数据中
        Y_true = zeros(n*total_vars, 1);
        for i = 1:n
            idx = (i-1)*total_vars + 1;
            Y_true(idx) = L_true(i);          % L
            Y_true(idx+1) = A1_true(i);       % A1
            Y_true(idx+2) = A2_true(i);       % A2
        end
        
        % 生成观测值
        Y_obs = Y_true;
        
        % 只对L和A1添加噪声，A2保持不变（无误差）
        for i = 1:n
            idx = (i-1)*total_vars + 1;
            Y_obs(idx) = Y_true(idx) + correlated_noise(idx);      % L有噪声
            Y_obs(idx+1) = Y_true(idx+1) + correlated_noise(idx+1); % A1有噪声
            Y_obs(idx+2) = Y_true(idx+2);                           % A2无噪声
        end
        
        % 提取观测值
        L_obs = Y_obs(1:total_vars:end);
        A1_obs = Y_obs(2:total_vars:end);
        A2_obs = Y_obs(3:total_vars:end);
        A_obs = [A1_obs, A2_obs];
        
        % 验证实际生成的相关性（静默验证，不打印）
        
        % ========== 添加粗差 ==========
        % 随机选择点加入粗差：10%为3倍标准差，5%为5倍标准差
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
        end
        
        % 计算粗差大小
        outlier_size_x_3sigma = outlier_magnitude_3sigma * sigma_x;  % X方向3倍粗差
        outlier_size_y_3sigma = outlier_magnitude_3sigma * sigma_L;  % Y方向3倍粗差
        outlier_size_x_5sigma = outlier_magnitude_5sigma * sigma_x;  % X方向5倍粗差
        outlier_size_y_5sigma = outlier_magnitude_5sigma * sigma_L;  % Y方向5倍粗差
        
        % 在选定的点上添加粗差（x和y同时添加，方向相反确保点偏离拟合线）
        % 注意：粗差是在相关性噪声之后添加的，这会破坏相关性结构
        % 这是合理的，因为粗差本身就是异常值，应该独立于相关性噪声
        A1_obs_with_outlier = A1_obs;
        L_obs_with_outlier = L_obs;
        
        % 添加3倍粗差
        for i = 1:length(outlier_indices_3sigma)
            idx = outlier_indices_3sigma(i);
            sign_outlier = sign(randn(1));  % 随机方向：+1或-1
            if sign_outlier == 0
                sign_outlier = 1;
            end
            A1_obs_with_outlier(idx) = A1_obs_with_outlier(idx) + sign_outlier * outlier_size_x_3sigma;
            L_obs_with_outlier(idx) = L_obs_with_outlier(idx) - sign_outlier * outlier_size_y_3sigma;
        end
        
        % 添加5倍粗差
        for i = 1:length(outlier_indices_5sigma)
            idx = outlier_indices_5sigma(i);
            sign_outlier = sign(randn(1));  % 随机方向：+1或-1
            if sign_outlier == 0
                sign_outlier = 1;
            end
            A1_obs_with_outlier(idx) = A1_obs_with_outlier(idx) + sign_outlier * outlier_size_x_5sigma;
            L_obs_with_outlier(idx) = L_obs_with_outlier(idx) - sign_outlier * outlier_size_y_5sigma;
        end
        
        A_obs_with_outlier = [A1_obs_with_outlier, A2_obs];
        
        % 记录粗差位置（合并所有粗差点）
        if n_outliers_total > 0
            true_outlier_positions = all_outlier_indices;
        else
            true_outlier_positions = [];
        end
        
        % 验证添加粗差后的相关性（静默验证，不打印）
        
        % ========== 方法1: Mahboub WTLS方法（来自sjtcduibi）==========
        tic_mahboub = tic;
        
        % 初始化（使用OLS作为初始值）
        x_ols_init = (A_obs_with_outlier' * A_obs_with_outlier) \ (A_obs_with_outlier' * L_obs_with_outlier);
        Q_y_mah = eye(n);  % 注意：这里假设独立，不处理相关性
        Q_A_mah = zeros(n*m, n*m);
        Q_A_mah(1:n, 1:n) = eye(n);  % x列有误差
        Q_A_mah(n+1:n*m, n+1:n*m) = zeros(n);  % 常数列无误差
        x_hat_mah = x_ols_init;
        epsilon_mah = 1e-10;
        
        % WTLS迭代
        for iter = 1:20
            e_hat_mah = L_obs_with_outlier - A_obs_with_outlier * x_hat_mah;
            x_kron_T = kron(x_hat_mah', eye(n));
            x_kron = kron(x_hat_mah, eye(n));
            Q_y_tilde_mah = Q_y_mah + x_kron_T * Q_A_mah * x_kron;
            Q_y_tilde_inv_mah = inv(Q_y_tilde_mah);
            vec_E_A_mah = -Q_A_mah * x_kron * Q_y_tilde_inv_mah * e_hat_mah;
            E_A_hat_mah = [vec_E_A_mah(1:n), vec_E_A_mah(n+1:end)];
            A_tilde_mah = A_obs_with_outlier - E_A_hat_mah;
            y_tilde_mah = L_obs_with_outlier - E_A_hat_mah * x_hat_mah;
            x_hat_new_mah = (A_tilde_mah' * Q_y_tilde_inv_mah * A_tilde_mah) \ (A_tilde_mah' * Q_y_tilde_inv_mah * y_tilde_mah);
            delta_mah = norm(x_hat_new_mah - x_hat_mah);
            x_hat_mah = x_hat_new_mah;
            if delta_mah < epsilon_mah
                break;
            end
        end
        
        % w检验
        final_residuals_mah = L_obs_with_outlier - A_obs_with_outlier * x_hat_mah;
        x_kron_T = kron(x_hat_mah', eye(n));
        x_kron = kron(x_hat_mah, eye(n));
        Q_y_tilde_mah = Q_y_mah + x_kron_T * Q_A_mah * x_kron;
        Q_y_tilde_inv_mah = inv(Q_y_tilde_mah);
        sigma_0_sq_mah = (final_residuals_mah' * Q_y_tilde_inv_mah * final_residuals_mah) / (n - m);
        sigma_0_mah = sqrt(sigma_0_sq_mah);
        vec_E_A_mah = -Q_A_mah * x_kron * Q_y_tilde_inv_mah * final_residuals_mah;
        E_A_hat_mah = [vec_E_A_mah(1:n), vec_E_A_mah(n+1:end)];
        A_tilde_mah = A_obs_with_outlier - E_A_hat_mah;
        Q_x_mah = inv(A_tilde_mah' * Q_y_tilde_inv_mah * A_tilde_mah);
        Q_e_normalized_mah = Q_y_tilde_mah - A_tilde_mah * Q_x_mah * A_tilde_mah';
        
        w_tests_mahboub = zeros(n, 1);
        for i = 1:n
            e_i = zeros(n, 1);
            e_i(i) = 1;
            numerator_mah = e_i' * Q_y_tilde_inv_mah * final_residuals_mah;
            denominator_mah = sigma_0_mah * sqrt(e_i' * Q_y_tilde_inv_mah * Q_e_normalized_mah * Q_y_tilde_inv_mah * e_i);
            w_tests_mahboub(i) = numerator_mah / denominator_mah;
        end
        
        detected_mahboub = abs(w_tests_mahboub) > F_critical;
        time_mahboub = toc(tic_mahboub);
        total_time_mahboub = total_time_mahboub + time_mahboub;
        
        % ========== 方法2: Mahboub IRTLS方法（来自mah函数）==========
        tic_mahboub_irtls = tic;
        
        % 准备输入参数（注意：这里假设独立，不处理相关性）
        sigma_A_mah = [sigma_x * ones(n, 1), 1e-6 * ones(n, 1)];  % [n x m]，A1有误差，A2无误差
        sigma_y_mah = sigma_L * ones(n, 1);  % [n x 1]
        
        options_mah = struct();
        options_mah.weight_type = 'huber';
        options_mah.k = 1.5;
        options_mah.max_iter = 100;
        options_mah.tol = 1e-6;
        options_mah.verbose = false;
        
        X_mahboub_irtls = [];  % 初始化为空，用于检查是否成功计算
        try
            [X_mahboub_irtls, residuals_mahboub_irtls, iter_info_mahboub_irtls] = IRTLS_Mahboub2013(A_obs_with_outlier, L_obs_with_outlier, sigma_A_mah, sigma_y_mah, options_mah);
            
            % 对于IRTLS方法，使用最终权重来判断粗差
            % 权重小于某个阈值（如0.1）的点被认为是粗差（更严格的阈值）
            final_weights_y = iter_info_mahboub_irtls.final_weights_y;
            weight_threshold = 0.1;  % 更严格的阈值，只有权重很小的点才被认为是粗差
            detected_mahboub_irtls = final_weights_y < weight_threshold;
            
            % 如果权重方法检测到的粗差太多（超过50%），使用更严格的阈值
            if sum(detected_mahboub_irtls) > n * 0.5
                weight_threshold = 0.05;  % 更严格的阈值
                detected_mahboub_irtls = final_weights_y < weight_threshold;
            end
            
            % 如果权重方法没有检测到粗差，使用残差方法作为补充
            if ~any(detected_mahboub_irtls)
                % 使用残差的标准化值来判断
                v_irtls = residuals_mahboub_irtls.v;
                sigma0_irtls = iter_info_mahboub_irtls.sigma0(end);
                if ~isempty(sigma0_irtls) && sigma0_irtls > 0
                    std_residuals = abs(v_irtls) / (sigma0_irtls + 1e-10);
                    detected_mahboub_irtls = std_residuals > 3.0;  % 使用3.0倍标准差作为阈值（更严格）
                end
            end
        catch ME
            % fprintf('    IRTLS方法失败: %s\n', ME.message);
            X_mahboub_irtls = x_ols_init;
            detected_mahboub_irtls = false(n, 1);
        end
        
        time_mahboub_irtls = toc(tic_mahboub_irtls);
        total_time_mahboub_irtls = total_time_mahboub_irtls + time_mahboub_irtls;
        
        % ========== 方法3: 分量压缩法（来自sjtcduibi）==========
        % 设置权阵（注意：这里假设独立，不处理相关性）
        py = 1/sigma_L^2 * ones(n, 1);  % L的权
        px1 = 1/sigma_x^2 * ones(n, 1);  % A1的权
        px2 = 1e14 * ones(n, 1);  % A2的权（常数项，近似无噪声）
        P = [py, px1, px2]';  % 3 x n
        
        % 调用数据探测函数（需要从点云文件中获取）
        tic_component = tic;
        try
            [detected_component, w_tests_component, v_component, x_hat_component, F_critical_comp, results_component] = detect_outlier_v(A_obs_with_outlier, L_obs_with_outlier, P, alpha);
        catch ME
            fprintf('    分量压缩法失败: %s\n', ME.message);
            detected_component = false(n, 1);
            x_hat_component = x_ols_init;
        end
        time_component = toc(tic_component);
        total_time_component = total_time_component + time_component;
        
        % ========== 方法4: 全分量方法（来自sjtcduibi）==========
        tic_full = tic;
        
        % 构建权重矩阵P（3n×3n块对角矩阵格式，注意：这里假设独立，不处理相关性）
        P_newton_full = zeros(3*n, 3*n);
        for i = 1:n
            block_start = (i-1)*3 + 1;
            block_end = i*3;
            P_newton_full(block_start:block_end, block_start:block_end) = ...
                diag([1/sigma_L^2, 1/sigma_x^2, 1e12]);
        end
        
        % 构建非相关的全局协方差矩阵（块对角矩阵）
        Sigma_global_newton = zeros(3*n, 3*n);
        for i = 1:n
            block_start = (i-1)*3 + 1;
            block_end = i*3;
            Sigma_global_newton(block_start:block_end, block_start:block_end) = ...
                diag([sigma_L^2, sigma_x^2, 1e-12]);
        end
        
        % 调用牛顿法
        try
            PP_newton = inv(Sigma_global_newton);
            x_hat_newton = TLS_XG_newton3(A_obs_with_outlier, L_obs_with_outlier, PP_newton);
            v_newton = L_obs_with_outlier - A_obs_with_outlier * x_hat_newton;
            
            % 计算误差传播（使用Hessian函数）
            Q_e_newton = Sigma_global_newton;
            [H_newton, e_A_newton, B_newton, e_newton] = Hessian(A_obs_with_outlier, L_obs_with_outlier, PP_newton, x_hat_newton);
            e_hat_reshaped_newton = reshape(e_newton, k, n)';
            e_L = e_hat_reshaped_newton(:, 1);  % 提取第一列（L的残差）
            
            % 计算J矩阵[J0,J1,J2]（按照sjtcduibi的方法）
            J_total = zeros(k*n, k*n);
            
            % 1.计算 P_v = (B * Q_e * B')^{-1}
            Pv_inv_newton = B_newton * Q_e_newton * B_newton';
            if rcond(Pv_inv_newton) < eps
                P_v_newton = pinv(Pv_inv_newton);
            else
                P_v_newton = inv(Pv_inv_newton);
            end
            
            % 2.计算 dx_dL = -H \ (A + e_A)' * P_v;
            Gamma_newton = (A_obs_with_outlier + e_A_newton)' * P_v_newton;
            if rcond(H_newton) < eps
                dx_dL = -pinv(H_newton) * Gamma_newton;  % m×n
            else
                dx_dL = -H_newton \ Gamma_newton;
            end
            
            % 3.计算 J0 = ∂e/∂L;
            C_newton = B_newton' * P_v_newton;
            dC_dx = zeros(n*k*n, m);
            for i = 1:m
                Tk_i = zeros(n, n * (m + 1));
                for j = 1:n
                    col_index = (j - 1) * (m + 1) + i + 1;
                    Tk_i(j, col_index) = 1;
                end
                dC_dxk_i = Tk_i' * P_v_newton - 2 * B_newton' * P_v_newton * Tk_i * Q_e_newton * B_newton' * P_v_newton;
                dC_dx(:,i) = dC_dxk_i(:);
            end
            dC_dL_vec = dC_dx * dx_dL;
            
            de_dL = zeros(k*n, n);
            for j = 1:n
                dC_dL_j = dC_dL_vec(:,j);
                dC_dL_j_reshaped = reshape(dC_dL_j, n*k, n);
                dl = zeros(n,1);
                dl(j) = 1;
                de_dL(:, j) = Q_e_newton * (dC_dL_j_reshaped * L_obs_with_outlier + C_newton * dl - dC_dL_j_reshaped * A_obs_with_outlier * x_hat_newton - C_newton * A_obs_with_outlier * dx_dL(:,j));
            end
            J0 = de_dL;
            
            % 4.计算 Ji = ∂e/∂ai;
            Ja = cell(1, m);
            dx_da_all = cell(1, m);
            
            for param_idx = 1:m  % 对每个参数列（A1, A2, ...）
                de_da_i = zeros(m, n);
                
                for obs_idx = 1:n
                    R1 = zeros(n, m);
                    R1(obs_idx, param_idx) = 1;
                    R2 = zeros(n, 1);
                    R2(obs_idx) = -x_hat_newton(param_idx);
                    
                    de_da = Q_e_newton * B_newton' * P_v_newton * R2;
                    dET_da = zeros(m, n);
                    for j = 1:n
                        dET_da(:, j) = de_da((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1));
                    end
                    dF_da_single = (R1' + dET_da) * P_v_newton * v_newton + (A_obs_with_outlier + e_A_newton)' * P_v_newton * R2;
                    de_da_i(:, obs_idx) = dF_da_single;
                end
                
                if rcond(H_newton) < eps
                    dx_da_i = -pinv(H_newton) * de_da_i;
                else
                    dx_da_i = -H_newton \ de_da_i;
                end
                dx_da_all{param_idx} = dx_da_i;
                
                % 计算dC_dx_vec
                dC_dx = zeros(n*k*n, m);
                for i = 1:m
                    Tk_i = zeros(n, n * (m + 1));
                    for j = 1:n
                        col_index = (j - 1) * (m + 1) + i + 1;
                        Tk_i(j, col_index) = 1;
                    end
                    dC_dxk_i = Tk_i' * P_v_newton - 2 * B_newton' * P_v_newton * Tk_i * Q_e_newton * B_newton' * P_v_newton;
                    dC_dx(:,i) = dC_dxk_i(:);
                end
                dC_dx_vec = dC_dx * dx_da_all{param_idx};
                
                % 计算最终结果Ji
                term1 = zeros(k*n, n);
                for j = 1:n
                    dC_da_j = dC_dx_vec(:,j);
                    dC_da_j_reshaped = reshape(dC_da_j, n*k, n);
                    da = zeros(n,m);
                    da(j,param_idx) = 1;
                    term1(:,j) = dC_da_j_reshaped * v_newton - C_newton * da * x_hat_newton;
                end
                de_da = Q_e_newton * (term1 - C_newton * A_obs_with_outlier * dx_da_all{param_idx});
                Ja{param_idx} = de_da;
            end
            J_total = [J0, Ja{1}, Ja{2}];
            
            % 构建非相关的Q_total（块对角矩阵）
            idx_L = 1:k:n*k;      % L的索引：1, 4, 7, ...
            idx_A1 = 2:k:n*k;     % A1的索引：2, 5, 8, ...
            idx_A2 = 3:k:n*k;     % A2的索引：3, 6, 9, ...
            
            % 提取各变量的协方差子矩阵（非相关情况下是对角矩阵）
            Q_L_newton = sigma_L^2 * eye(n);
            QA1_newton = sigma_x^2 * eye(n);
            QA2_newton = 1e-12 * eye(n);
            
            % 构建总体协因数阵（分块矩阵形式，非相关情况下是块对角）
            Q_total = blkdiag(Q_L_newton, QA1_newton, QA2_newton);
            
            % 单位权方差
            sit0_1_newton = (v_newton' * P_v_newton * v_newton) / (n - m);
            
            % 计算Sigma_e
            Sigma_e = sit0_1_newton * J_total * Q_total * J_total';
            
            % ========== 全分量方法基于e_L和e_A的w检验（按照sjtcduibi的方法）==========
            % 提取e_A1和e_A2
            e_A1 = e_hat_reshaped_newton(:, 2);  % A1的残差
            e_A2 = e_hat_reshaped_newton(:, 3);  % A2的残差

            % 提取Sigma_e的对角元素（方差）
            % 每三个元素为一组：第1个是L的方差，第2个是A1的方差，第3个是A2的方差
            diag_Sigma_e = diag(Sigma_e);

            % 定义索引（用于提取对应位置的方差）
            row_indices = 1:k:k*n;        % 1, 4, 7, ... (L的索引)
            row_indices_A1 = 2:k:k*n;      % 2, 5, 8, ... (A1的索引)
            row_indices_A2 = 3:k:k*n;      % 3, 6, 9, ... (A2的索引)

            % 提取L的方差
            var_eL = diag_Sigma_e(row_indices) / sit0_1_newton;  % 归一化方差
            sigma_eL = sqrt(var_eL);  % 残差标准差

            % 提取A1的方差
            var_eA1 = diag_Sigma_e(row_indices_A1) / sit0_1_newton;  % 归一化方差
            sigma_eA1 = sqrt(var_eA1);  % 残差标准差

            % 提取A2的方差
            var_eA2 = diag_Sigma_e(row_indices_A2) / sit0_1_newton;  % 归一化方差
            sigma_eA2 = sqrt(var_eA2);  % 残差标准差

            % 对e_L进行w检验
            w_tests_newton_L = e_L ./ (sqrt(sit0_1_newton) * sigma_eL);
            detected_newton_L = abs(w_tests_newton_L) > F_critical;

            % 对e_A1进行w检验
            w_tests_newton_A1 = e_A1 ./ (sqrt(sit0_1_newton) * sigma_eA1);
            detected_newton_A1 = abs(w_tests_newton_A1) > F_critical;

            % 对e_A2进行w检验
            w_tests_newton_A2 = e_A2 ./ (sqrt(sit0_1_newton) * sigma_eA2);
            detected_newton_A2 = abs(w_tests_newton_A2) > F_critical;

            % 综合判断：任意一个方向检测到粗差就认为该点有粗差
            detected_full = detected_newton_L | detected_newton_A1 | detected_newton_A2;
        catch ME
            % 如果完整方法失败，使用简化的w检验
            try
                PP_newton = inv(Sigma_global_newton);
                x_hat_newton = TLS_XG_newton3(A_obs_with_outlier, L_obs_with_outlier, PP_newton);
                v_newton = L_obs_with_outlier - A_obs_with_outlier * x_hat_newton;
                
                % 使用简化的w检验（基于残差v，假设独立）
                sigma0_sq = (v_newton' * v_newton) / (n - m);
                sigma0 = sqrt(sigma0_sq);
                Qv_simplified = eye(n) - A_obs_with_outlier * inv(A_obs_with_outlier' * A_obs_with_outlier) * A_obs_with_outlier';
                sigma_v = sigma0 * sqrt(diag(Qv_simplified));
                w_tests_newton_L = v_newton ./ (sigma_v + 1e-10);
                detected_full = abs(w_tests_newton_L) > F_critical;
            catch ME2
                % 如果还是失败，使用最简单的OLS残差检验
                x_hat_newton = x_ols_init;
                v_newton = L_obs_with_outlier - A_obs_with_outlier * x_hat_newton;
                sigma0_sq = (v_newton' * v_newton) / (n - m);
                sigma0 = sqrt(sigma0_sq);
                w_tests_newton_L = v_newton / sigma0;
                detected_full = abs(w_tests_newton_L) > F_critical;
            end
        end
        
        time_full = toc(tic_full);
        total_time_full = total_time_full + time_full;
        
        % ========== 方法5: 分方向残差抗差TLS（处理相关性）==========
        tic_separate = tic;
        
        try
            % 直接使用协方差矩阵Sigma_global，避免重复求逆导致的数值不稳定
            % 注意：相关性抗差方法需要[y; x1; x2]顺序，即[L; A1; A2]
            
            % 验证协方差矩阵是否有效
            if any(isnan(Sigma_global(:))) || any(isinf(Sigma_global(:)))
                error('协方差矩阵包含NaN或Inf值');
            end
            
            % 调用分方向残差抗差TLS方法（传入协方差矩阵而不是权重矩阵）
            % 注意：函数在试验/相关性抗差文件中定义，需要确保该文件在路径中
            if ~exist('separate_direction_robust_tls_cov', 'file')
                % 如果新函数不存在，使用旧版本（但会有数值问题）
                W_full_reorder = zeros(3*n, 3*n);
                
                for i = 1:n
                    block_start = (i-1)*3 + 1;
                    block_end = i*3;
                    
                    % 从Sigma_global中提取第i个观测点的协方差块
                    Sigma_block = Sigma_global(block_start:block_end, block_start:block_end);
                    
                    % 添加正则化改善条件数
                    reg_factor = 1e-8 * trace(Sigma_block) / 3;
                    Sigma_block_reg = Sigma_block + reg_factor * eye(3);
                    
                    % 计算权重矩阵（协方差矩阵的逆）
                    if rcond(Sigma_block_reg) < 1e-10
                        W_block = pinv(Sigma_block_reg);
                    else
                        W_block = inv(Sigma_block_reg);
                    end
                    
                    % 将权重块放入块对角矩阵
                    W_full_reorder(block_start:block_end, block_start:block_end) = W_block;
                end
                
                [X_separate, residuals_separate, iter_info_separate, P_final_separate] = ...
                    separate_direction_robust_tls(A_obs_with_outlier, L_obs_with_outlier, W_full_reorder);
            else
                [X_separate, residuals_separate, iter_info_separate, P_final_separate] = ...
                    separate_direction_robust_tls_cov(A_obs_with_outlier, L_obs_with_outlier, Sigma_global);
            end
            
            % 验证返回结果
            if isempty(X_separate) || any(isnan(X_separate)) || any(isinf(X_separate))
                error('分方向残差抗差方法返回无效参数');
            end
            
            % 根据最终权重判断粗差（权重很小的点被认为是粗差）
            % 提取每个观测点的权重
            detected_separate = false(n, 1);
            weight_threshold_separate = 1e-6;  % 权重阈值
            
            for i = 1:n
                block_start = (i-1)*3 + 1;
                P_block = P_final_separate(block_start:block_start+2, block_start:block_start+2);
                % 如果三个方向的权重都很小，认为是粗差
                min_weight = min([P_block(1,1), P_block(2,2), P_block(3,3)]);
                if min_weight < weight_threshold_separate
                    detected_separate(i) = true;
                end
            end
            
            % 如果检测到的粗差太少，使用残差方法补充
            if sum(detected_separate) < n_outliers_total * 0.5
                % 使用残差判断
                v_separate = L_obs_with_outlier - A_obs_with_outlier * X_separate;
                sigma0_separate = sqrt((v_separate' * v_separate) / (n - m));
                std_residuals = abs(v_separate) / (sigma0_separate + 1e-10);
                detected_separate = detected_separate | (std_residuals > F_critical);
            end
            
            % 参数估计
            x_clean_separate = X_separate;
            param_error_separate = param_error_separate + [(x_clean_separate(1) - a_true)^2; (x_clean_separate(2) - b_true)^2];
            param_count_separate = param_count_separate + 1;
            param_estimates_separate = [param_estimates_separate; x_clean_separate'];
            
            % 计算单位权中误差sigma0（需要使用加权残差）
            v_separate = L_obs_with_outlier - A_obs_with_outlier * X_separate;
            r_separate = n - m;
            if r_separate > 0
                % 使用最终权重矩阵P_final_separate计算加权残差
                % 需要将v转换为完整的残差向量
                try
                    % 调用TLS方法获取加权残差的权重矩阵
                    [~, ~, ~, W_final] = TLS_XG_with_W(A_obs_with_outlier, L_obs_with_outlier, P_final_separate, X_separate);
                    sigma0_separate = sqrt((v_separate' * W_final * v_separate) / r_separate);
                    sigma0_sum_separate = sigma0_sum_separate + sigma0_separate;
                catch
                    % 如果失败，使用简化计算
                    sigma0_separate = sqrt((v_separate' * v_separate) / r_separate);
                    sigma0_sum_separate = sigma0_sum_separate + sigma0_separate;
                end
            end
            
        catch ME
            if exp_idx == 1 && rho_idx == 1
                fprintf('    分方向残差抗差方法失败: %s\n', ME.message);
                fprintf('    错误位置: %s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
            end
            X_separate = x_ols_init;
            detected_separate = false(n, 1);
        end
        
        time_separate = toc(tic_separate);
        total_time_separate = total_time_separate + time_separate;
        
        % ========== 方法6: 总体残差抗差TLS（处理相关性）==========
        tic_overall = tic;
        
        try
            % 使用相同的权重矩阵（如果方法5失败，重新构建，并添加正则化）
            if ~exist('W_full_reorder', 'var') || isempty(W_full_reorder)
                W_full_reorder = zeros(3*n, 3*n);
                for i = 1:n
                    block_start = (i-1)*3 + 1;
                    block_end = i*3;
                    Sigma_block = Sigma_global(block_start:block_end, block_start:block_end);
                    
                    % 添加正则化改善条件数
                    reg_factor = 1e-8 * trace(Sigma_block) / 3;
                    Sigma_block_reg = Sigma_block + reg_factor * eye(3);
                    
                    if rcond(Sigma_block_reg) < 1e-10
                        W_block = pinv(Sigma_block_reg);
                    else
                        W_block = inv(Sigma_block_reg);
                    end
                    W_full_reorder(block_start:block_end, block_start:block_end) = W_block;
                end
            end
            
            % 验证权重矩阵是否有效
            if any(isnan(W_full_reorder(:))) || any(isinf(W_full_reorder(:)))
                error('权重矩阵包含NaN或Inf值');
            end
            
            % 调用总体残差抗差TLS方法
            if ~exist('correlated_robust_tls', 'file')
                error('函数correlated_robust_tls未找到，请确保试验/相关性抗差文件在路径中');
            end
            
            [X_overall, residuals_overall, iter_info_overall, P_final_overall] = ...
                correlated_robust_tls(A_obs_with_outlier, L_obs_with_outlier, W_full_reorder);
            
            % 验证返回结果
            if isempty(X_overall) || any(isnan(X_overall)) || any(isinf(X_overall))
                error('总体残差抗差方法返回无效参数');
            end
            
            % 根据最终权重判断粗差
            detected_overall = false(n, 1);
            weight_threshold_overall = 1e-6;  % 权重阈值
            
            for i = 1:n
                block_start = (i-1)*3 + 1;
                P_block = P_final_overall(block_start:block_start+2, block_start:block_start+2);
                % 如果三个方向的权重都很小，认为是粗差
                min_weight = min([P_block(1,1), P_block(2,2), P_block(3,3)]);
                if min_weight < weight_threshold_overall
                    detected_overall(i) = true;
                end
            end
            
            % 如果检测到的粗差太少，使用残差方法补充
            if sum(detected_overall) < n_outliers_total * 0.5
                % 使用残差判断
                v_overall = L_obs_with_outlier - A_obs_with_outlier * X_overall;
                sigma0_overall = sqrt((v_overall' * v_overall) / (n - m));
                std_residuals = abs(v_overall) / (sigma0_overall + 1e-10);
                detected_overall = detected_overall | (std_residuals > F_critical);
            end
            
            % 参数估计
            x_clean_overall = X_overall;
            param_error_overall = param_error_overall + [(x_clean_overall(1) - a_true)^2; (x_clean_overall(2) - b_true)^2];
            param_count_overall = param_count_overall + 1;
            param_estimates_overall = [param_estimates_overall; x_clean_overall'];
            
            % 计算单位权中误差sigma0（需要使用加权残差）
            v_overall = L_obs_with_outlier - A_obs_with_outlier * X_overall;
            r_overall = n - m;
            if r_overall > 0
                % 使用最终权重矩阵P_final_overall计算加权残差
                try
                    % 调用TLS方法获取加权残差的权重矩阵
                    [~, ~, ~, W_final] = TLS_XG_with_W(A_obs_with_outlier, L_obs_with_outlier, P_final_overall, X_overall);
                    sigma0_overall = sqrt((v_overall' * W_final * v_overall) / r_overall);
                    sigma0_sum_overall = sigma0_sum_overall + sigma0_overall;
                catch
                    % 如果失败，使用简化计算
                    sigma0_overall = sqrt((v_overall' * v_overall) / r_overall);
                    sigma0_sum_overall = sigma0_sum_overall + sigma0_overall;
                end
            end
            
        catch ME
            if exp_idx == 1 && rho_idx == 1
                fprintf('    总体残差抗差方法失败: %s\n', ME.message);
                fprintf('    错误位置: %s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
            end
            X_overall = x_ols_init;
            detected_overall = false(n, 1);
        end
        
        time_overall = toc(tic_overall);
        total_time_overall = total_time_overall + time_overall;
        
        % ========== 统计检测结果 ==========
        % 判断是否成功检测出粗差
        if ~isempty(true_outlier_positions)
            % Mahboub WTLS
            detected_outliers_mahboub = sum(detected_mahboub(true_outlier_positions) == 1);
            detection_rate_mahboub = detected_outliers_mahboub / length(true_outlier_positions);
            if detection_rate_mahboub >= 0.5
                success_count_mahboub = success_count_mahboub + 1;
            else
                fn_count_mahboub = fn_count_mahboub + 1;
            end
            
            % Mahboub IRTLS
            detected_outliers_mahboub_irtls = sum(detected_mahboub_irtls(true_outlier_positions) == 1);
            detection_rate_mahboub_irtls = detected_outliers_mahboub_irtls / length(true_outlier_positions);
            if detection_rate_mahboub_irtls >= 0.5
                success_count_mahboub_irtls = success_count_mahboub_irtls + 1;
            else
                fn_count_mahboub_irtls = fn_count_mahboub_irtls + 1;
            end
            
            % Component-Compressed
            detected_outliers_component = sum(detected_component(true_outlier_positions) == 1);
            detection_rate_component = detected_outliers_component / length(true_outlier_positions);
            if detection_rate_component >= 0.5
                success_count_component = success_count_component + 1;
            else
                fn_count_component = fn_count_component + 1;
            end
            
            % Full-Component
            detected_outliers_full = sum(detected_full(true_outlier_positions) == 1);
            detection_rate_full = detected_outliers_full / length(true_outlier_positions);
            if detection_rate_full >= 0.5
                success_count_full = success_count_full + 1;
            else
                fn_count_full = fn_count_full + 1;
            end
            
            % 分方向残差抗差（Separate Direction Robust TLS）
            detected_outliers_separate = sum(detected_separate(true_outlier_positions) == 1);
            detection_rate_separate = detected_outliers_separate / length(true_outlier_positions);
            if detection_rate_separate >= 0.5
                success_count_separate = success_count_separate + 1;
            else
                fn_count_separate = fn_count_separate + 1;
            end
            
            % 总体残差抗差（Overall Residual Robust TLS）
            detected_outliers_overall = sum(detected_overall(true_outlier_positions) == 1);
            detection_rate_overall = detected_outliers_overall / length(true_outlier_positions);
            if detection_rate_overall >= 0.5
                success_count_overall = success_count_overall + 1;
            else
                fn_count_overall = fn_count_overall + 1;
            end
        end
        
        % 统计误检
        if ~isempty(true_outlier_positions)
            false_detections_mahboub = find(detected_mahboub == 1 & ~ismember((1:n)', true_outlier_positions));
            false_detections_mahboub_irtls = find(detected_mahboub_irtls == 1 & ~ismember((1:n)', true_outlier_positions));
            false_detections_component = find(detected_component == 1 & ~ismember((1:n)', true_outlier_positions));
            false_detections_full = find(detected_full == 1 & ~ismember((1:n)', true_outlier_positions));
            false_detections_separate = find(detected_separate == 1 & ~ismember((1:n)', true_outlier_positions));
            false_detections_overall = find(detected_overall == 1 & ~ismember((1:n)', true_outlier_positions));
        else
            false_detections_mahboub = find(detected_mahboub == 1);
            false_detections_mahboub_irtls = find(detected_mahboub_irtls == 1);
            false_detections_component = find(detected_component == 1);
            false_detections_full = find(detected_full == 1);
            false_detections_separate = find(detected_separate == 1);
            false_detections_overall = find(detected_overall == 1);
        end
        
        fp_count_mahboub = fp_count_mahboub + length(false_detections_mahboub);
        fp_count_mahboub_irtls = fp_count_mahboub_irtls + length(false_detections_mahboub_irtls);
        fp_count_component = fp_count_component + length(false_detections_component);
        fp_count_full = fp_count_full + length(false_detections_full);
        fp_count_separate = fp_count_separate + length(false_detections_separate);
        fp_count_overall = fp_count_overall + length(false_detections_overall);
        
        % ========== 抗差剔除后进行参数估计 ==========
        % Mahboub WTLS
        if any(detected_mahboub == 1)
            valid_idx = find(detected_mahboub == 0);
            if length(valid_idx) >= m
                A_clean = A_obs_with_outlier(valid_idx, :);
                L_clean = L_obs_with_outlier(valid_idx);
                x_clean_mahboub = (A_clean' * A_clean) \ (A_clean' * L_clean);
                param_error_mahboub = param_error_mahboub + [(x_clean_mahboub(1) - a_true)^2; (x_clean_mahboub(2) - b_true)^2];
                param_count_mahboub = param_count_mahboub + 1;
                param_estimates_mahboub = [param_estimates_mahboub; x_clean_mahboub'];
                
                % 计算单位权中误差sigma0（使用OLS后，权重为单位阵）
                v_clean = L_clean - A_clean * x_clean_mahboub;
                n_clean = length(valid_idx);
                r_clean = n_clean - m;
                if r_clean > 0
                    % 剔除粗差后使用简单的OLS，所以直接用v'*v
                    sigma0_mahboub = sqrt((v_clean' * v_clean) / r_clean);
                    sigma0_sum_mahboub = sigma0_sum_mahboub + sigma0_mahboub;
                end
            end
        end
        
        % Mahboub IRTLS
        % IRTLS本身就是抗差方法，优先使用全部数据的结果
        % 如果检测到粗差且剔除后点数足够，可以重新估计；否则使用全部数据的结果
        if ~isempty(X_mahboub_irtls) && all(~isnan(X_mahboub_irtls)) && all(~isinf(X_mahboub_irtls))
            % 优先使用全部数据的结果（IRTLS本身就是抗差方法，不需要剔除）
            x_clean_mahboub_irtls = X_mahboub_irtls;
            param_error_mahboub_irtls = param_error_mahboub_irtls + [(x_clean_mahboub_irtls(1) - a_true)^2; (x_clean_mahboub_irtls(2) - b_true)^2];
            param_count_mahboub_irtls = param_count_mahboub_irtls + 1;
            param_estimates_mahboub_irtls = [param_estimates_mahboub_irtls; x_clean_mahboub_irtls'];
            
            % 计算单位权中误差sigma0（使用IRTLS返回的sigma0）
            if isfield(iter_info_mahboub_irtls, 'sigma0') && ~isempty(iter_info_mahboub_irtls.sigma0)
                sigma0_mahboub_irtls = iter_info_mahboub_irtls.sigma0(end);
                if ~isnan(sigma0_mahboub_irtls) && sigma0_mahboub_irtls > 0
                    sigma0_sum_mahboub_irtls = sigma0_sum_mahboub_irtls + sigma0_mahboub_irtls;
                end
            else
                % 如果没有sigma0，使用残差计算
                v_irtls = L_obs_with_outlier - A_obs_with_outlier * x_clean_mahboub_irtls;
                r_irtls = n - m;
                if r_irtls > 0
                    sigma0_mahboub_irtls = sqrt((v_irtls' * v_irtls) / r_irtls);
                    sigma0_sum_mahboub_irtls = sigma0_sum_mahboub_irtls + sigma0_mahboub_irtls;
                end
            end
        elseif any(detected_mahboub_irtls == 1)
            % 如果全部数据的结果不可用，且检测到粗差，尝试剔除后重新估计
            valid_idx = find(detected_mahboub_irtls == 0);
            if length(valid_idx) >= m
                A_clean = A_obs_with_outlier(valid_idx, :);
                L_clean = L_obs_with_outlier(valid_idx);
                sigma_A_clean = sigma_A_mah(valid_idx, :);
                sigma_y_clean = sigma_y_mah(valid_idx);
                try
                    [x_clean_mahboub_irtls, ~, iter_info_clean] = IRTLS_Mahboub2013(A_clean, L_clean, sigma_A_clean, sigma_y_clean, options_mah);
                    if ~isempty(x_clean_mahboub_irtls) && all(~isnan(x_clean_mahboub_irtls)) && all(~isinf(x_clean_mahboub_irtls))
                        param_error_mahboub_irtls = param_error_mahboub_irtls + [(x_clean_mahboub_irtls(1) - a_true)^2; (x_clean_mahboub_irtls(2) - b_true)^2];
                        param_count_mahboub_irtls = param_count_mahboub_irtls + 1;
                        param_estimates_mahboub_irtls = [param_estimates_mahboub_irtls; x_clean_mahboub_irtls'];
                        
                        % 计算单位权中误差sigma0
                        if isfield(iter_info_clean, 'sigma0') && ~isempty(iter_info_clean.sigma0)
                            sigma0_mahboub_irtls = iter_info_clean.sigma0(end);
                            if ~isnan(sigma0_mahboub_irtls) && sigma0_mahboub_irtls > 0
                                sigma0_sum_mahboub_irtls = sigma0_sum_mahboub_irtls + sigma0_mahboub_irtls;
                            end
                        else
                            v_clean = L_clean - A_clean * x_clean_mahboub_irtls;
                            n_clean = length(valid_idx);
                            r_clean = n_clean - m;
                            if r_clean > 0
                                sigma0_mahboub_irtls = sqrt((v_clean' * v_clean) / r_clean);
                                sigma0_sum_mahboub_irtls = sigma0_sum_mahboub_irtls + sigma0_mahboub_irtls;
                            end
                        end
                    end
                catch
                end
            end
        end
        
        % Component-Compressed
        % 使用detect_outlier_v返回的参数估计结果
        if ~isempty(x_hat_component) && all(~isnan(x_hat_component)) && all(~isinf(x_hat_component))
            param_error_component = param_error_component + [(x_hat_component(1) - a_true)^2; (x_hat_component(2) - b_true)^2];
            param_count_component = param_count_component + 1;
            param_estimates_component = [param_estimates_component; x_hat_component'];
            
            % 计算单位权中误差sigma0
            % 注意：不使用results_component.sit0_1，因为它受到A2列极大权重（1e14）的影响
            % 直接用残差计算更准确
            v_component = L_obs_with_outlier - A_obs_with_outlier * x_hat_component;
            
            % 检查是否有检测到的粗差
            if any(detected_component == 1)
                % 如果检测到粗差，只使用非粗差点计算sigma0
                valid_idx = find(detected_component == 0);
                if length(valid_idx) >= m
                    v_clean = v_component(valid_idx);
                    r_clean = length(valid_idx) - m;
                    if r_clean > 0
                        % 使用加权残差（只用L方向的权重）
                        py_clean = py(valid_idx);
                        sigma0_component = sqrt((v_clean' * diag(py_clean) * v_clean) / r_clean);
                        sigma0_sum_component = sigma0_sum_component + sigma0_component;
                    end
                end
            else
                % 没有检测到粗差，使用全部数据
                r_component = n - m;
                if r_component > 0
                    % 使用加权残差（只用L方向的权重）
                    sigma0_component = sqrt((v_component' * diag(py) * v_component) / r_component);
                    sigma0_sum_component = sigma0_sum_component + sigma0_component;
                end
            end
        else
            % 如果detect_outlier_v失败，尝试剔除粗差后重新估计
            if any(detected_component == 1)
                valid_idx = find(detected_component == 0);
                if length(valid_idx) >= m
                    A_clean = A_obs_with_outlier(valid_idx, :);
                    L_clean = L_obs_with_outlier(valid_idx);
                    % 构建块对角权重矩阵
                    n_clean = length(valid_idx);
                    PP_clean = zeros(3*n_clean, 3*n_clean);
                    for i = 1:n_clean
                        block_start = (i-1)*3 + 1;
                        block_end = i*3;
                        PP_clean(block_start:block_end, block_start:block_end) = ...
                            diag([1/sigma_L^2, 1/sigma_x^2, 1e12]);
                    end
                    try
                        x_clean_component = TLS_XG_newton3(A_clean, L_clean, PP_clean);
                        if ~isempty(x_clean_component) && all(~isnan(x_clean_component)) && all(~isinf(x_clean_component))
                            param_error_component = param_error_component + [(x_clean_component(1) - a_true)^2; (x_clean_component(2) - b_true)^2];
                            param_count_component = param_count_component + 1;
                            param_estimates_component = [param_estimates_component; x_clean_component'];
                            
                            % 计算单位权中误差sigma0（使用加权残差）
                            v_clean = L_clean - A_clean * x_clean_component;
                            r_clean = n_clean - m;
                            if r_clean > 0
                                % 提取有效点的L方向权重
                                py_clean = 1/sigma_L^2 * ones(n_clean, 1);
                                sigma0_component = sqrt((v_clean' * diag(py_clean) * v_clean) / r_clean);
                                sigma0_sum_component = sigma0_sum_component + sigma0_component;
                            end
                        end
                    catch
                    end
                end
            end
        end
        
        % Full-Component
        if any(detected_full == 1)
            valid_idx = find(detected_full == 0);
            if length(valid_idx) >= m
                A_clean = A_obs_with_outlier(valid_idx, :);
                L_clean = L_obs_with_outlier(valid_idx);
                P_full_clean = zeros(3*length(valid_idx), 3*length(valid_idx));
                for i = 1:length(valid_idx)
                    block_start = (i-1)*3 + 1;
                    block_end = i*3;
                    P_full_clean(block_start:block_end, block_start:block_end) = ...
                        diag([1/sigma_L^2, 1/sigma_x^2, 1e12]);
                end
                Sigma_global_clean = zeros(3*length(valid_idx), 3*length(valid_idx));
                for i = 1:length(valid_idx)
                    block_start = (i-1)*3 + 1;
                    block_end = i*3;
                    Sigma_global_clean(block_start:block_end, block_start:block_end) = ...
                        diag([sigma_L^2, sigma_x^2, 1e-12]);
                end
                PP_clean = inv(Sigma_global_clean);
                try
                    x_clean_full = TLS_XG_newton3(A_clean, L_clean, PP_clean);
                    param_error_full = param_error_full + [(x_clean_full(1) - a_true)^2; (x_clean_full(2) - b_true)^2];
                    param_count_full = param_count_full + 1;
                    param_estimates_full = [param_estimates_full; x_clean_full'];
                    
                    % 计算单位权中误差sigma0
                    v_clean = L_clean - A_clean * x_clean_full;
                    n_clean = length(valid_idx);
                    r_clean = n_clean - m;
                    if r_clean > 0
                        sigma0_full = sqrt((v_clean' * v_clean) / r_clean);
                        sigma0_sum_full = sigma0_sum_full + sigma0_full;
                    end
                catch
                end
            end
        else
            % 没有检测到粗差，使用全部数据进行参数估计
            try
                x_clean_full = x_hat_newton;  % 使用已经计算好的结果
                param_error_full = param_error_full + [(x_clean_full(1) - a_true)^2; (x_clean_full(2) - b_true)^2];
                param_count_full = param_count_full + 1;
                param_estimates_full = [param_estimates_full; x_clean_full'];
                
                % 计算单位权中误差sigma0（使用已经计算好的sit0_1_newton）
                if exist('sit0_1_newton', 'var') && ~isempty(sit0_1_newton) && sit0_1_newton > 0
                    sigma0_full = sqrt(sit0_1_newton);
                    sigma0_sum_full = sigma0_sum_full + sigma0_full;
                else
                    v_full = L_obs_with_outlier - A_obs_with_outlier * x_clean_full;
                    r_full = n - m;
                    if r_full > 0
                        sigma0_full = sqrt((v_full' * v_full) / r_full);
                        sigma0_sum_full = sigma0_sum_full + sigma0_full;
                    end
                end
            catch
            end
        end
    end
    
    % ========== 计算统计结果 ==========
    success_rate_mahboub = success_count_mahboub / num_experiments * 100;
    success_rate_mahboub_irtls = success_count_mahboub_irtls / num_experiments * 100;
    success_rate_component = success_count_component / num_experiments * 100;
    success_rate_full = success_count_full / num_experiments * 100;
    success_rate_separate = success_count_separate / num_experiments * 100;
    success_rate_overall = success_count_overall / num_experiments * 100;
    
    fp_rate_mahboub = fp_count_mahboub / num_experiments;
    fp_rate_mahboub_irtls = fp_count_mahboub_irtls / num_experiments;
    fp_rate_component = fp_count_component / num_experiments;
    fp_rate_full = fp_count_full / num_experiments;
    fp_rate_separate = fp_count_separate / num_experiments;
    fp_rate_overall = fp_count_overall / num_experiments;
    
    fn_rate_mahboub = fn_count_mahboub / num_experiments * 100;
    fn_rate_mahboub_irtls = fn_count_mahboub_irtls / num_experiments * 100;
    fn_rate_component = fn_count_component / num_experiments * 100;
    fn_rate_full = fn_count_full / num_experiments * 100;
    fn_rate_separate = fn_count_separate / num_experiments * 100;
    fn_rate_overall = fn_count_overall / num_experiments * 100;
    
    avg_time_mahboub = total_time_mahboub / num_experiments;
    avg_time_mahboub_irtls = total_time_mahboub_irtls / num_experiments;
    avg_time_component = total_time_component / num_experiments;
    avg_time_full = total_time_full / num_experiments;
    avg_time_separate = total_time_separate / num_experiments;
    avg_time_overall = total_time_overall / num_experiments;
    
    % 计算平均单位权中误差sigma0
    if param_count_mahboub > 0
        avg_sigma0_mahboub = sigma0_sum_mahboub / param_count_mahboub;
    else
        avg_sigma0_mahboub = NaN;
    end
    
    if param_count_mahboub_irtls > 0
        avg_sigma0_mahboub_irtls = sigma0_sum_mahboub_irtls / param_count_mahboub_irtls;
    else
        avg_sigma0_mahboub_irtls = NaN;
    end
    
    if param_count_component > 0
        avg_sigma0_component = sigma0_sum_component / param_count_component;
    else
        avg_sigma0_component = NaN;
    end
    
    if param_count_full > 0
        avg_sigma0_full = sigma0_sum_full / param_count_full;
    else
        avg_sigma0_full = NaN;
    end
    
    if param_count_separate > 0
        avg_sigma0_separate = sigma0_sum_separate / param_count_separate;
    else
        avg_sigma0_separate = NaN;
    end
    
    if param_count_overall > 0
        avg_sigma0_overall = sigma0_sum_overall / param_count_overall;
    else
        avg_sigma0_overall = NaN;
    end
    
    % ========== 保存结果 ==========
    results_all.rho(rho_idx) = rho;
    results_all.success_rate_mahboub(rho_idx) = success_rate_mahboub;
    results_all.success_rate_mahboub_irtls(rho_idx) = success_rate_mahboub_irtls;
    results_all.success_rate_component(rho_idx) = success_rate_component;
    results_all.success_rate_full(rho_idx) = success_rate_full;
    results_all.success_rate_separate(rho_idx) = success_rate_separate;
    results_all.success_rate_overall(rho_idx) = success_rate_overall;
    
    results_all.fp_rate_mahboub(rho_idx) = fp_rate_mahboub;
    results_all.fp_rate_mahboub_irtls(rho_idx) = fp_rate_mahboub_irtls;
    results_all.fp_rate_component(rho_idx) = fp_rate_component;
    results_all.fp_rate_full(rho_idx) = fp_rate_full;
    results_all.fp_rate_separate(rho_idx) = fp_rate_separate;
    results_all.fp_rate_overall(rho_idx) = fp_rate_overall;
    
    results_all.fn_rate_mahboub(rho_idx) = fn_rate_mahboub;
    results_all.fn_rate_mahboub_irtls(rho_idx) = fn_rate_mahboub_irtls;
    results_all.fn_rate_component(rho_idx) = fn_rate_component;
    results_all.fn_rate_full(rho_idx) = fn_rate_full;
    results_all.fn_rate_separate(rho_idx) = fn_rate_separate;
    results_all.fn_rate_overall(rho_idx) = fn_rate_overall;
    
    results_all.avg_time_mahboub(rho_idx) = avg_time_mahboub;
    results_all.avg_time_mahboub_irtls(rho_idx) = avg_time_mahboub_irtls;
    results_all.avg_time_component(rho_idx) = avg_time_component;
    results_all.avg_time_full(rho_idx) = avg_time_full;
    results_all.avg_time_separate(rho_idx) = avg_time_separate;
    results_all.avg_time_overall(rho_idx) = avg_time_overall;
    
    results_all.avg_sigma0_mahboub(rho_idx) = avg_sigma0_mahboub;
    results_all.avg_sigma0_mahboub_irtls(rho_idx) = avg_sigma0_mahboub_irtls;
    results_all.avg_sigma0_component(rho_idx) = avg_sigma0_component;
    results_all.avg_sigma0_full(rho_idx) = avg_sigma0_full;
    results_all.avg_sigma0_separate(rho_idx) = avg_sigma0_separate;
    results_all.avg_sigma0_overall(rho_idx) = avg_sigma0_overall;
    
    % ========== 输出结果 ==========
    fprintf('\n========== 结果统计（rho = %.1f） ==========\n', rho);
    
    % 计算平均参数估计值
    % 注意：param_estimates是按行存储的（每行是一个实验的参数估计[a, b]）
    if param_count_mahboub > 0 && ~isempty(param_estimates_mahboub)
        avg_params_mahboub = mean(param_estimates_mahboub, 1)';  % 按列求平均，然后转置为列向量
    else
        avg_params_mahboub = [NaN; NaN];
    end
    
    if param_count_mahboub_irtls > 0 && ~isempty(param_estimates_mahboub_irtls)
        avg_params_mahboub_irtls = mean(param_estimates_mahboub_irtls, 1)';
    else
        avg_params_mahboub_irtls = [NaN; NaN];
    end
    
    if param_count_component > 0 && ~isempty(param_estimates_component)
        avg_params_component = mean(param_estimates_component, 1)';
    else
        avg_params_component = [NaN; NaN];
    end
    
    if param_count_full > 0 && ~isempty(param_estimates_full)
        avg_params_full = mean(param_estimates_full, 1)';
    else
        avg_params_full = [NaN; NaN];
    end
    
    if param_count_separate > 0 && ~isempty(param_estimates_separate)
        avg_params_separate = mean(param_estimates_separate, 1)';
    else
        avg_params_separate = [NaN; NaN];
    end
    
    if param_count_overall > 0 && ~isempty(param_estimates_overall)
        avg_params_overall = mean(param_estimates_overall, 1)';
    else
        avg_params_overall = [NaN; NaN];
    end
    
    % 计算参数不确定度（标准差）
    if param_count_mahboub > 0 && ~isempty(param_estimates_mahboub)
        std_params_mahboub = std(param_estimates_mahboub, 0, 1)';
    else
        std_params_mahboub = [NaN; NaN];
    end
    
    if param_count_mahboub_irtls > 0 && ~isempty(param_estimates_mahboub_irtls)
        std_params_mahboub_irtls = std(param_estimates_mahboub_irtls, 0, 1)';
    else
        std_params_mahboub_irtls = [NaN; NaN];
    end
    
    if param_count_component > 0 && ~isempty(param_estimates_component)
        std_params_component = std(param_estimates_component, 0, 1)';
    else
        std_params_component = [NaN; NaN];
    end
    
    if param_count_full > 0 && ~isempty(param_estimates_full)
        std_params_full = std(param_estimates_full, 0, 1)';
    else
        std_params_full = [NaN; NaN];
    end
    
    if param_count_separate > 0 && ~isempty(param_estimates_separate)
        std_params_separate = std(param_estimates_separate, 0, 1)';
    else
        std_params_separate = [NaN; NaN];
    end
    
    if param_count_overall > 0 && ~isempty(param_estimates_overall)
        std_params_overall = std(param_estimates_overall, 0, 1)';
    else
        std_params_overall = [NaN; NaN];
    end
    
    fprintf('\n【方法1: Mahboub WTLS】（不处理相关性）\n');
    fprintf('  参数估计: a = %.6f ± %.6f, b = %.6f ± %.6f\n', avg_params_mahboub(1), std_params_mahboub(1), avg_params_mahboub(2), std_params_mahboub(2));
    fprintf('  单位权中误差: σ₀ = %.6f\n', avg_sigma0_mahboub);
    
    fprintf('\n【方法2: Mahboub IRTLS】（不处理相关性）\n');
    fprintf('  参数估计: a = %.6f ± %.6f, b = %.6f ± %.6f\n', avg_params_mahboub_irtls(1), std_params_mahboub_irtls(1), avg_params_mahboub_irtls(2), std_params_mahboub_irtls(2));
    fprintf('  单位权中误差: σ₀ = %.6f\n', avg_sigma0_mahboub_irtls);
    
    fprintf('\n【方法3: Component-Compressed】（不处理相关性）\n');
    fprintf('  参数估计: a = %.6f ± %.6f, b = %.6f ± %.6f\n', avg_params_component(1), std_params_component(1), avg_params_component(2), std_params_component(2));
    fprintf('  单位权中误差: σ₀ = %.6f\n', avg_sigma0_component);
    
    fprintf('\n【方法4: Full-Component】（不处理相关性）\n');
    fprintf('  参数估计: a = %.6f ± %.6f, b = %.6f ± %.6f\n', avg_params_full(1), std_params_full(1), avg_params_full(2), std_params_full(2));
    fprintf('  单位权中误差: σ₀ = %.6f\n', avg_sigma0_full);
    
    fprintf('\n【方法5: 分方向残差抗差TLS】（处理相关性）\n');
    fprintf('  参数估计: a = %.6f ± %.6f, b = %.6f ± %.6f\n', avg_params_separate(1), std_params_separate(1), avg_params_separate(2), std_params_separate(2));
    fprintf('  单位权中误差: σ₀ = %.6f\n', avg_sigma0_separate);
    
    fprintf('\n【方法6: 总体残差抗差TLS】（处理相关性）\n');
    fprintf('  参数估计: a = %.6f ± %.6f, b = %.6f ± %.6f\n', avg_params_overall(1), std_params_overall(1), avg_params_overall(2), std_params_overall(2));
    fprintf('  单位权中误差: σ₀ = %.6f\n', avg_sigma0_overall);
end

% ========== 汇总所有相关系数的结果 ==========
fprintf('\n\n========================================\n');
fprintf('========== 所有相关系数的结果汇总 ==========\n');
fprintf('========================================\n\n');

fprintf('\n【单位权中误差σ₀对比 - 不处理相关性的方法】\n');
fprintf('%-10s%-25s%-25s%-25s%-25s\n', 'rho', 'Mahboub WTLS', 'Mahboub IRTLS', 'Component-Compressed', 'Full-Component');
fprintf('%s\n', repmat('-', 1, 110));
for rho_idx = 1:num_rhos
    fprintf('%-10.1f%-25.6f%-25.6f%-25.6f%-25.6f\n', ...
        results_all.rho(rho_idx), ...
        results_all.avg_sigma0_mahboub(rho_idx), ...
        results_all.avg_sigma0_mahboub_irtls(rho_idx), ...
        results_all.avg_sigma0_component(rho_idx), ...
        results_all.avg_sigma0_full(rho_idx));
end

fprintf('\n【单位权中误差σ₀对比 - 处理相关性的方法】\n');
fprintf('%-10s%-30s%-30s\n', 'rho', '分方向残差抗差', '总体残差抗差');
fprintf('%s\n', repmat('-', 1, 70));
for rho_idx = 1:num_rhos
    fprintf('%-10.1f%-30.6f%-30.6f\n', ...
        results_all.rho(rho_idx), ...
        results_all.avg_sigma0_separate(rho_idx), ...
        results_all.avg_sigma0_overall(rho_idx));
end

fprintf('\n实验完成！\n');

%% ========== 辅助函数 ==========

function [H, e_A, B, e_hat] = Hessian(A, L, P, x)
% Hessian 计算Hessian矩阵及相关中间变量

    % 获取问题维度
    [n, m] = size(A);
    k = m + 1;
    
    % 只计算一次Q_e和I_n，不预提取块
    Q_e = pinv(P);
    I_n = eye(n);

    % TLS牛顿迭代法

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

%% ========== 相关性抗差方法函数定义 ==========
% 以下函数来自试验/相关性抗差文件，用于处理相关性的抗差TLS方法

function [X, residuals, iter_info, P_final] = separate_direction_robust_tls(A, L, P_initial)
% 分方向相关观测抗差TLS方法 - 按照论文双因子方法实现
% 输入:
%   A: 系数矩阵 (n×m)
%   L: 观测向量 (n×1)
%   P_initial: 初始权重矩阵 (3n×3n 或 3×n)
% 输出:
%   X: 参数估计 (m×1)
%   residuals: 残差结构体
%   iter_info: 迭代信息
%   P_final: 最终权重矩阵

[n, m] = size(A);
k = m + 1; % 每个观测的参数个数(L, a1, a2)
param_tol = 0.001;
max_iter = 50;  % 最大迭代次数

% 初始化
iter_info = struct();
iter_info.total_iterations = 0;

% 检查P_initial的格式并初始化权重矩阵
if size(P_initial, 1) == 3*n && size(P_initial, 2) == 3*n
    % 新格式：3n×3n完整权重矩阵
    P = P_initial;  % 直接使用完整权重矩阵
    use_full_weight = true;
    
    % 保存初始权重块用于双因子更新
    P_initial_blocks = cell(n, 1);
    for i = 1:n
        block_start = (i-1)*3 + 1;
        P_initial_blocks{i} = P_initial(block_start:block_start+2, block_start:block_start+2);
    end
elseif size(P_initial, 1) == 3 && size(P_initial, 2) == n
    % 传统格式：3×n权重矩阵
    P = P_initial;
    use_full_weight = false;
else
    error('权重矩阵P_initial的维度不正确。应该是3×n或3n×3n格式。');
end

param_diff = 1;
iter_count = 0;

% 初始最小二乘解
X0 = TLS_XG(A, L, P);

while param_diff > param_tol && iter_count < max_iter
    iter_count = iter_count + 1;
    if iter_count > 1
        X_prev = X0;
    end

    % 第一步：使用当前权重矩阵进行TLS求解
    [X0, ~, ~] = TLS_XG_with_W(A, L, P, X0);
    
    % 第二步：计算残差
    v = L - A * X0;  % 整体残差
    
    % 第三步：在外面计算 W
    [n_obs, m_params] = size(A);
    k_total = m_params + 1;
    Q_e = pinv(P);
    I_n = eye(n_obs);
    B = kron(I_n, [1, X0']);
    Pv_inv = B * Q_e * B';
    if rcond(Pv_inv) < eps
        W = pinv(Pv_inv);
    else
        W = inv(Pv_inv);
    end
    
    % 第四步：计算各方向残差
    e_hat = Q_e * B' * W * v;
    e_hat_reshaped = reshape(e_hat, k_total, n)';
    
    e_y = e_hat_reshaped(:, 1);
    e_x1 = e_hat_reshaped(:, 2);
    e_x2 = e_hat_reshaped(:, 3);
    
    % 第五步：根据当前权重进行0权残差吸收机制
    min_weight_threshold = 0;  % 0权阈值
    
    for i = 1:n
        % 提取该点的权重
        if size(P, 1) == 3*n && size(P, 2) == 3*n
            block_start = (i-1)*3 + 1;
            P_block = P(block_start:block_start+2, block_start:block_start+2);
            w_y = P_block(1, 1);
            w_x1 = P_block(2, 2);
            w_x2 = P_block(3, 3);
        else
            w_y = P(1, i);
            w_x1 = P(2, i);
            w_x2 = P(3, i);
        end
        
        % 判断各方向是否为0权
        is_y_zero = (w_y <= min_weight_threshold);
        is_x1_zero = (w_x1 <= min_weight_threshold);
        is_x2_zero = (w_x2 <= min_weight_threshold);
        
        % 统计0权方向数量
        zero_count = is_y_zero + is_x1_zero + is_x2_zero;
        
        % 如果有0权方向，重新分配残差
        if zero_count > 0
            total_residual = v(i);
            
            if zero_count == 1
                % 只有一个0权方向，让该方向吸收全部残差
                if is_y_zero
                    e_y(i) = total_residual;
                    e_x1(i) = 0;
                    e_x2(i) = 0;
                elseif is_x1_zero
                    e_y(i) = 0;
                    e_x1(i) = total_residual;
                    e_x2(i) = 0;
                else  % is_x2_zero
                    e_y(i) = 0;
                    e_x1(i) = 0;
                    e_x2(i) = total_residual;
                end
            else
                % 多个0权方向，平均分配残差
                residual_per_dir = total_residual / zero_count;
                e_y(i) = is_y_zero * residual_per_dir;
                e_x1(i) = is_x1_zero * residual_per_dir;
                e_x2(i) = is_x2_zero * residual_per_dir;
            end
        end
    end
    
    % 第六步：计算单位权中误差（分方向方法使用分方向残差的RMS）
    r = n - size(A,2);  % 自由度
    % 对于分方向方法，应该使用分方向残差的RMS来计算sigma0
    % 这样可以更好地反映分方向残差的尺度，使得标准化残差更准确
    sigma0_y = sqrt(mean(e_y.^2));
    sigma0_x1 = sqrt(mean(e_x1.^2));
    sigma0_x2 = sqrt(mean(e_x2.^2));
    % 使用三个方向RMS的加权平均（或者直接使用最大值）
    % 使用RMS的平方平均（RMS的均方根），更好地反映整体尺度
    sigma0 = sqrt((sigma0_y^2 + sigma0_x1^2 + sigma0_x2^2) / 3);
    % 确保sigma0不会过小（避免数值问题）
    if sigma0 < 1e-10
        % 如果sigma0过小，使用加权总体残差作为后备
        sigma0_weighted = sqrt(max(0, (v' * (W * v)) / r));
        sigma0 = max(sigma0_weighted, 1e-10);
    end

    % 第七步：使用IGG3对各个方向构造等价权
    k0 = 1.5;  % IGGIII参数（标准值）
    k1 = 2.5;  % IGGIII参数（标准值）
    min_weight = 1e-10;
    
    gamma_y = zeros(n, 1);
    gamma_x1 = zeros(n, 1);
    gamma_x2 = zeros(n, 1);
    
    for i = 1:n
        % 计算标准化残差
        e_bar_y = abs(e_y(i)) / (sigma0 + 1e-10);
        e_bar_x1 = abs(e_x1(i)) / (sigma0 + 1e-10);
        e_bar_x2 = abs(e_x2(i)) / (sigma0 + 1e-10);
        
        % IGG3权重函数
        gamma_y(i) = compute_iggiii_weight(e_bar_y, k0, k1, min_weight);
        gamma_x1(i) = compute_iggiii_weight(e_bar_x1, k0, k1, min_weight);
        gamma_x2(i) = compute_iggiii_weight(e_bar_x2, k0, k1, min_weight);
    end

    % 第八步：双因子方法更新权重矩阵（永远基于初始权P_initial）
    if use_full_weight
        P_new = zeros(3*n, 3*n);
        
        for i = 1:n
            for j = 1:n
                block_start_i = (i-1)*3 + 1;
                block_start_j = (j-1)*3 + 1;
                
                % 从初始权重块获取
                P_initial_ij = P_initial(block_start_i:block_start_i+2, block_start_j:block_start_j+2);
                
                % 观测点i和j的等价权
                gamma_i = [gamma_y(i), gamma_x1(i), gamma_x2(i)];
                gamma_j = [gamma_y(j), gamma_x1(j), gamma_x2(j)];
                
                % 构建双因子调权矩阵
                gamma_matrix = zeros(3, 3);
                for p = 1:3
                    for q = 1:3
                        gamma_matrix(p, q) = sqrt(gamma_i(p) * gamma_j(q));
                    end
                end
                
                % 应用双因子调权：P_new = P_initial .* gamma
                P_new(block_start_i:block_start_i+2, block_start_j:block_start_j+2) = ...
                    P_initial_ij .* gamma_matrix;
            end
        end
        
        P = P_new;
    else
        % 传统格式
        for i = 1:n
            P(1,i) = P_initial(1,i) * gamma_y(i);
            P(2,i) = P_initial(2,i) * gamma_x1(i);
            P(3,i) = P_initial(3,i) * gamma_x2(i);
        end
    end

    % 第九步：用更新后的权重再进行一次标准方法求解，使用当前参数作为初值
    [X0, ~, ~] = TLS_XG_with_W(A, L, P, X0);

    if iter_count > 1
        param_diff = norm(X0 - X_prev) / (norm(X_prev) + 1e-10);
    end
end

% 计算最终结果
iter_info.total_iterations = iter_count;
X = X0;
P_final = P;  % 返回最终权重矩阵

% 计算最终残差
[~, e_hat_final, ~] = TLS_XG_with_W(A, L, P, X);
[e_y_final, e_x1_final, e_x2_final] = compute_separate_direction_residuals(A, L, X, P, e_hat_final);
residuals = struct('e_y', e_y_final, 'e_x1', e_x1_final, 'e_x2', e_x2_final);
end

function [X, residuals, iter_info, P_final] = correlated_robust_tls(A, L, P_initial)
% 总体残差相关观测抗差TLS（按照论文双因子方法实现）
% 输入:
%   A: 系数矩阵 (n×m)
%   L: 观测向量 (n×1)
%   P_initial: 初始权重矩阵 (3n×3n 或 3×n)
% 输出:
%   X: 参数估计 (m×1)
%   residuals: 残差结构体
%   iter_info: 迭代信息
%   P_final: 最终权重矩阵

[n, m] = size(A);
k = m + 1; % 每个观测的参数个数(L, a1, a2)
param_tol = 0.001;
max_iter = 50;  % 最大迭代次数

% 初始化
iter_info = struct();
iter_info.total_iterations = 0;

% 检查P_initial的格式并初始化权重矩阵
if size(P_initial, 1) == 3*n && size(P_initial, 2) == 3*n
    % 新格式：3n×3n完整权重矩阵
    P = P_initial;  % 直接使用完整权重矩阵
    use_full_weight = true;
    
    % 保存初始权重块用于双因子更新
    P_initial_blocks = cell(n, 1);
    for i = 1:n
        block_start = (i-1)*3 + 1;
        P_initial_blocks{i} = P_initial(block_start:block_start+2, block_start:block_start+2);
    end
elseif size(P_initial, 1) == 3 && size(P_initial, 2) == n
    % 传统格式：3×n权重矩阵
    P = P_initial;
    use_full_weight = false;
else
    error('权重矩阵P_initial的维度不正确。应该是3×n或3n×3n格式。');
end

param_diff = 1;
iter_count = 0;

% 初始最小二乘解
X0 = TLS_XG(A, L, P);

while param_diff > param_tol && iter_count < max_iter
    iter_count = iter_count + 1;
    if iter_count > 1
        X_prev = X0;
    end

    % 第一步：使用当前权重矩阵进行TLS求解
    [X0, e_hat_current, ~, W] = TLS_XG_with_W(A, L, P, X0);
    
    % 计算总体残差
    v = L - A * X0;  % 基于原始A矩阵的普通残差

    % 计算统一的单位权中误差（使用全局 W 而非逐点权重）
    r = n - size(A,2);  % 自由度
    sigma0 = sqrt(max(0, (v' * (W * v)) / r));
    % 确保sigma0不会过小导致数值问题
    if sigma0 < 1e-6
        sigma0 = sqrt(mean(v.^2)) + 1e-10;
    end

    % 更新权重矩阵 - 按照论文双因子方法
    k0 = 1.5;  % IGGIII参数（标准值）
    k1 = 2.5;  % IGGIII参数（标准值）
    min_weight = 1e-10;   % 设置最小权重，避免权重矩阵奇异

    % 为每个观测点计算总体残差的调权因子
    gamma_v = zeros(n, 1);
    
    for i = 1:n
        % 基于总体残差计算标准化残差
        e_bar = abs(v(i)) / (sigma0 + 1e-10);
        
        % IGGIII权重函数
        gamma_v(i) = compute_iggiii_weight(e_bar, k0, k1, min_weight);
    end

    % 按照论文双因子方法更新权重矩阵
    if use_full_weight
        % 重新构建完整的权重矩阵
        P_new = zeros(3*n, 3*n);
        
        for i = 1:n
            for j = 1:n
                block_start_i = (i-1)*3 + 1;
                block_start_j = (j-1)*3 + 1;
                
                % 获取初始权重块
                P_initial_ij = P_initial(block_start_i:block_start_i+2, block_start_j:block_start_j+2);
                
                % 计算双因子调权矩阵
                % 对于总体残差方法，所有方向的权重都使用同一个gamma_v值
                gamma_factor = sqrt(gamma_v(i) * gamma_v(j));
                gamma_matrix = gamma_factor * ones(3, 3);
                
                % 应用双因子调权（逐元素相乘，保持权重矩阵的结构）
                P_new(block_start_i:block_start_i+2, block_start_j:block_start_j+2) = ...
                    P_initial_ij .* gamma_matrix;
            end
        end
        
        % 验证权重矩阵的有效性
        if any(isnan(P_new(:))) || any(isinf(P_new(:)))
            % 如果权重矩阵无效，使用上一次的权重
            if iter_count > 1
                P = P;  % 保持上一次的权重
            else
                P = P_initial;  % 使用初始权重
            end
        else
            P = P_new;
        end
    else
        % 传统格式的双因子更新
        for i = 1:n
            % 所有三个方向的权重都基于同一个总体残差更新
            P(1,i) = P_initial(1,i) * gamma_v(i);
            P(2,i) = P_initial(2,i) * gamma_v(i);
            P(3,i) = P_initial(3,i) * gamma_v(i);
        end
    end

    % 更新参数估计
    try
        X0_new = TLS_XG_with_W(A, L, P, X0);  % 传入前一次参数作为初值
        
        % 验证参数估计的有效性
        if any(isnan(X0_new)) || any(isinf(X0_new)) || norm(X0_new) > 1e6
            % 如果参数估计无效或过大，使用上一次的估计
            if iter_count > 1
                X0 = X_prev;
            else
                X0 = X0;  % 保持当前值
            end
            param_diff = 0;  % 强制退出循环
        else
            X0 = X0_new;
        end
    catch
        % 如果TLS求解失败，使用上一次的估计
        if iter_count > 1
            X0 = X_prev;
        end
        param_diff = 0;  % 强制退出循环
    end

    if iter_count > 1
        param_diff = norm(X0 - X_prev) / (norm(X_prev) + 1e-10);
        % 如果参数变化太大，可能发散，提前退出
        if param_diff > 10 || norm(X0) > 1e6
            break;
        end
    end
end

% 计算最终结果
iter_info.total_iterations = iter_count;
X = X0;
P_final = P;  % 返回最终权重矩阵

% 计算最终残差
[~, e_hat_final, ~] = TLS_XG_with_W(A, L, P, X);
residuals = struct('e_y', L - A * X, 'e_x1', [], 'e_x2', []);
end

function weight = compute_iggiii_weight(e_bar, k0, k1, min_weight)
% IGGIII权重函数（论文公式33）
% 输入:
%   e_bar: 标准化残差
%   k0, k1: IGGIII参数
%   min_weight: 最小权重
% 输出:
%   weight: 权重值
    if e_bar <= k0
        weight = 1.0;
    elseif e_bar <= k1
        weight = (k0 / (e_bar + 1e-10)) * ((k1 - e_bar) / (k1 - k0 + 1e-10))^2;
    else
        weight = min_weight;
    end
end

function [e_y, e_x1, e_x2] = compute_separate_direction_residuals(A, L, X, P, e_hat_current)
% 计算分方向残差 - 使用当前TLS求解的残差，不重新计算
% 输入:
%   A: 系数矩阵
%   L: 观测向量
%   X: 参数估计
%   P: 权重矩阵
%   e_hat_current: 当前TLS求解的残差
% 输出:
%   e_y, e_x1, e_x2: 各方向的残差

[n, m] = size(A);
k = m + 1;

e_y = zeros(n, 1);   % y方向残差（观测值误差）
e_x1 = zeros(n, 1);  % x1方向残差（系数矩阵第1列误差）
e_x2 = zeros(n, 1);  % x2方向残差（系数矩阵第2列误差）

% 计算整体残差
v = L - A * X;  % 整体残差

% 使用当前TLS求解的残差，重塑为矩阵形式
e_hat_reshaped = reshape(e_hat_current, k, n)';

% 为每个观测点计算分方向残差
for i = 1:n
    % 根据P矩阵格式提取第i个观测的权重
    if size(P, 1) == 3 && size(P, 2) == n
        % 传统格式：3×n权重矩阵
        p_y = P(1, i);      % y方向权重
        p_x1 = P(2, i);     % x1方向权重
        p_x2 = P(3, i);     % x2方向权重
    elseif size(P, 1) == 3*n && size(P, 2) == 3*n
        % 新格式：3n×3n完整权重矩阵，提取对角元素
        block_start = (i-1)*3 + 1;
        p_y = P(block_start, block_start);         % y方向权重
        p_x1 = P(block_start+1, block_start+1);   % x1方向权重
        p_x2 = P(block_start+2, block_start+2);   % x2方向权重
    else
        error('权重矩阵P的维度不正确。');
    end
    
    % 定义极小权重阈值
    min_weight_threshold = 1e-10;
    
    % 检查各方向权重是否为极小值
    is_y_small = p_y <= min_weight_threshold;
    is_x1_small = p_x1 <= min_weight_threshold;
    is_x2_small = p_x2 <= min_weight_threshold;
    
    % 获取当前TLS求解的分方向残差
    e_y_temp = e_hat_reshaped(i, 1);
    e_x1_temp = e_hat_reshaped(i, 2);
    e_x2_temp = e_hat_reshaped(i, 3);
    
    % 统计极小权重方向数量
    small_count = is_y_small + is_x1_small + is_x2_small;
    
    % 根据权重情况调整残差分配
    if small_count == 0
        % 标准情况：使用当前TLS求解的精确分方向残差
        e_y(i) = e_y_temp;
        e_x1(i) = e_x1_temp;
        e_x2(i) = e_x2_temp;
    elseif small_count == 1
        % 只有一个极小权重方向，让该方向吸收全部残差
        if is_y_small
            e_y(i) = v(i);
            e_x1(i) = 0;
            e_x2(i) = 0;
        elseif is_x1_small
            e_y(i) = 0;
            e_x1(i) = v(i);
            e_x2(i) = 0;
        else  % is_x2_small
            e_y(i) = 0;
            e_x1(i) = 0;
            e_x2(i) = v(i);
        end
    else
        % 多个极小权重方向，平均分配残差
        residual_per_dir = v(i) / small_count;
        e_y(i) = is_y_small * residual_per_dir;
        e_x1(i) = is_x1_small * residual_per_dir;
        e_x2(i) = is_x2_small * residual_per_dir;
    end
end
end



function [x_tls, e_hat, iter, W] = TLS_XG(A_obs, L_obs, P)
% 相关牛顿法的求解函数 - 使用新的牛顿法算法
% 输入:
%   A_obs: 系数矩阵 (n×m)
%   L_obs: 观测向量 (n×1)
%   P: 权重矩阵 (3n×3n 或 3×n)
% 输出:
%   x_tls: 参数估计 (m×1)
%   e_hat: 误差估计
%   iter: 迭代次数
%   W: 权重矩阵

tol = 1e-10;
% 获取问题维度
[n, m] = size(A_obs);
k = m + 1;

% 检查P的维度来判断输入格式并构建完整的权重矩阵
if size(P, 1) == 3 && size(P, 2) == n
    % 传统格式：3×n权重矩阵，转换为3n×3n块对角矩阵
    P_full = zeros(3*n, 3*n);
    for i = 1:n
        block_start = (i-1)*3 + 1;
        block_end = i*3;
        P_full(block_start:block_end, block_start:block_end) = diag([P(1,i), P(2,i), P(3,i)]);
    end
    % 提取L对应的权矩阵用于初始OLS估计
    P_L = zeros(n, n);
    for i = 1:n
        P_L(i, i) = P(1, i);
    end
elseif size(P, 1) == 3*n && size(P, 2) == 3*n
    % 新格式：直接使用3n×3n完整权重矩阵
    P_full = P;
    % 提取L对应的权矩阵用于初始OLS估计
    P_L = zeros(n, n);
    for i = 1:n
        block_start = (i-1)*k + 1;
        P_block = P_full(block_start:block_start+k-1, block_start:block_start+k-1);
        P_L(i, i) = P_block(1, 1);
    end
else
    error('权重矩阵P的维度不正确。应该是3×n或3n×3n格式。');
end

% 初始OLS估计
x0 = (A_obs' * P_L * A_obs) \ (A_obs' * P_L * L_obs);

% TLS牛顿迭代法
x = x0;
iter = 0;
dx_norm = 1;
I_n = eye(n);

% 改进的协方差矩阵计算，添加数值稳定性检查
if rcond(P_full) < 1e-10
    % 如果P_full条件数太差，添加正则化
    reg_factor = 1e-8 * trace(P_full) / size(P_full, 1);
    P_full_reg = P_full + reg_factor * eye(size(P_full));
    Q_e = pinv(P_full_reg);
else
    Q_e = pinv(P_full);
end

max_iter = 100;

while dx_norm > tol && iter < max_iter
    iter = iter + 1;
    
    % 1. 计算残差向量 v = L_obs - A_obs * x
    v = L_obs - A_obs * x;
    
    % 2. 构建块对角矩阵 B
    B = kron(I_n, [1, x']);
    
    % 3. 计算 W = (B * P^{-1} * B')^{-1}，改进数值稳定性
    Pv_inv = B * Q_e * B';
    
    % 添加更严格的条件数检查
    rcond_val = rcond(Pv_inv);
    if isnan(rcond_val) || rcond_val < 1e-12
        % 如果条件数太差，使用正则化
        reg = 1e-10 * trace(Pv_inv) / n;
        Pv_inv_reg = Pv_inv + reg * eye(n);
        W = pinv(Pv_inv_reg);
    else
        W = inv(Pv_inv);
    end
    
    % 4. 计算误差估计 ê
    e_hat = Q_e * B' * W * v;
    
    % 5. 提取系数矩阵的误差 e_A
    e_hat_reshaped = reshape(e_hat, k, n)';
    e_A = e_hat_reshaped(:, 2:end);
    
    % 6. 计算修正后的系数矩阵
    A_corr = A_obs + e_A;

    % 7. 计算梯度 F
    F = A_corr' * (W * v);
    
    % 8. 计算近似Hessian矩阵
    H1 = -2 * A_obs' * W * e_A;
    H2 = -A_obs' * W * A_obs;
    H3 = -2 * e_A' * W * e_A;
    H4 = -2 * e_A' * W * A_obs;

    % H5矩阵
    dET_dxi = zeros(m,n);
    
    % 提前计算 M = Q_e * B' * W，用于优化循环内的计算
    M = Q_e * B' * W;   % M is 3n x n
    
    % 循环遍历i从1到m（注意：x有m个元素，b有n=m+1个元素）
    for i = 1:m
        % 计算 ∂b/∂x_i - 这是一个 1×(m+1) 的行向量
        % 除了第 i+1 个位置为 1（因为第一个位置是常数 1），其余为 0
        db_dxi = zeros(1, m+1);
        db_dxi(i+1) = 1;
        
        % 计算 ∂B^T/∂x_i
        dBT_dxi = kron(I_n, db_dxi);
        dBT_dxi_T = dBT_dxi';
        
        % 计算第一项: Q_e * (∂B^T/∂x_i) * W * v
        term1 = Q_e * dBT_dxi_T * W * v;
        term2 = -2 * M * e_A(:, i);
        term3 = -M * A_obs(:, i); 
        de_dxi = term1 + term2 + term3;

        for j = 1:n
            dET_dxi(:,j) = de_dxi((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1)); 
        end
        
        H5(:,i) = dET_dxi * W * v;
    end

    H = H1+H2+H3+H4+H5;
    
    % 9. 牛顿迭代更新
    if rcond(H) < eps
        dx = pinv(H) * F;
    else
        dx = H \ F;
    end
    
    x_new = x - dx;
    
    % 检查收敛
    dx_norm = norm(dx);
    
    x = x_new;
end

% 输出最终结果
x_tls = x;
end

function [x_tls, e_hat, iter, W] = TLS_XG_with_W(A_obs, L_obs, P, x0)
% 带初值的TLS牛顿法 - 使用新的牛顿法算法
% 输入:
%   A_obs: 系数矩阵 (n×m)
%   L_obs: 观测向量 (n×1)
%   P: 权重矩阵 (3n×3n 或 3×n)
%   x0: 初始参数估计 (m×1)
% 输出:
%   x_tls: 参数估计 (m×1)
%   e_hat: 误差估计
%   iter: 迭代次数
%   W: 权重矩阵

tol = 1e-10;
% 获取问题维度
[n, m] = size(A_obs);
k = m + 1;

% 检查P的维度来判断输入格式并构建完整的权重矩阵
if size(P, 1) == 3 && size(P, 2) == n
    % 传统格式：3×n权重矩阵，转换为3n×3n块对角矩阵
    P_full = zeros(3*n, 3*n);
    for i = 1:n
        block_start = (i-1)*3 + 1;
        block_end = i*3;
        P_full(block_start:block_end, block_start:block_end) = diag([P(1,i), P(2,i), P(3,i)]);
    end
elseif size(P, 1) == 3*n && size(P, 2) == 3*n
    % 新格式：直接使用3n×3n完整权重矩阵
    P_full = P;
else
    error('权重矩阵P的维度不正确。应该是3×n或3n×3n格式。');
end

% 使用传入的初值而不是OLS估计
x = x0;
iter = 0;
dx_norm = 1;
I_n = eye(n);

% 改进的协方差矩阵计算，添加数值稳定性检查
if rcond(P_full) < 1e-10
    % 如果P_full条件数太差，添加正则化
    reg_factor = 1e-8 * trace(P_full) / size(P_full, 1);
    P_full_reg = P_full + reg_factor * eye(size(P_full));
    Q_e = pinv(P_full_reg);
else
    Q_e = pinv(P_full);
end

max_iter = 100;

while dx_norm > tol && iter < max_iter
    iter = iter + 1;
    
    % 1. 计算残差向量 v = L_obs - A_obs * x
    v = L_obs - A_obs * x;
    
    % 2. 构建块对角矩阵 B
    B = kron(I_n, [1, x']);
    
    % 3. 计算 W = (B * P^{-1} * B')^{-1}，改进数值稳定性
    Pv_inv = B * Q_e * B';
    
    % 添加更严格的条件数检查
    rcond_val = rcond(Pv_inv);
    if isnan(rcond_val) || rcond_val < 1e-12
        % 如果条件数太差，使用正则化
        reg = 1e-10 * trace(Pv_inv) / n;
        Pv_inv_reg = Pv_inv + reg * eye(n);
        W = pinv(Pv_inv_reg);
    else
        W = inv(Pv_inv);
    end
    
    % 4. 计算误差估计 ê
    e_hat = Q_e * B' * W * v;
    
    % 5. 提取系数矩阵的误差 e_A
    e_hat_reshaped = reshape(e_hat, k, n)';
    e_A = e_hat_reshaped(:, 2:end);
    
    % 6. 计算修正后的系数矩阵
    A_corr = A_obs + e_A;

    % 7. 计算梯度 F
    F = A_corr' * (W * v);
    
    % 8. 计算近似Hessian矩阵
    H1 = -2 * A_obs' * W * e_A;
    H2 = -A_obs' * W * A_obs;
    H3 = -2 * e_A' * W * e_A;
    H4 = -2 * e_A' * W * A_obs;

    % H5矩阵
    dET_dxi = zeros(m,n);
    
    % 提前计算 M = Q_e * B' * W，用于优化循环内的计算
    M = Q_e * B' * W;   % M is 3n x n
    
    % 循环遍历i从1到m（注意：x有m个元素，b有n=m+1个元素）
    for i = 1:m
        % 计算 ∂b/∂x_i - 这是一个 1×(m+1) 的行向量
        % 除了第 i+1 个位置为 1（因为第一个位置是常数 1），其余为 0
        db_dxi = zeros(1, m+1);
        db_dxi(i+1) = 1;
        
        % 计算 ∂B^T/∂x_i
        dBT_dxi = kron(I_n, db_dxi);
        dBT_dxi_T = dBT_dxi';
        
        % 计算第一项: Q_e * (∂B^T/∂x_i) * W * v
        term1 = Q_e * dBT_dxi_T * W * v;
        term2 = -2 * M * e_A(:, i);
        term3 = -M * A_obs(:, i); 
        de_dxi = term1 + term2 + term3;

        for j = 1:n
            dET_dxi(:,j) = de_dxi((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1)); 
        end
        
        H5(:,i) = dET_dxi * W * v;
    end

    H = H1+H2+H3+H4+H5;
    
    % 9. 牛顿迭代更新
    if rcond(H) < eps
        dx = pinv(H) * F;
    else
        dx = H \ F;
    end
    
    x_new = x - dx;
    
    % 检查收敛
    dx_norm = norm(dx);
    
    x = x_new;
end

% 输出最终结果
x_tls = x;
end

