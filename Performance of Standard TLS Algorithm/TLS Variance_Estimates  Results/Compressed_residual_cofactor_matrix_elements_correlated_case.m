clear all; clc; close all;

% 参数设置
n = 10;                % 观测值数量
m = 2;                 % 参数个数
k = m + 1;             % 每个观测点的变量数(L + A1 + A2)
num_monte_carlo = 1000; % 蒙特卡洛模拟次数

% 生成真实数据
A1_true = (0:9)';
A2_true = ones(n, 1);
L_true = 2*A1_true - A2_true;
A_true = [A1_true, A2_true];

% 设置噪声方差
sigma_x = 0.001;         % x的噪声方差
sigma_L = 0.001;         % L的噪声方差

% 定义三个相关系数
rho_values = [0.3, 0.6, 0.9];
num_rhos = length(rho_values);

% 初始化存储结果的数组
diag_newton_all = zeros(n, num_rhos);
diag_traditional_all = zeros(n, num_rhos);
diag_numerical_all = zeros(n, num_rhos);

% 循环处理三个相关系数
for rho_idx = 1:num_rhos
    rho = rho_values(rho_idx);
    
    fprintf('=== 实验: 使用相关系数 rho = %.1f ===\n', rho);
    
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
    
    % 构建不同观测点间的相关性矩阵（简单结构）
    % 所有不同观测点间的相关系数都相同
    Sigma_between = ones(n, n) * rho;  % 所有元素初始化为rho
    Sigma_between(logical(eye(n))) = 1;  % 对角线设为1
    
    % 构建完整的协方差矩阵（Kronecker乘积方式）
    Sigma_global = kron(Sigma_between, Sigma_point);
    
    % 生成观测数据
    % 注意：A2没有误差，所以只对L和A1添加噪声
    R = chol(Sigma_global);  % Sigma_global = R' * R
    noise = randn(n*total_vars, 1);
    correlated_noise = R' * noise;
    
    % 将噪声添加到真实数据中
    % 首先构造真实数据向量（按观测点排列：L1, A11, A21, L2, A12, A22, ...）
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
        % A2保持原值，没有噪声
        Y_obs(idx+2) = Y_true(idx+2);                           % A2无噪声
    end
    
    % 提取观测值
    L_obs = Y_obs(1:total_vars:end);
    A1_obs = Y_obs(2:total_vars:end);
    A2_obs = Y_obs(3:total_vars:end);
    A_obs = [A1_obs, A2_obs];
    
    % 索引定义
    idx_L = 1:total_vars:n*total_vars;      % L的索引：1, 4, 7, ...
    idx_A1 = 2:total_vars:n*total_vars;     % A1的索引：2, 5, 8, ...
    idx_A2 = 3:total_vars:n*total_vars;     % A2的索引：3, 6, 9, ...
    
    % 提取各变量的协方差子矩阵
    Q_L = Sigma_global(idx_L, idx_L);
    QA1 = Sigma_global(idx_A1, idx_A1);
    QA2 = Sigma_global(idx_A2, idx_A2);
    
    % 提取互协方差子矩阵
    Qcov_A1L = Sigma_global(idx_A1, idx_L);
    Qcov_A2L = Sigma_global(idx_A2, idx_L);
    Qcov_A1A2 = Sigma_global(idx_A1, idx_A2);
    
    % 构建总体协因数阵（分块矩阵形式）
    Q_total = [Q_L, Qcov_A1L', Qcov_A2L';
               Qcov_A1L, QA1, Qcov_A1A2';
               Qcov_A2L, Qcov_A1A2, QA2];
    
    %% ============== 牛顿法 ==============
    fprintf('  计算牛顿法...\n');
    PP = inv(Sigma_global);
    X_newton = TLS_XG_newton3(A_obs, L_obs, PP);
    v_newton = L_obs - A_obs * X_newton;
    
    % 计算误差传播
    Q_e = Sigma_global;
    [H, e_A, B, e] = Hessian(A_obs, L_obs, PP, X_newton);
    e_hat_reshaped = reshape(e, k, n)';
    
    % 计算J矩阵
    J_total = zeros(m*n, k*n);
    
    % 1.计算 P_v = (B * Q_e * B')^{-1}
    Pv_inv = B * Q_e * B';
    if rcond(Pv_inv) < eps
        P_v = pinv(Pv_inv);
    else
        P_v = inv(Pv_inv);
    end
    
    % 2.计算 dx_dL = -H \ (A + e_A)' * P_v;
    Gamma = (A_obs + e_A)' * P_v;
    if rcond(H) < eps
        dx_dL = -pinv(H) * Gamma;  % m×n
    else
        dx_dL = -H \ Gamma;        % 更高效，等价于 -inv(H) * Gamma
    end
    
    % 3.计算J0
    In = eye(n);
    J0 = In - A_obs * dx_dL;
    
    % 4.计算Ji
    Ja = cell(1, m);
    dx_da_all = cell(1, m);
    for param_idx = 1:m
        dF_da_i = zeros(m, n);
        term1 = zeros(n, 1);
        for obs_idx = 1:n
            % 计算R1和R2矩阵
            R1 = zeros(n, m);
            R1(obs_idx, param_idx) = 1;
            
            R2 = zeros(n, 1);
            R2(obs_idx) = -X_newton(param_idx);
            
            % 计算∂e/∂a_{param_idx,obs_idx}
            de_da = Q_e * B' * P_v * R2;
            dET_da = zeros(m, n);
            for j = 1:n
                dET_da(:, j) = de_da((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1));
            end
            dF_da_single = (R1' + dET_da) * P_v * v_newton + (A_obs + e_A)' * P_v * R2;
            dF_da_i(:, obs_idx) = dF_da_single;
            da = zeros(n,m);
            da(obs_idx, param_idx) = 1;
            term1(:, obs_idx) = da * X_newton;
        end
        
        if rcond(H) < eps
            dx_da_i = -pinv(H) * dF_da_i;  % m×n
        else
            dx_da_i = -H \ dF_da_i;
        end
        dx_da_all{param_idx} = dx_da_i;
        
        dv_da = -(term1 + A_obs * dx_da_all{param_idx});
        Ja{param_idx} = dv_da;
    end
    J_total = [J0, Ja{1}, Ja{2}];
    
    % 单位权方差
    sit0_1 = (v_newton' * P_v * v_newton) / (n - m);
    Dx_newton =  J_total * Q_total * J_total';
    diag_newton_all(:, rho_idx) = diag(Dx_newton);
    
    %% ============== jazaeri法 ==============
    fprintf('  计算传统方法 (WTLS)...\n');
    % 初始化：使用普通最小二乘(OLS)
    Py = inv(Q_L);
    x_ols = (A_obs' * Py * A_obs) \ (A_obs' * Py * L_obs);
    
    % 构建协因数阵 - 与传统法一致的块对角结构
    Q_y = Q_L;          % L的协因数阵，n×n
    Q_A1 = QA1;        % A1的协因数阵，n×n
    Q_A2 = QA2;        % A2的协因数阵，n×n
    Q_A = blkdiag(Q_A1, Q_A2);   % 块对角矩阵
    
    % WTLS迭代
    x_hat = x_ols;
    epsilon = 1e-10;
    
    for iter = 1:20
        % 计算残差
        e_hat = L_obs - A_obs * x_hat;
        
        % 计算 Q_y_tilde = Q_y + (x^T ⊗ I_n) Q_A (x ⊗ I_n)
        x_kron_T = kron(x_hat', eye(n));  % n × (m*n)
        x_kron = kron(x_hat, eye(n));      % (m*n) × n
        Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
        Q_y_tilde_inv = inv(Q_y_tilde);
        
        % 计算 E_A_hat = -vec^{-1}(Q_A (x ⊗ I_n) Q_y_tilde^{-1} e_hat)
        vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
        E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];
        
        % 修正后的设计矩阵和观测值
        A_tilde = A_obs - E_A_hat;
        L_tilde = L_obs - E_A_hat * x_hat;
        
        % 更新参数
        x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * L_tilde);
        
        % 收敛检查
        delta = norm(x_hat_new - x_hat);
        
        x_hat = x_hat_new;
        
        if delta < epsilon
            break;
        end
    end
    
    % 计算最终残差
    final_residuals = L_obs - A_obs * x_hat;
    
    % 重新计算Q_y_tilde用于最终参数
    x_kron_T = kron(x_hat', eye(n));
    x_kron = kron(x_hat, eye(n));
    Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
    Q_y_tilde_inv = inv(Q_y_tilde);
    
    % 计算E_A_hat
    vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * final_residuals;
    E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];
    
    % 计算A_tilde
    A_tilde = A_obs - E_A_hat;
    
    % 参数协方差矩阵（不含σ_0²）
    Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);
    sigma_0_sq = (final_residuals' * Q_y_tilde_inv * final_residuals) / (n - m);
    
    % 残差协方差矩阵(归一化，不含σ_0²)
    Q_vv = (Q_y_tilde - A_tilde * Q_x * A_tilde');
    diag_traditional_all(:, rho_idx) = diag(Q_vv);
    
    %% ============== 数值法 (蒙特卡洛) ==============
    fprintf('  计算数值法 (蒙特卡洛)...\n');
    a_sim = zeros(num_monte_carlo, n);
    
    % 构建L和A1的联合协方差矩阵（2n×2n）
    % 从完整的协方差矩阵中提取L和A1的子矩阵
    idx_L = 1:total_vars:n*total_vars;
    idx_A1 = 2:total_vars:n*total_vars;
    
    % 提取L和A1的协方差子矩阵
    Q_L_only = Sigma_global(idx_L, idx_L);
    Q_A1_only = Sigma_global(idx_A1, idx_A1);
    Q_cross_LA1 = Sigma_global(idx_L, idx_A1);
    
    % 构建联合协方差矩阵
    Sigma_LA1 = [Q_L_only, Q_cross_LA1;
                Q_cross_LA1', Q_A1_only];
    
    % 计算Cholesky分解
    R_LA1 = chol(Sigma_LA1);
    
    % 蒙特卡洛模拟
    for i = 1:num_monte_carlo
        % 生成相关随机噪声
        noise_LA1_sim = randn(2*n, 1);
        correlated_noise_sim = R_LA1' * noise_LA1_sim;
        
        % 分离L和A1的噪声
        L_noise_sim = correlated_noise_sim(1:n);
        A1_noise_sim = correlated_noise_sim(n+1:end);
        
        % 生成模拟观测值
        L_sim = L_true + L_noise_sim;
        A1_sim = A1_true + A1_noise_sim;
        A2_sim = A2_true;  % A2没有误差
        A_sim = [A1_sim, A2_sim];
        X_sim = TLS_XG_newton3(A_sim, L_sim, PP);
        e_hat = L_sim - A_sim * X_sim;
        a_sim(i, :) = e_hat;
    end
    
    cov_num = cov(a_sim);  % 数值法协方差矩阵
    diag_numerical_all(:, rho_idx) = diag(cov_num);
    
    fprintf('  相关系数 rho = %.1f 完成\n\n', rho);
end

%% ============== 第一个图框：三种方法的对角线元素对比 ==============
fprintf('绘制第一个图框：三种方法的对角线元素对比...\n');
% 设置16:9比例图框
fig_width = 1200;
fig_height = fig_width * 7/16; % 16:9比例
fig1 = figure('Position', [100, 100, fig_width, 350]);
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% 为每个相关系数创建子图
for rho_idx = 1:num_rhos
    subplot(1, 3, rho_idx); % 1行3列
    
    rho = rho_values(rho_idx);
    
    % 提取当前相关系数的结果
    diag_newton = diag_newton_all(:, rho_idx);
    diag_traditional = diag_traditional_all(:, rho_idx);
    diag_numerical = diag_numerical_all(:, rho_idx);
    
    % 绘制三种方法的对角线元素对比
    plot(1:n, diag_newton, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Newton');
    hold on;
    plot(1:n, diag_traditional, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Jazaeri');
    plot(1:n, diag_numerical, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'MCS');
    hold off;
    
    xlabel('Obs.No.', 'FontSize', 18);
    ylabel('Diagonal Elements', 'FontSize', 18);
    % 添加子图标题并获取句柄
    h1 = title(sprintf('ρ = %.1f', rho), 'FontSize', 18, 'FontWeight', 'bold');



    % 获取当前坐标轴范围并调整标题位置
    ax1 = gca;
    % 使用归一化坐标调整标题位置
    set(h1, 'Units', 'normalized');
    titlePos1 = get(h1, 'Position');

    
    % 只在第一个子图设置图例
    if rho_idx == 1
        legend('Location', 'northwest', 'FontSize', 16);
    end
    
    grid on;
    
    % 设置y轴范围一致以便比较
    y_max = max([diag_newton; diag_traditional; diag_numerical]) * 1.1;
    y_min = min([diag_newton; diag_traditional; diag_numerical]) * 0.9;
    ylim([y_min, y_max]);
end

