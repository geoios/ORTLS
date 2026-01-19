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
diag_newton_all = zeros(m, num_rhos);
diag_Jaz_all = zeros(m, num_rhos);
diag_Wang_all = zeros(m, num_rhos);
diag_Xu_all = zeros(m, num_rhos);
diag_numerical_all = zeros(m, num_rhos);

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

    Py = inv(Q_L);
    PA = inv(QA1);
    PA2 = inv(QA2);
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
    e_hat = e_hat_reshaped(:);
    
    % 计算J矩阵
    J_total = zeros(m, k*n);

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

    J0 =dx_dL;
    % 4.计算Ji
    for param_idx = 1:m  % 对每个参数列（A1, A2, ...）
        dF_da_i = zeros(m, n);
        for obs_idx = 1:n  % 对每个观测值
            % 计算R1和R2矩阵（针对特定参数和观测值）
            R1 = zeros(n, m);
            R1(obs_idx, param_idx) = 1;  % 第obs_idx行，第param_idx列设为1

            R2 = zeros(n, 1);
            R2(obs_idx) = -X_newton(param_idx);  % -X(param_idx)

            % 计算∂e/∂a_{param_idx,obs_idx}
            de_da = Q_e * B' * P_v * R2;
            dET_da = zeros(m, n);
            for j = 1:n
                dET_da(:, j) = de_da((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1));
            end
            dF_da_single = (R1' + dET_da) * P_v * v_newton + (A_obs + e_A)' * P_v * R2;
            dF_da_i(:, obs_idx) = dF_da_single;
        end
        if rcond(H) < eps
            dx_da_i = -pinv(H) * dF_da_i;  % m×n
        else
            dx_da_i = -H \ dF_da_i;        % 更高效，等价于 -inv(H) * de_da_i
        end
        Ja{param_idx} = dx_da_i;
    end
    J_total = [J0, Ja{1}, Ja{2}];

    % 单位权方差
    sit0_1 = (e' * PP * e) / (n - m);
    Dx_newton =   J_total * Q_total * J_total';
    diag_newton_all(1:m, rho_idx) = diag(Dx_newton);
    
    %% ============== Wang法 ==============
    [X, sigma2, Qxx, ci, e_W] = Wang(n, A_obs, L_obs, Qcov_A1L, QA1, Q_L, 0);
    Dx_Wang = Qxx;
    diag_Wang_all(1:m, rho_idx) = diag(Dx_Wang);

    %% ============== Xu法 ==============
    P(1,:)=diag(Py);
    P(2,:)=diag(PA);
    P(3,:)=diag(PA2);
    [beta_hat, a_hat, iter, conv_flag, S0, A0, PA] = Xu_m(A_obs, L_obs, P);
    A = A_obs;
    L = L_obs;
    a_hat0 = A(:,1);
    r_a = a_hat0(:) - a_hat;            % 随机元素的残差（使用观测值a）
    r_y = L_obs - A0 * beta_hat;       % 观测向量的残差

    % 计算单位权方差
    sigma2_hat = (r_a' * PA * r_a + r_y' * Py * r_y) / (n - m);
    sigma2_hat = 1;
    
    % 计算法方程矩阵
    N11 = A0' * Py * A0;
    N12 = A0' * Py * S0;
    N21 = S0' * Py * A0;
    N22 = S0' * Py * S0 + PA;

    % 计算参数β的方差-协方差矩阵
    Dx_xu = inv(N11 - N12 * inv(N22) * N21) ;
    diag_Xu_all(1:m, rho_idx) = diag(Dx_xu);

    %% ============== jazaeri法 ==============
    fprintf('  计算传统方法 (WTLS)...\n');
    % 初始化：使用普通最小二乘(OLS)
    
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
    Q_vv =  Q_x;
    diag_Jaz_all(1:m, rho_idx) = diag(Q_vv);
    
    %% ============== 数值法 (蒙特卡洛) ==============
    fprintf('  计算数值法 (蒙特卡洛)...\n');
    a_sim = zeros(num_monte_carlo, m);
    
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
        a_sim(i, :) = X_sim;
    end
    
    cov_num = cov(a_sim);  % 数值法协方差矩阵
    diag_numerical_all(1:m, rho_idx) = diag(cov_num);
    
    fprintf('  相关系数 rho = %.1f 完成\n\n', rho);
end

%% ============== Plotting Section (Two Separate Figures) ==============
fprintf('=== Starting Plotting (Two Separate Figures) ===\n');

% Method names (now including MCS)
method_names_5 = {'Newton', 'Wang', 'Xu', 'Jazaeri', 'MCS'};

%% Figure 1: Parameter 1 (Slope x_1) Variance Comparison
fig1 = figure('Position', [100, 100, 1200, 400], 'Name', 'Parameter 1 (Slope) Variance Comparison');
set(fig1, 'Color', 'w');

% Create subplots for each correlation coefficient
for rho_idx = 1:num_rhos
    rho = rho_values(rho_idx);
    
    subplot(1, 3, rho_idx);
    
    % Prepare data for parameter 1 (5 methods including MCS)
    data_methods = [diag_newton_all(1, rho_idx); 
                    diag_Wang_all(1, rho_idx);
                    diag_Xu_all(1, rho_idx);
                    diag_Jaz_all(1, rho_idx);
                    diag_numerical_all(1, rho_idx)];  % Add MCS
    
    % Create grouped bar chart for 5 methods
    x_positions = 1:5;
    hb = bar(x_positions, data_methods);
    set(gca, 'Color', 'w');
    hold on;
    
    % 设置不同方法的颜色：前4种方法使用蓝色，MCS使用红色
    colors = zeros(5, 3);  % 创建颜色矩阵
    colors(1:4, :) = repmat([0.2, 0.6, 0.8], 4, 1);  % 前4种方法使用蓝色
    colors(5, :) = [0.8, 0.2, 0.2];  % MCS使用红色
    
    % 设置每个柱子的颜色
    hb.FaceColor = 'flat';  % 启用独立颜色设置
    hb.CData = colors;      % 设置每个柱子的颜色
    hb.EdgeColor = 'k';
    hb.LineWidth = 1;
    
    % Customize plot
    title(sprintf('ρ = %.1f', rho), 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Variance', 'FontSize', 11);
    xticks(1:5);
    xticklabels(method_names_5);
    xtickangle(45);
    grid on;
    box on;
    
    % Add value labels on bars
    for i = 1:5
        text(i, data_methods(i), sprintf('%.2e', data_methods(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
    end
    
    hold off;
end

% Add overall title
sgtitle('Parameter 1 (Slope x₁) Variance Estimation Comparison', 'FontSize', 14, 'FontWeight', 'bold');



%% Figure 2: Parameter 2 (Intercept x_0) Variance Comparison
fig2 = figure('Position', [100, 100, 1200, 400], 'Name', 'Parameter 2 (Intercept) Variance Comparison');
set(fig2, 'Color', 'w');

% Create subplots for each correlation coefficient
for rho_idx = 1:num_rhos
    rho = rho_values(rho_idx);
    
    subplot(1, 3, rho_idx);
    
    % Prepare data for parameter 2 (5 methods including MCS)
    data_methods = [diag_newton_all(2, rho_idx); 
                    diag_Wang_all(2, rho_idx);
                    diag_Xu_all(2, rho_idx);
                    diag_Jaz_all(2, rho_idx);
                    diag_numerical_all(2, rho_idx)];  % Add MCS
    
    % Create grouped bar chart for 5 methods
    x_positions = 1:5;
    hb = bar(x_positions, data_methods);
    set(gca, 'Color', 'w');
    hold on;
    
    % 设置不同方法的颜色：前4种方法使用蓝色，MCS使用红色
    colors = zeros(5, 3);  % 创建颜色矩阵
    colors(1:4, :) = repmat([0.2, 0.6, 0.8], 4, 1);  % 前4种方法使用蓝色
    colors(5, :) = [0.8, 0.2, 0.2];  % MCS使用红色
    
    % 设置每个柱子的颜色
    hb.FaceColor = 'flat';  % 启用独立颜色设置
    hb.CData = colors;      % 设置每个柱子的颜色
    hb.EdgeColor = 'k';
    hb.LineWidth = 1;
    
    % Customize plot
    title(sprintf('ρ = %.1f', rho), 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Variance', 'FontSize', 11);
    xticks(1:5);
    xticklabels(method_names_5);
    xtickangle(45);
    grid on;
    box on;
    
    % Add value labels on bars
    for i = 1:5
        text(i, data_methods(i), sprintf('%.2e', data_methods(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
    end
    
    hold off;
end

% Add overall title
sgtitle('Parameter 2 (Intercept x₀) Variance Estimation Comparison', 'FontSize', 14, 'FontWeight', 'bold');

function [X, sigma2, Qxx, ci, e_full] = Wang(n, A, L, Qcov_A1L, QA1, Q_L, is_random_col2)
% 修正的相关观测Partial EIV模型求解
% Q_full: 完整协因数矩阵，包含相关性信息

Q_full=[Q_L, Qcov_A1L';
       Qcov_A1L, QA1];

% 提取随机元素
if is_random_col2
    a = reshape(A, 2*n, 1);  % 所有元素随机
    t = 2*n;
    h = zeros(2*n, 1);
    B = eye(2*n);
else
    a = A(:,1);  % 仅第一列随机
    t = n;
    h = [zeros(n,1); ones(n,1)];
    B = kron([1;0], eye(n));
end

% 权矩阵
W = inv(Q_full);

% 最小二乘初值
Py = inv(Q_L);  % 观测向量权阵
X0 = pinv(A' * Py * A) * A' * Py * L;

% 初始化
I = eye(n);
C1 = eye(n + t);
ci = 0;
cita = 1;

% 迭代求解
while cita > 1e-10
    ci = ci + 1;
    
    % 构造C2
    C2 = [-kron(X0', I) * B, I] * C1;
    
    % 计算Q1
    Q1 = C2 * Q_full * C2';
    
    % 计算拉格朗日乘子
    lambda = Q1 \ (L - kron(X0', I) * (h + B * a));
    
    % 计算残差
    e_temp = -Q_full * C2' * lambda;
    eA = e_temp(1:t);
    
    % 更新
    EA = reshape(B * eA, n, 2);
    A0 = A + EA;
    y1 = L + EA * X0;
    
    % 新参数估计
    X = pinv(A0' / Q1 * A0) * A0' / Q1 * y1;
    
    % 收敛判断
    cita = norm(X - X0);
    X0 = X;
end

% 最终残差计算
e_full = e_temp;

% 单位权方差
r = (n + t) - (t + 2);  % 自由度
sigma2 = (e_full' * W * e_full) / r;

% 参数协因数阵
Qxx = inv(A0' / Q1 * A0);
end