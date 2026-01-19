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
diag_newton_all = zeros(n*2, num_rhos);
diag_numerical_all = zeros(n*2, num_rhos);

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
    X = TLS_XG_newton3(A_obs, L_obs, PP);
    v = L_obs - A_obs * X;
    % 计算误差传播
    Q_e = Sigma_global;
    [H, e_A, B, e] = Hessian(A_obs, L_obs, PP, X);
    e_hat_reshaped = reshape(e, k, n)';
    e_L = e_hat_reshaped(:, 1);
    

    % 计算J矩阵[J0,J1,J2]10*30
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
    % 3.计算 J0 = ∂e/∂L;
    C = B' * P_v;
    dC_dx = zeros(n*k*n,m);
    for i = 1:m

        Tk_i = zeros(n, n * (m + 1));
        for j = 1:n
            col_index = (j - 1) * (m + 1) + i + 1;
            Tk_i(j, col_index) = 1;
        end

        dC_dxk{i} = Tk_i' * P_v - 2 * B' * P_v * Tk_i * Q_e * B' * P_v;
        dC_dxkk = dC_dxk{i};
        dC_dx(:,i) = dC_dxkk(:);
    end
    dC_dL_vec = dC_dx * dx_dL;


    for j = 1:n
        dC_dL_j = dC_dL_vec(:,j);
        dC_dL{j} = reshape(dC_dL_j, n*k, n);
        dl=zeros(n,1);
        dl(j)=1;
        de_dL(:, j) = Q_e * (dC_dL{j} * L_obs + C * dl  - dC_dL{j} * A_obs * X - C * A_obs * dx_dL(:,j));
    end
    J0 = de_dL;

    % 4.计算 Ji = ∂e/∂ai;
    % 4.1. 设置初始存储    
    Sigma_e_a_total = zeros(k*n, k*n);
    Ja = cell(1, m);
    dx_da_all = cell(1, m);
    dC_dxk = cell(1, m);

    % 4.2. 对每个设计矩阵参数ai进行循环计算
    for param_idx = 1:m  % 对每个参数列（A1, A2, ...）

        % 对该参数的每个观测值进行循环
        de_da_i = zeros(m, n);  % ∂e/∂a_{param_idx,j} for j=1 to n

        for obs_idx = 1:n  % 对每个观测值
            % 计算R1和R2矩阵（针对特定参数和观测值）
            R1 = zeros(n, m);
            R1(obs_idx, param_idx) = 1;  % 第obs_idx行，第param_idx列设为1

            R2 = zeros(n, 1);
            R2(obs_idx) = -X(param_idx);  % -X(param_idx)

            % 计算∂e/∂a_{param_idx,obs_idx}
            de_da = Q_e * B' * P_v * R2;
                dET_da = zeros(m, n);
            for j = 1:n
                dET_da(:, j) = de_da((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1));
            end
            dF_da_single = (R1' + dET_da) * P_v * v + (A_obs + e_A)' * P_v * R2;
            de_da_i(:, obs_idx) = dF_da_single;
            if rcond(H) < eps
                dx_da_i = -pinv(H) * de_da_i;  % m×n
            else
                dx_da_i = -H \ de_da_i;        % 更高效，等价于 -inv(H) * de_da_i
            end
        end
        dx_da_all{param_idx} = dx_da_i;
        % 计算dC_dx_vec，即对偏导dC_dx进行线性化
        dC_dx = zeros(n*k*n,m);
         for i = 1:m

             Tk_i = zeros(n, n * (m + 1));
             for j = 1:n
                 col_index = (j - 1) * (m + 1) + i + 1;
                 Tk_i(j, col_index) = 1;
             end

             dC_dxk{i} = Tk_i' * P_v - 2 * B' * P_v * Tk_i * Q_e * B' * P_v;
             dC_dxkk = dC_dxk{i};
             dC_dx(:,i) = dC_dxkk(:);
         end
            dC_dx_vec = dC_dx * dx_da_all{param_idx};
        % 计算最终结果Ji
        for j = 1:n
            dC_da_j = dC_dx_vec(:,j);
            dC_da{j} = reshape(dC_da_j, n*k, n);
            da=zeros(n,m);
            da(j,param_idx)=1;
            term1(:,j) = dC_da{j} * v -  C * da * X ;              
        end
        de_da = Q_e * ( term1 - C * A_obs * dx_da_all{param_idx});

        % 保存该参数的所有偏导
        Ja{param_idx} = de_da;

    end
    J_total = [J0, Ja{1}, Ja{2}];

    % 单位权方差
    sit0_1 = (v' * P_v * v) / (n - m);
    sit0_1 = 1;
    Sigma_e = sit0_1 * J_total * Q_total * J_total';
    row_indices = 1:3:3*n;

    % 创建列索引：1, 4, 7, ..., 3n
    col_indices = 1:3:3*n;

    % 提取子矩阵
    Sigma_e_L = Sigma_e(row_indices, col_indices);
    Sigma_e_a1 = Sigma_e(row_indices+1, col_indices+1);
    diag_newton_all(:, rho_idx) = [diag(Sigma_e_L);diag(Sigma_e_a1)];
    

    
    %% ============== 数值法 (蒙特卡洛) ==============
    fprintf('  计算数值法 (蒙特卡洛)...\n');
    a_sim = zeros(num_monte_carlo, k*n);
    % 构建L和A1的联合协方差矩阵（2n×2n）
    % 从完整的协方差矩阵中提取L和A1的子矩阵
    idx_L = 1:total_vars:n*total_vars;      % L的索引：1, 4, 7, ...
    idx_A1 = 2:total_vars:n*total_vars;     % A1的索引：2, 5, 8, ...

    % 提取L和A1的协方差子矩阵
    Q_L_only = Sigma_global(idx_L, idx_L);
    Q_A1_only = Sigma_global(idx_A1, idx_A1);
    Q_cross_LA1 = Sigma_global(idx_L, idx_A1); % L和A1的互协方差

    % 构建联合协方差矩阵
    Sigma_LA1 = [Q_L_only, Q_cross_LA1;
        Q_cross_LA1', Q_A1_only];

    % 计算Cholesky分解
    R_LA1 = chol(Sigma_LA1);
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
        e_hat = E_hat(A_sim, L_sim, PP, X_sim);
        a_sim(i, :) = e_hat;
    end

    cov_num = cov(a_sim);  % 数值法协方差矩阵
    cov_num_e_L = cov_num(row_indices, col_indices);
    cov_num_e_a1 = cov_num(row_indices+1, col_indices+1);
    diag_numerical_all(:, rho_idx) = [diag(cov_num_e_L);diag(cov_num_e_a1)];
    
    fprintf('  相关系数 rho = %.1f 完成\n\n', rho);
end

%% ============== 第一个图框：对比eL ==============
fprintf('绘制第一个图框：对比eL...\n');

fig_width = 1200;

fig1 = figure('Position', [100, 100, 1200, 350]);
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% 为每个相关系数创建子图
for rho_idx = 1:num_rhos
    subplot(1, 3, rho_idx); % 1行3列
    
    rho = rho_values(rho_idx);
    
    % 提取当前相关系数的结果（eL部分）
    diag_newton_eL = diag_newton_all(1:n, rho_idx);
    diag_numerical_eL = diag_numerical_all(1:n, rho_idx);
    
    % 绘制两种方法的eL对比
    plot(1:n, diag_newton_eL, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Newton');
    hold on;
    plot(1:n, diag_numerical_eL, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'MCS');
    hold off;
    
    xlabel('Obs.No.', 'FontSize', 18);
    ylabel('Diagonal Elements of e_L', 'FontSize', 18);
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
    y_max = max([diag_newton_eL; diag_numerical_eL]) * 1.1;
    y_min = min([diag_newton_eL; diag_numerical_eL]) * 0.9;
    if y_min > 0
        y_min = y_min * 0.9;
    else
        y_min = y_min * 1.1;
    end
    ylim([y_min, y_max]);
end



%% ============== 第二个图框：对比eA1 ==============
fprintf('绘制第二个图框：对比eA1...\n');

fig2 = figure('Position', [100, 100, 1200, 350]);
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% 为每个相关系数创建子图
for rho_idx = 1:num_rhos
    subplot(1, 3, rho_idx); % 1行3列
    
    rho = rho_values(rho_idx);
    
    % 提取当前相关系数的结果（eA1部分）
    diag_newton_eA1 = diag_newton_all(n+1:2*n, rho_idx);
    diag_numerical_eA1 = diag_numerical_all(n+1:2*n, rho_idx);
    
    % 绘制两种方法的eA1对比
    plot(1:n, diag_newton_eA1, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Newton');
    hold on;
    plot(1:n, diag_numerical_eA1, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'MCS');
    hold off;
    
    xlabel('Obs.No.', 'FontSize', 18);
    ylabel('Diagonal Elements of e_{A1}', 'FontSize', 18);
    % 添加子图标题并获取句柄
    h2 = title(sprintf('ρ = %.1f', rho), 'FontSize', 18, 'FontWeight', 'bold');


    % 获取当前坐标轴范围并调整标题位置
    ax2 = gca;
    % 使用归一化坐标调整标题位置
    set(h2, 'Units', 'normalized');
    titlePos2 = get(h2, 'Position');

    
    % 只在第一个子图设置图例
    if rho_idx == 1
        legend('Location', 'northwest', 'FontSize', 16);
    end
    
    grid on;
    
    % 设置y轴范围一致以便比较
    y_max = max([diag_newton_eA1; diag_numerical_eA1]) * 1.1;
    y_min = min([diag_newton_eA1; diag_numerical_eA1]) * 0.9;
    if y_min > 0
        y_min = y_min * 0.9;
    else
        y_min = y_min * 1.1;
    end
    ylim([y_min, y_max]);
end
