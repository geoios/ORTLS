clear all; clc; close all;

% 实验参数设置
n_values = [10, 100:100:1000]; % 观测值数量数组
m = 2;                          % 参数个数
k = m + 1;                      % 每个观测点的变量数(L + A1 + A2)
num_n = length(n_values);       % 观测维度数量

% 设置噪声方差
sigma_x = 0.001;                % x的噪声方差
sigma_L = 0.001;                % L的噪声方差

% 初始化存储运行时间的数组
time_newton = zeros(num_n, 1);      % 牛顿法
time_schaffrin = zeros(num_n, 1);   % Schaffrin法
time_wang = zeros(num_n, 1);        % Wang法
time_xu = zeros(num_n, 1);          % Xu法
time_fang = zeros(num_n, 1);        % Fang法
time_jazaeri = zeros(num_n, 1);     % Jazaeri法

% 循环处理不同观测维度
for n_idx = 1:num_n
    n = n_values(n_idx);
    fprintf('=== 实验: 观测维度 n = %d ===\n', n);
    
    % 生成真实数据
    A1_true = 10 * rand(n, 1);      % [0,10]区间均匀分布
    A2_true = ones(n, 1);          % 常数项
    L_true = 2 * A1_true - A2_true;
    
    % 生成带噪声的观测数据
    A1_obs = A1_true + sigma_x * randn(n, 1);
    L_obs = L_true + sigma_L * randn(n, 1);
    A2_obs = A2_true;  % A2没有误差
    A_obs = [A1_obs, A2_obs];
    
    % 设置权阵
    py = (1/sigma_L) * ones(n, 1);  % L的权
    Py = diag(py);
    px1 = (1/sigma_x) * ones(n, 1); % A1的权
    px2 = 10^12 * ones(n, 1);       % A2的权（近似无噪声）
    pA = [px1, px2];
    P = [py, px1, px2]';            % 3 x n
    PP = diag(P(:));
    
    %% ============== 牛顿法计时 ==============
    fprintf('  计算牛顿法...\n');
    tic;
    X_newton = TLS_XG_newton3(A_obs, L_obs, PP);
    time_newton(n_idx) = toc;
    
    %% ============== Schaffrin法计时 ==============
    fprintf('  计算Schaffrin法...\n');
    % 构建观测向量的协因数矩阵 Q_y
    p_y = py;  % L的权值
    p_A = px1; % A第一列的权值
    Q_y = diag(1 ./ p_y);  % n×n 对角矩阵
    
    % 构建设计矩阵A的协因数矩阵 Q_A
    p_A_vec = [1 ./ p_A; 1e-15 * ones(n, 1)];  % 前n个是x的协因数，后n个是常数项的协因数
    Q_A = diag(p_A_vec);  % 2n×2n 对角矩阵
    
    % 交叉协方差矩阵（假设为零）
    Q_yA = zeros(n, 2*n);
    
    max_iter = 100;
    tol = 1e-10;
    
    tic;
    [X_schaffrin, ~, ~, ~, ~, ~] = schaffrin_line_fit_with_residuals(A_obs, L_obs, Q_y, Q_A, Q_yA, max_iter, tol);
    time_schaffrin(n_idx) = toc;
    
    %% ============== Wang法计时 ==============
    fprintf('  计算Wang法...\n');
    tic;
    [X_wang, ~] = Wang(n, A1_obs, L_obs, py, px1);
    time_wang(n_idx) = toc;
    
    %% ============== Xu法计时 ==============
    fprintf('  计算Xu法...\n');
    tic;
    [X_xu, ~, ~, ~, ~] = Xu_m(A_obs, L_obs, P);
    time_xu(n_idx) = toc;
    
    %% ============== Fang法计时 ==============
    fprintf('  计算Fang法...\n');
    tic;
    [X_fang, ~] = Fang(A1_obs, L_obs, py, px1);
    time_fang(n_idx) = toc;
    
    %% ============== Jazaeri法计时 ==============
    fprintf('  计算Jazaeri法...\n');
    tic;
    % 初始化权矩阵
    A = A_obs;
    L = L_obs;
    Q_y = diag(1./py);
    Q_A = zeros(m*n, m*n);
    Q_A(1:n, 1:n) = diag(1./px1);
    Q_A(n+1:m*n, n+1:m*n) = zeros(n);
    
    % 从加权最小二乘开始
    P0 = diag(py);
    x_jaz = inv(A' * P0 * A) * (A' * P0 * L);
    epsilon = 1e-10;
    
    % 开始迭代
    for iter = 1:20
        % 估计残差
        e_hat = L - A * x_jaz;
        
        % 计算 Q_y_tilde = Q_y + (x^T ⊗ I_m) Q_A (x ⊗ I_m)
        x_kron_T = kron(x_jaz', eye(n));  % m × mn
        x_kron = kron(x_jaz, eye(n));      % mn × m
        Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
        Q_y_tilde_inv = inv(Q_y_tilde);
        
        % 计算 E_A_hat = -vec^{-1}(Q_A (x ⊗ I_m) Q_y_tilde^{-1} e_hat)
        vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
        E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];
        
        % 预测值
        A_tilde = A - E_A_hat;
        y_tilde = L - E_A_hat * x_jaz;
        
        % 更新参数
        x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * y_tilde);
        
        % 收敛检查
        delta = norm(x_hat_new - x_jaz);
        
        x_jaz = x_hat_new;
        
        if delta < epsilon
            break;
        end
    end
    time_jazaeri(n_idx) = toc;
    
    fprintf('  观测维度 n = %d 完成\n\n', n);
end

%% ============== 绘制运行时间对比图 ==============
fprintf('绘制运行时间对比图...\n');

% 创建图形窗口
figure('Position', [100, 100, 1200, 800]);

% 设置图形参数
markerSize = 10;
lineWidth = 2;
fontSize = 20;
tickFontSize = 18;  % 新增：坐标刻度字号为8

% 使用单纵轴绘制所有曲线
hold on;

% 牛顿法 - 蓝色
plot(n_values, time_newton, 'b-o', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'Newton Method');

% 其他方法 - 红色，使用不同标记区分
plot(n_values, time_schaffrin, 'r-s', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'Schaffrin Method');
plot(n_values, time_wang, 'r-^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'Wang Method');
plot(n_values, time_xu, 'r-d', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'Xu Method');
plot(n_values, time_fang, 'r-v', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'Fang Method');
plot(n_values, time_jazaeri, 'r-*', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'Jazaeri Method');

% 设置图形属性
xlabel('Observation Dimension (n)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Runtime / s', 'FontSize', 20, 'FontWeight', 'bold');

% 设置坐标刻度字号为8
set(gca, 'FontSize', tickFontSize);

% 添加图例
legend('Location', 'northwest', 'FontSize', fontSize-2, 'NumColumns', 2);

% 添加网格
grid on;
box on;

% 设置y轴范围，确保所有数据可见
ylim([0, max([time_newton; time_schaffrin; time_wang; time_xu; time_fang; time_jazaeri]) * 1.1]);



%% ============== 显示运行时间统计信息 ==============
fprintf('\n运行时间统计信息:\n');
fprintf('观测维度: %s\n', mat2str(n_values));
fprintf('牛顿法运行时间: %s 秒\n', mat2str(time_newton'));
fprintf('Schaffrin法运行时间: %s 秒\n', mat2str(time_schaffrin'));
fprintf('Wang法运行时间: %s 秒\n', mat2str(time_wang'));
fprintf('Xu法运行时间: %s 秒\n', mat2str(time_xu'));
fprintf('Fang法运行时间: %s 秒\n', mat2str(time_fang'));
fprintf('Jazaeri法运行时间: %s 秒\n', mat2str(time_jazaeri'));

fprintf('\n所有实验完成！\n');