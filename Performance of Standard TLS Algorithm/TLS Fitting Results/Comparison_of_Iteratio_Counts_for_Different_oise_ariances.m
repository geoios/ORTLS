%% ============== 系统测试四种方法的迭代次数（复杂问题） ==============
clear all; close all; clc;

% 设置随机数种子以保证可重复性
rng(123);

% 设置问题维度
n = 100;  % 观测数量增加到100
m = 2;    % 参数数量

% 真实参数 - 生成随机数据
A1_true = 10 * rand(n, 1);  % 0-10均匀分布的随机数
A2_true = ones(n, 1);
x_true = [2.5; -1.2];  % 随机真实参数值
L_true = A1_true * x_true(1) + A2_true * x_true(2);

% 噪声方差水平设置（从小到大）
sigma_levels = [1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.5, 1, 2, 5, 10];
n_levels = length(sigma_levels);

% 收敛容差设置
epsilon = 1e-10;
tol = 1e-10;

% 每种噪声水平运行10次以获得统计结果
num_runs = 10;  % 增加到10次实验

% 初始化结果存储
iter_jazaeri = zeros(n_levels, 1);
iter_newton = zeros(n_levels, 1);
iter_wang = zeros(n_levels, 1);
iter_fang = zeros(n_levels, 1);

% 初始化收敛状态存储（1=收敛，0=不收敛）
converge_jazaeri = ones(n_levels, num_runs);
converge_newton = ones(n_levels, num_runs);
converge_wang = ones(n_levels, num_runs);
converge_fang = ones(n_levels, num_runs);

% 初始化每次实验的迭代次数存储
iter_details_jazaeri = zeros(n_levels, num_runs);
iter_details_newton = zeros(n_levels, num_runs);
iter_details_wang = zeros(n_levels, num_runs);
iter_details_fang = zeros(n_levels, num_runs);

fprintf('开始测试不同噪声水平下的迭代次数...\n');
fprintf('问题维度: n=%d, m=%d\n', n, m);
fprintf('噪声水平: '); fprintf('%g ', sigma_levels); fprintf('\n\n');
fprintf('每个噪声水平进行 %d 次实验\n\n', num_runs);

% 主测试循环
for sigma_idx = 1:n_levels
    sigma = sigma_levels(sigma_idx);
    
    % 初始化本次噪声水平的迭代次数
    temp_iter_jazaeri = zeros(num_runs, 1);
    temp_iter_newton = zeros(num_runs, 1);
    temp_iter_wang = zeros(num_runs, 1);
    temp_iter_fang = zeros(num_runs, 1);
    
    fprintf('测试噪声水平: σ² = %g\n', sigma);
    fprintf('   进度: ');
    
    for run = 1:num_runs
        % 每完成一次运行显示进度
        if mod(run, ceil(num_runs/10)) == 0 || run == num_runs
            fprintf('=');
        end
        
        % 设置当前噪声水平
        sigma_x = sigma;         % x的噪声方差
        sigma_L = sigma;         % L的噪声方差
        
        % 生成带噪声的观测值
        L_obs = L_true + sqrt(sigma_L) * randn(n, 1);
        A1_obs = A1_true + sqrt(sigma_x) * randn(n, 1);
        A_obs = [A1_obs, A2_true];
        
        % 设置权阵
        py = (1/sigma_L) * ones(n, 1);  % L的权
        px1 = (1/sigma_x) * ones(n, 1); % A1的权
        px2 = 10^12 * ones(n, 1);       % A2的权（近似无噪声）
        pA = px1;                       % Fang和Wang法使用
        
        % 构造权矩阵P（3n×3n对角阵）
        P = zeros(3*n, 3*n);
        for i = 1:n
            P((i-1)*3+1, (i-1)*3+1) = py(i);     % L的权
            P((i-1)*3+2, (i-1)*3+2) = px1(i);    % A1的权
            P((i-1)*3+3, (i-1)*3+3) = px2(i);    % A2的权
        end
        
        %% ============== Jazaeri法 ==============
        % 初始化权矩阵
        Q_y = diag(1./py);
        Q_A = zeros(n*m, n*m);
        Q_A(1:n, 1:n) = diag(1./px1);           % A1列有误差
        Q_A(n+1:n*m, n+1:n*m) = diag(1./px2);   % A2列有误差（实际上很小）
        
        % 从加权最小二乘开始
        P0 = diag(py);
        x_hat = (A_obs' * P0 * A_obs) \ (A_obs' * P0 * L_obs);
        
        % 开始迭代
        iter = 0;
        max_iter = 1000;  % 最大迭代次数
        converged = 0;    % 收敛标志
        
        for iter = 1:max_iter
            % 估计残差
            e_hat = L_obs - A_obs * x_hat;
            
            % 计算 Q_y_tilde
            x_kron_T = kron(x_hat', eye(n));
            x_kron = kron(x_hat, eye(n));
            Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
            
            % 检查矩阵条件数，避免病态问题
            if rcond(Q_y_tilde) < eps
                Q_y_tilde_inv = pinv(Q_y_tilde);
            else
                Q_y_tilde_inv = inv(Q_y_tilde);
            end
            
            % 计算 E_A_hat
            vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
            E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];
            
            % 预测值
            A_tilde = A_obs - E_A_hat;
            y_tilde = L_obs - E_A_hat * x_hat;
            
            % 更新参数
            x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ ...
                       (A_tilde' * Q_y_tilde_inv * y_tilde);
            
            % 收敛检查
            delta = norm(x_hat_new - x_hat);
            x_hat = x_hat_new;
            
            if delta < epsilon
                converged = 1;
                break;
            end
        end
        
        % 记录结果
        temp_iter_jazaeri(run) = iter;
        iter_details_jazaeri(sigma_idx, run) = iter;
        converge_jazaeri(sigma_idx, run) = converged;
        
        %% ============== Newton法 ==============
        [~, ~, iter] = TLS_XG_newton3(A_obs, L_obs, P);
        temp_iter_newton(run) = iter;
        iter_details_newton(sigma_idx, run) = iter;
        converge_newton(sigma_idx, run) = (iter < max_iter);
        
        %% ============== Fang法 ==============
        % 注意：Fang法需要将观测数据组织为x和y
        x_data = A1_obs;  % A1列
        y_data = L_obs;   % L观测值
        
        % 调用Fang法
        try
            [~, cf] = Fang(x_data, y_data, py, pA);
            temp_iter_fang(run) = cf;
            iter_details_fang(sigma_idx, run) = cf;
            converge_fang(sigma_idx, run) = (cf < 1000);
        catch
            warning('Fang法在噪声水平%g的第%d次运行中失败', sigma, run);
            temp_iter_fang(run) = max_iter;
            iter_details_fang(sigma_idx, run) = max_iter;
            converge_fang(sigma_idx, run) = 0;
        end
        
        %% ============== Wang法 ==============
        % 注意：Wang法只处理A1有误差的情况
        try
            [~, ci] = Wang(n, A1_obs, L_obs, py, pA);
            temp_iter_wang(run) = ci;
            iter_details_wang(sigma_idx, run) = ci;
            converge_wang(sigma_idx, run) = (ci < 1000);
        catch
            warning('Wang法在噪声水平%g的第%d次运行中失败', sigma, run);
            temp_iter_wang(run) = max_iter;
            iter_details_wang(sigma_idx, run) = max_iter;
            converge_wang(sigma_idx, run) = 0;
        end
    end
    fprintf('\n');
    
    % 计算平均迭代次数并取整（舍去小数）
    iter_jazaeri(sigma_idx) =mean(temp_iter_jazaeri);
    iter_newton(sigma_idx) = mean(temp_iter_newton);
    iter_wang(sigma_idx) = mean(temp_iter_wang);
    iter_fang(sigma_idx) = mean(temp_iter_fang);
    
    % 计算收敛率
    jazaeri_converge_rate = sum(converge_jazaeri(sigma_idx, :)) / num_runs * 100;
    newton_converge_rate = sum(converge_newton(sigma_idx, :)) / num_runs * 100;
    wang_converge_rate = sum(converge_wang(sigma_idx, :)) / num_runs * 100;
    fang_converge_rate = sum(converge_fang(sigma_idx, :)) / num_runs * 100;
    
    fprintf('   迭代次数: Jazaeri=%d, Newton=%d, Wang=%d, Fang=%d\n', ...
            iter_jazaeri(sigma_idx), iter_newton(sigma_idx), ...
            iter_wang(sigma_idx), iter_fang(sigma_idx));
    fprintf('   收敛率: Jazaeri=%.1f%%, Newton=%.1f%%, Wang=%.1f%%, Fang=%.1f%%\n\n', ...
            jazaeri_converge_rate, newton_converge_rate, ...
            wang_converge_rate, fang_converge_rate);
end

%% ============== 结果显示 ==============
fprintf('\n============== 平均迭代次数结果汇总 ==============\n');
fprintf('Noise Variance\tJazaeri\t\tNewton\t\tWang\t\tFang\n');
fprintf('----------------------------------------------------\n');
for i = 1:n_levels
    fprintf('%12.2e\t%8d\t%8d\t%8d\t%8d\n', ...
            sigma_levels(i), iter_jazaeri(i), iter_newton(i), ...
            iter_wang(i), iter_fang(i));
end

%% ============== 详细收敛情况分析 ==============
fprintf('\n============== 详细收敛情况分析 ==============\n');
fprintf('每个噪声水平下10次实验的收敛情况:\n');
fprintf('(C=收敛, NC=不收敛, 数字=迭代次数)\n\n');

for sigma_idx = 1:n_levels
    sigma = sigma_levels(sigma_idx);
    fprintf('噪声水平 σ² = %g:\n', sigma);
    fprintf('Experiment  Jazaeri   Newton    Wang      Fang\n');
    fprintf('----------  -------   ------    ----      ----\n');
    
    for run = 1:num_runs
        % Jazaeri法
        jazaeri_status = 'C';
        if converge_jazaeri(sigma_idx, run) == 0
            jazaeri_status = 'NC';
        end
        
        % Newton法
        newton_status = 'C';
        if converge_newton(sigma_idx, run) == 0
            newton_status = 'NC';
        end
        
        % Wang法
        wang_status = 'C';
        if converge_wang(sigma_idx, run) == 0
            wang_status = 'NC';
        end
        
        % Fang法
        fang_status = 'C';
        if converge_fang(sigma_idx, run) == 0
            fang_status = 'NC';
        end
        
        fprintf('  %2d        %s(%3d)   %s(%3d)   %s(%3d)   %s(%3d)\n', ...
                run, ...
                jazaeri_status, iter_details_jazaeri(sigma_idx, run), ...
                newton_status, iter_details_newton(sigma_idx, run), ...
                wang_status, iter_details_wang(sigma_idx, run), ...
                fang_status, iter_details_fang(sigma_idx, run));
    end
    
    % 统计本噪声水平的收敛率
    jazaeri_converge_rate = sum(converge_jazaeri(sigma_idx, :)) / num_runs * 100;
    newton_converge_rate = sum(converge_newton(sigma_idx, :)) / num_runs * 100;
    wang_converge_rate = sum(converge_wang(sigma_idx, :)) / num_runs * 100;
    fang_converge_rate = sum(converge_fang(sigma_idx, :)) / num_runs * 100;
    
    fprintf('收敛率:      %.1f%%        %.1f%%        %.1f%%        %.1f%%\n\n', ...
            jazaeri_converge_rate, newton_converge_rate, ...
            wang_converge_rate, fang_converge_rate);
end

%% ============== 总体收敛率统计 ==============
fprintf('\n============== 总体收敛率统计 ==============\n');
fprintf('Method\t\tConvergence Rate\tSuccess Count\tTotal Runs\n');
fprintf('-----------------------------------------------------------\n');

overall_jazaeri = mean(converge_jazaeri(:)) * 100;
overall_newton = mean(converge_newton(:)) * 100;
overall_wang = mean(converge_wang(:)) * 100;
overall_fang = mean(converge_fang(:)) * 100;

fprintf('Jazaeri\t\t%.1f%%\t\t\t%d/%d\n', ...
        overall_jazaeri, sum(converge_jazaeri(:)), numel(converge_jazaeri));
fprintf('Newton\t\t%.1f%%\t\t\t%d/%d\n', ...
        overall_newton, sum(converge_newton(:)), numel(converge_newton));
fprintf('Wang\t\t%.1f%%\t\t\t%d/%d\n', ...
        overall_wang, sum(converge_wang(:)), numel(converge_wang));
fprintf('Fang\t\t%.1f%%\t\t\t%d/%d\n', ...
        overall_fang, sum(converge_fang(:)), numel(converge_fang));

%% ============== 可视化结果 ==============
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% 设置字体为新罗马字体，16号
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');
set(groot, 'DefaultLegendFontName', 'Times New Roman');

% 只绘制一个图：迭代次数对比
plot(sigma_levels, iter_jazaeri, 'ro-', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Jazaeri');
hold on;
plot(sigma_levels, iter_newton, 'bs-', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Newton');
plot(sigma_levels, iter_wang, 'g^-', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Wang');
plot(sigma_levels, iter_fang, 'md-', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Fang');
grid on;
box on;

% 设置坐标轴标签和标题（英文）
xlabel('Noise Variance \sigma^2', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Average Iteration Count (Rounded)', 'FontSize', 16, 'FontWeight', 'bold');
title('Comparison of Iteration Counts for Different Methods', 'FontSize', 18, 'FontWeight', 'bold');

% 设置图例
legend('Location', 'best', 'FontSize', 14);

% 设置坐标轴属性
set(gca, 'FontSize', 16, 'XScale', 'log', 'XMinorTick', 'on', 'YMinorTick', 'on');
set(gca, 'LineWidth', 1.5, 'GridLineStyle', '--', 'GridAlpha', 0.3);

% 添加网格
grid on;
grid minor;

% 调整坐标轴范围（如有需要）
xlim([min(sigma_levels)/10, max(sigma_levels)*10]);



%% ============== 定义Fang法函数 ==============
function [X, cf] = Fang(x, y, py, px)
    [n, m] = size(x);
    
    % 构建权矩阵
    Py = diag(py);
    PA = diag(px);
    QA = inv(PA);
    Qy = inv(Py);
    Qxy = zeros(n, n);
    Q = [QA, zeros(n), Qxy; zeros(n, 3*n); Qxy, zeros(n), Qy];
    
    % 最小二乘初值
    A = [x, ones(n, 1)];
    X0 = pinv(A' * Py * A) * A' * Py * y;
    
    cf = 0;
    cita = 1;
    max_iter = 1000;
    
    while cita > 1e-10 && cf < max_iter
        cf = cf + 1;
        I = eye(n);
        B = [kron(X0', I), -I];
        
        % 计算残差
        r = pinv(B * Q * B') * (y - A * X0);
        
        % 计算误差估计
        vA = [QA, zeros(n), Qxy; zeros(n, 3*n)] * B' * r;
        
        % 更新参数
        X = pinv(A' * pinv(B * Q * B') * A) * ...
            (kron(eye(2), r') * vA + A' * pinv(B * Q * B') * y);
        
        % 检查收敛
        cita = norm(X - X0);
        X0 = X;
    end
end

%% ============== 定义Wang法函数 ==============
function [X, ci] = Wang(n, A1, L, py, pA)
    % 构建权矩阵
    Py = diag(py);
    PA = diag(pA);
    QA = inv(PA);
    Qy = inv(Py);
    Qxy = zeros(n, n);
    Q = [QA, Qxy; Qxy, Qy];
    
    % 最小二乘初值
    A = [A1, ones(n, 1)];
    X0 = pinv(A' * Py * A) * A' * Py * L;
    
    % 初始化变量
    I = eye(n);
    C1 = eye(2*n);
    h = [zeros(n, 1); ones(n, 1)];
    b = [1; 0];
    B = kron(b, I);
    
    ci = 0;
    cita = 1;
    max_iter = 1000;
    
    while cita > 1e-10 && ci < max_iter
        ci = ci + 1;
        C2 = [-kron(X0', I) * B, I] * C1;
        
        % 计算拉格朗日乘子
        lambda = (C2 * Q * C2') \ (L - kron(X0', I) * (h + B * A1));
        
        % 计算误差
        e = -Q * C2' * lambda;
        eA = e(1:n);
        EA = reshape(B * eA, n, 2);
        
        % 更新系数矩阵和观测值
        A0 = A + EA;
        y1 = L + EA * X0;
        
        % 更新参数
        X = pinv(A0' * pinv(C2 * Q * C2') * A0) * A0' * pinv(C2 * Q * C2') * y1;
        
        % 检查收敛
        cita = norm(X - X0);
        X0 = X;
    end
end

%% ============== 定义Newton法函数 ==============
function [x_tls, e_hat, iter] = TLS_XG_newton3(A_obs, L_obs, P)
    tol = 1e-10;
    [n, m] = size(A_obs);
    k = m + 1;
    
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
    max_iter = 1000;
    
    while dx_norm > tol && iter < max_iter
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
        A_corr2 = A_obs + 2 * e_A;
        
        % 7. 计算梯度 F
        F = A_corr' * vp;
        
        % 8. 计算Hessian矩阵
        P_v_e_A = P_v * e_A;
        P_v_A_obs = P_v * A_obs;
        
        H1 = -A_corr' * P_v * A_corr2;
        
        % H5矩阵
        H5 = zeros(m, m);
        for i = 1:m
            v_tk = zeros(k*n, 1);
            v_tk(i+1:k:(k*n)) = vp;
            
            E_bk = kron(P_v_e_A(:, i)', b);
            A_bk = kron(P_v_A_obs(:, i)', b);
            
            de_dxi = Q_e * (v_tk - 2 * E_bk' - A_bk');
            
            de_dxi_reshaped = reshape(de_dxi, k, n);
            dET_dxi = de_dxi_reshaped(2:k, :);
            
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