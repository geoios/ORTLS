
clear all; clc;

% 实验参数设置
num_experiments = 100;  % 实验次数
n = 10;                % 观测值数量
m = 2;                  % 参数个数
k = m+1;
num_monte_carlo = 1000; % 蒙特卡洛模拟次数

% 初始化结果存储数组
newton_errors = zeros(num_experiments, 1);
wang_errors = zeros(num_experiments, 1);
xu_errors = zeros(num_experiments, 1);
jaz_errors = zeros(num_experiments, 1);

newton_errors_sum = zeros(num_experiments, 1);
wang_errors_sum = zeros(num_experiments, 1);
xu_errors_sum = zeros(num_experiments, 1);
jaz_errors_sum = zeros(num_experiments, 1);

newton_closer_count = 0;
wang_closer_count = 0;
xu_closer_count = 0;
jaz_closer_count = 0;

% 生成真实数据
A1_true = (0:9)';
A2_true = ones(n, 1);
% A2_true = (15:24)';
L_true = 2*A1_true - A2_true;
A_true = [A1_true, A2_true];

% 设置噪声方差
sigma_x = 1;         % x的噪声方差
sigma_L = 1;         % L的噪声方差

% 设置权阵
py = (1/sigma_L) * ones(n, 1);  % L的权
Py = diag(py);
px1 = (1/sigma_x) * ones(n, 1); % A1的权
px2 = 10^12 * ones(n, 1);        % A2的权（近似无噪声）
pA = [px1, px2];
P = [py, px1, px2]';       % 4 x n
PP = diag(P(:));

for exp_idx = 1:num_experiments
    fprintf('实验 %d/%d 进行中...\n', exp_idx, num_experiments);



    % 生成带噪声的x和y
    L_obs = L_true + sqrt(sigma_L) * randn(n, 1);
    A1 = A1_true + sqrt(sigma_x) * randn(n, m-1);
%     A2 = A2_true + sqrt(sigma_x) * randn(n, m-1);
    A_obs = [A1, A2_true];
    %% ============== newton法 ==============
    X = TLS_XG_newton3(A_obs, L_obs, PP);
    v = L_obs - A_obs * X;
    % 计算误差传播
    Q_e = inv(PP);
    Q_total = [P(1,:),P(2,:),P(3,:)];
    Q_total = inv(diag(Q_total));
    [H, e_A_newton, B, e_newton] = Hessian(A_obs, L_obs, PP, X);
    newton_errors_sum(exp_idx) = mean(e_newton);
    e_A=-e_A_newton;

    % 计算J矩阵[J0,J1,J2]2*30
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
        %         Sigma_a_i = sit0_1 * inv(Pa{param_idx});
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
            dF_da_i(:, obs_idx) = dF_da_single;
            if rcond(H) < eps
                dx_da_i = -pinv(H) * dF_da_i;  % m×n
            else
                dx_da_i = -H \ dF_da_i;        % 更高效，等价于 -inv(H) * de_da_i
            end
            da=zeros(n,m);
            da(obs_idx,param_idx)=1;
            term1(:,obs_idx) = da * X ;
        end
        Ja{param_idx} = dx_da_i;

    end
    J_total = [J0, Ja{1}, Ja{2}];

    % 单位权方差
    sit0_1 = (e_newton' * PP * e_newton) / (n - m);
    mm_newton=sqrt(sit0_1);
    Dx_newton = sit0_1 *  J_total * Q_total * J_total';

    %% ============== Wang法 ==============

    [X, sigma2, Qxx, ci] = Wang_sigma(n, A_obs, L_obs, py, pA(:,1), 0);
%     wang_errors_sum(exp_idx) = mean(e_wang);
    Dx_Wang = sigma2*Qxx;

    %% ============== Xu法 ==============
    [beta_hat, a_hat, iter, conv_flag, S0, A0, PA] = Xu_m(A_obs, L_obs, P);

    % 修正：使用观测值（而非真实值）计算单位权方差
    % 注意：实际应用中无法获得真实值A1_true，只能使用观测值a
    a_hat0 = A_obs(:,1);
    r_a = a_hat0(:) - a_hat;            % 随机元素的残差（使用观测值a）
    r_y = L_obs - A0 * beta_hat;       % 观测向量的残差

    xu_errors_sum(exp_idx) = mean(r_a + r_y);
    % 计算单位权方差
    sigma2_hat = (r_a' * PA * r_a + r_y' * Py * r_y) / (n - m);
    mm = sqrt(sigma2_hat);
    aa_hat=[a_hat,ones(n,1)];
    % 计算法方程矩阵
    N11 = A0' * Py * A0;
    N12 = A0' * Py * S0;
    N21 = S0' * Py * A0;
    N22 = S0' * Py * S0 + PA;

    % 计算参数β的方差-协方差矩阵
    Dx_xu = inv(N11 - N12 * inv(N22) * N21) * sigma2_hat;

    %% ============== jazaeri法 ==============
    % 初始化权矩阵
    Q_y = diag(1./py);
    Q_A = zeros(m*n, m*n);
    % Q_A结构：vec(A)按列堆叠，所以前m个元素是x列，后m个是常数列
    Q_A(1:n, 1:n) = diag(1./px1);           % x列有误差
    Q_A(n+1:m*n, n+1:m*n) = zeros(n); % 常数列无误差

    % 从加权最小二乘开始
    P0 = diag(py);
    x_hat = inv(A_obs' * P0 * A_obs) * (A_obs' * P0 * L_obs);
    epsilon = 1e-10;
    for iter = 1:20
        % 估计残差
        e_hat = L_obs - A_obs * x_hat;

        % 计算 Q_y_tilde = Q_y + (x^T ⊗ I_m) Q_A (x ⊗ I_m)
        x_kron_T = kron(x_hat', eye(n));  % m × mn
        x_kron = kron(x_hat, eye(n));      % mn × m
        Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
        Q_y_tilde_inv = inv(Q_y_tilde);

        % 计算 E_A_hat = -vec^{-1}(Q_A (x ⊗ I_m) Q_y_tilde^{-1} e_hat)
        vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
        E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];

        % 预测值
        A_tilde = A_obs - E_A_hat;
        y_tilde = L_obs - E_A_hat * x_hat;

        % 更新参数
        x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * y_tilde);

        % 收敛检查
        delta = norm(x_hat_new - x_hat);

        x_hat = x_hat_new;

        if delta < epsilon
            break;
        end
    end

    %% 计算单位权方差和参数协方差矩阵
    final_residuals = L_obs - A_obs * x_hat;

    % 重新计算Q_y_tilde用于最终参数
    x_kron_T = kron(x_hat', eye(n));
    x_kron = kron(x_hat, eye(n));
    Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
    Q_y_tilde_inv = inv(Q_y_tilde);

    sigma_0_sq = (final_residuals' * Q_y_tilde_inv * final_residuals) / (n - m);
    sigma_0 = sqrt(sigma_0_sq);

    % 参数协方差矩阵计算
    vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * final_residuals;
    E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];
    A_tilde = A_obs - E_A_hat;
    Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);

    Dx_Jaz = sigma_0_sq * Q_x;

    %% ============== 数值法 (LS) ==============

    for i = 1:num_monte_carlo

        L_noise_sim = sqrt(sigma_L) * randn(n, 1);
        A1_noise_sim = sqrt(sigma_x) * randn(n, 1);
%         A2_noise_sim = sqrt(sigma_x) * randn(n, 1);
        % 生成模拟观测值
        L_sim = L_true + L_noise_sim;
        A1_sim = A1_true + A1_noise_sim;
        A2_sim = A2_true;  % A2没有误差
        A_sim = [A1_sim, A2_sim];
        X_sim = TLS_XG_newton3(A_sim, L_sim, PP);
        a_sim(i, :) = X_sim;
    end

    cov_num = cov(a_sim);  % 数值法协方差矩阵

   %% ============== 结果比较 ==============
    % 计算Frobenius范数差异
    newton_error = norm(Dx_newton, 'fro')-norm( cov_num, 'fro');
    wang_error = norm(Dx_Wang)-norm(cov_num, 'fro');
    xu_error = norm(Dx_xu, 'fro')-norm( cov_num, 'fro');
    jaz_error = norm(Dx_Jaz, 'fro')-norm(cov_num, 'fro');

    % 存储误差
    newton_errors(exp_idx) = newton_error;
    wang_errors(exp_idx) = wang_error;
    xu_errors(exp_idx) = xu_error;
    jaz_errors(exp_idx) = jaz_error;

    % 判断哪种方法更接近数值基准
    [min_error, min_idx] = min([newton_error, wang_error, xu_error, jaz_error]);

    if min_idx == 1
        newton_closer_count = newton_closer_count + 1;
    elseif min_idx == 2
        wang_closer_count = wang_closer_count + 1;
    elseif min_idx == 3
        xu_closer_count = xu_closer_count + 1;
    else
        jaz_closer_count = jaz_closer_count + 1;
    end
end

%% ============== 结果分析 ==============
fprintf('\n===== 实验结果统计 =====\n');
fprintf('总实验次数: %d\n', num_experiments);
fprintf('牛顿法更接近的次数: %d (%.2f%%)\n', newton_closer_count, 100*newton_closer_count/num_experiments);
fprintf('Wang法更接近的次数: %d (%.2f%%)\n', wang_closer_count, 100*wang_closer_count/num_experiments);
fprintf('Xu法更接近的次数: %d (%.2f%%)\n', xu_closer_count, 100*xu_closer_count/num_experiments);
fprintf('Jazaeri法更接近的次数: %d (%.2f%%)\n', jaz_closer_count, 100*jaz_closer_count/num_experiments);

% 计算每个点处四种方法的平均值
pointwise_avg_errors = (newton_errors + wang_errors + xu_errors + jaz_errors) / 4;
% pointwise_avg_errors = (newton_errors  + jaz_errors) / 2;

% 将每种方法在每个点处的范数减去该点的平均值
newton_errors_pointwise_centered = newton_errors - pointwise_avg_errors;
wang_errors_pointwise_centered = wang_errors - pointwise_avg_errors;
xu_errors_pointwise_centered = xu_errors - pointwise_avg_errors;
jaz_errors_pointwise_centered = jaz_errors - pointwise_avg_errors;

fprintf('\n每个点处中心化误差统计:\n');
fprintf('牛顿法中心化误差范围: [%.4e, %.4e]\n', min(newton_errors_pointwise_centered), max(newton_errors_pointwise_centered));
fprintf('Wang法中心化误差范围: [%.4e, %.4e]\n', min(wang_errors_pointwise_centered), max(wang_errors_pointwise_centered));
fprintf('Xu法中心化误差范围: [%.4e, %.4e]\n', min(xu_errors_pointwise_centered), max(xu_errors_pointwise_centered));
fprintf('Jazaeri法中心化误差范围: [%.4e, %.4e]\n', min(jaz_errors_pointwise_centered), max(jaz_errors_pointwise_centered));

% 计算每种方法的原始范数均值（用于参考）
avg_newton_error = mean(newton_errors);
avg_wang_error = mean(wang_errors);
avg_xu_error = mean(xu_errors);
avg_jaz_error = mean(jaz_errors);

fprintf('\n原始范数均值（参考）:\n');
fprintf('牛顿法: %.4e\n', avg_newton_error);
fprintf('Wang法: %.4e\n', avg_wang_error);
fprintf('Xu法: %.4e\n', avg_xu_error);
fprintf('Jazaeri法: %.4e\n', avg_jaz_error);

%% ============== Plotting ==============
% 第一个图：四种方法原始误差的100次变化值
figure('Position', [100, 100, 800, 400]);
hold on;
plot(1:num_experiments, newton_errors, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(1:num_experiments, wang_errors, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(1:num_experiments, xu_errors, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(1:num_experiments, jaz_errors, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 4);

xlabel('Experiment No.', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(' Δ*', 'FontSize', 12, 'FontName', 'Times New Roman');
legend('Newton', 'Wang', 'Xu', 'Jazaeri', 'Location', 'northwest', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gcf, 'Color', 'w');
