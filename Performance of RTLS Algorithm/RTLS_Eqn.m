%% 抗差总体最小二乘（Robust Total Least-Squares, RTLS）估计器
% 基于论文：Robust Total Least-Squares Estimator for Seafloor Geodetic Control Point Positioning
% 实现RTLS_Eqn算法
%
% ========== RTLS算法流程说明 ==========
% 
% 核心思想：通过两次标准化总体残差，然后基于标准化残差计算R矩阵进行抗差
%
% 算法流程：
% 1. 计算总体残差（total residuals）e_t（公式23）
%    e_t = V - A^T δζ，其中V是观测残差
%
% 2. 计算总体残差的协因数矩阵 Q_e（公式24）
%    Q_e = Q_c - A^T (A Q_c^(-1) A^T)^(-1) A^T
%
% 3. 第一次标准化：计算单位权后验标准差 δ_0（公式26）
%    ratio_i = |e_ti| / sqrt(Q_e_ti)  （对每个观测方程的初步标准化）
%    δ_0 = 1.4286 * median(ratio_i)   （取中位数，提高抗差性）
%
% 4. 第二次标准化：计算标准化总残差 v_t（公式25）
%    v_ti = |e_ti| / (δ_0 * sqrt(Q_e_ti))
%    这是最终的标准化残差，用于判断观测的异常程度
%
% 5. 基于标准化总残差计算IGGIII权函数得到R矩阵（公式27）
%    - 保守区（|v_ti| ≤ k0）：R_i = 1，保持原权重
%    - 关注区（k0 < |v_ti| ≤ k1）：R_i < 1，降低权重
%    - 排除区（|v_ti| > k1）：R_i ≈ 0，近似排除
%
% 6. 使用R矩阵进行抗差：计算等价协因数矩阵（公式28）
%    Q̄_c = √(R Q_c R)
%    通过调整协因数矩阵，降低异常值的影响
%
% 7. 使用等价协因数矩阵进行下一次TLS迭代，重复步骤1-6直到收敛
%
% ======================================

function [x_est, sigma0_sq, Q_x, iter_info] = RTLS_Eqn(A, y, Q_c, x0, options)
% RTLS_Eqn - 抗差总体最小二乘估计器
%
% 输入参数:
%   A       - 观测方程系数矩阵 (m x n)
%   y       - 观测向量 (m x 1)
%   Q_c     - 观测值协因数矩阵 (m x m)
%   x0      - 初始参数估计值 (n x 1)，可选
%   options - 结构体，包含以下可选参数:
%             .k0        - IGGIII权函数阈值k0 (默认: 2.5)
%             .k1        - IGGIII权函数阈值k1 (默认: 6.0)
%             .max_iter  - 最大迭代次数 (默认: 100)
%             .tol       - 收敛容差 (默认: 1e-6)
%             .max_inner_iter - 内部循环最大迭代次数 (默认: 50)
%
% 输出参数:
%   x_est      - 估计的参数向量 (n x 1)
%   sigma0_sq  - 单位权方差因子后验估计
%   Q_x        - 参数估计的协因数矩阵
%   iter_info  - 迭代信息结构体

%% 参数设置和初始化
if nargin < 4 || isempty(x0)
    % 初始LS估计
    P_c = inv(Q_c);
    x0 = (A' * P_c * A) \ (A' * P_c * y);
end

if nargin < 5
    options = struct();
end

% 默认参数
k0 = getfield(options, 'k0', 1.5);
k1 = getfield(options, 'k1', 2.5);
max_iter = getfield(options, 'max_iter', 20);
tol = getfield(options, 'tol', 1e-2);
max_inner_iter = getfield(options, 'max_inner_iter', 5);

[m, n] = size(A);
x = x0;

% 迭代信息
iter_info = struct();
iter_info.outer_iter = 0;
iter_info.inner_iter = [];
iter_info.converged = false;

%% 第一步：初始TLS估计
% 构造初始矩阵 A^0, G^0, I^0
% 忽略GNSS定位误差，设 e_x^0 = 0
A_current = A;
x_prev = x;

% 计算初始TLS解 x^1
[x, ~, tls_iter_initial] = TLS_solution(A_current, y, Q_c, x);
tls_iter_total = tls_iter_initial;  % 累积TLS迭代次数

iter_info.outer_iter = 1;

%% 外部循环
for j = 1:max_iter
    iter_info.outer_iter = j;
    
    % 显示进度（每5次迭代显示一次，或第一次和最后一次）
    if j == 1 || mod(j, 5) == 0 || j == max_iter
        fprintf('  外部迭代 %d/%d...\n', j, max_iter);
    end
    
    % 内部循环：重加权操作
    inner_iter_count = 0;
    Q_c_current = Q_c;
    
    for inner_iter = 1:max_inner_iter
        inner_iter_count = inner_iter;
        
        % 显示内部循环进度（仅第一次外部迭代时显示，避免输出过多）
        if j == 1 && (inner_iter == 1 || mod(inner_iter, 5) == 0 || inner_iter == max_inner_iter)
            fprintf('    内部迭代 %d/%d...\n', inner_iter, max_inner_iter);
        end
        
        % 步骤b: 计算总体残差（total residuals）和协因数矩阵（公式23和24）
        % 公式23: ~e_t^(j+1) = V - A^T δζ^(j+1) （总体残差）
        % 公式24: Q_e^(j+1) = Q_c - A^T (A Q_c^(-1) A^T)^(-1) A^T
        [e_t, Q_e] = compute_total_residuals(A_current, y, x, Q_c_current);
        
        % 步骤c: 计算单位权后验标准差
        % 只在第一次外部迭代的第一次内部迭代时使用中位数法（公式26）
        % 后续使用传统的加权残差平方和形式
        use_median = (j == 1 && inner_iter == 1);
        if use_median
            % 初值：使用中位数法（公式26）
            sigma0 = compute_posterior_std_median(e_t, Q_e);
        else
            % 后续：使用传统加权残差形式
            sigma0 = compute_posterior_std_traditional(e_t, Q_c_current, R, m, n);
        end
        
        % 步骤d: 第二次标准化 - 计算标准化总残差（公式25）
        % 公式25: ~v_ti^(j+1) = |~e_ti^(j+1)| / (δ_0 * sqrt(Q_e_ti^(j+1)))
        % 这是第二次标准化：使用δ_0对总体残差进行标准化
        v_t = compute_standardized_residuals(e_t, Q_e, sigma0);
        
        % 步骤e: 基于标准化总残差计算IGGIII等价权函数得到R矩阵（公式27）
        % 公式27: 三段函数，根据标准化总残差v_ti的值确定权值R_i
        R = compute_IGGIII_weights(v_t, k0, k1);
        
        % 步骤f: 使用R矩阵进行抗差 - 计算等价协因数矩阵（公式28）
        % 公式28: Q̄_c = √(R Q_c R)
        % 这是抗差的核心步骤：通过R矩阵调整协因数矩阵，降低异常值的影响
        Q_c_bar = compute_equivalent_cofactor_matrix(Q_c_current, R);
        
        % 更新Q_c用于下一次迭代
        Q_c_current = Q_c_bar;
        
        % 步骤g: 更新矩阵并计算RTLS解
        [x_new, A_new, tls_iter] = update_matrices_and_solve(A, y, Q_c_current, x);
        tls_iter_total = tls_iter_total + tls_iter;  % 累积TLS迭代次数
        
        % 检查内部循环收敛（可选）
        if norm(x_new - x) < tol * 0.1
            break;
        end
        
        x = x_new;
        A_current = A_new;
    end
    
    iter_info.inner_iter(j) = inner_iter_count;
    
    % 步骤h: 检查外部循环收敛条件
    delta_x = x - x_prev;
    if norm(delta_x) <= tol
        iter_info.converged = true;
        break;
    end
    
    x_prev = x;
    A_current = A_new;
end

%% 最终结果计算
x_est = x;

% 计算观测残差（用于粗差检测）
V = y - A * x_est;

% 使用原始Q_c
Q_v = Q_c;
Q_v_diag = diag(Q_v);

% 使用稳健方法估计sigma0
% 策略：使用残差的低分位数估计，避免被异常cluster污染

% 计算标准化残差（用于稳健估计）
V_normalized = V ./ sqrt(Q_v_diag);
V_normalized(~isfinite(V_normalized)) = [];

% 方法1：使用50%分位数（中位数）
median_estimate = median(abs(V_normalized));
sigma0_median = 1.4826 * median_estimate;

% 方法2：使用60%分位数（排除最大40%的异常值）
V_abs_sorted = sort(abs(V_normalized));
n_data = length(V_abs_sorted);
percentile_60_idx = ceil(0.60 * n_data);
sigma0_p60 = V_abs_sorted(percentile_60_idx) / norminv(0.80);  % 60%分位数对应正态分布的倍数

% 方法3：使用最小二乘回归的残差标准差（较为保守）
std_estimate = std(V_normalized);
sigma0_std = std_estimate * 0.6;  % 使用60%的标准差

% 综合三种估计，取最小值（更激进的粗差检测）
sigma0_estimates = [sigma0_median, sigma0_p60, sigma0_std];
sigma0_final = min(sigma0_estimates);  % 改用最小值而不是中位数

% 确保sigma0在合理范围内
% 至少是中位数估计的80%
sigma0_min = sigma0_median * 0.8;
% 不超过标准差的70%
sigma0_max = std_estimate * 0.7;
sigma0_final = max(min(sigma0_final, sigma0_max), sigma0_min);

% 第二次标准化：计算标准化观测残差（用于粗差检测）
v_t_final = abs(V) ./ (sigma0_final * sqrt(Q_v_diag));
v_t_final(Q_v_diag <= 0) = 0;

% 基于标准化残差计算R矩阵（公式27）
R_final = compute_IGGIII_weights(v_t_final, k0, k1);

% 识别粗差点
rejected_idx = find(abs(v_t_final) > k1);  % 剔除点
downweighted_idx = find((abs(v_t_final) > k0) & (abs(v_t_final) <= k1));  % 降权点

% 计算落入排除区的观测方程数量
p = length(rejected_idx);

% 使用稳健估计的单位权中误差
sigma0_sq = sigma0_final^2;

% 计算参数估计的协因数矩阵
Q_x = compute_parameter_cofactor_matrix(A_current, Q_c_current);

% 另外计算总残差（用于算法理论分析）
[e_t_final, Q_e_final] = compute_total_residuals(A, y, x_est, Q_c);

% 计算内层迭代总数（使用实际累积的TLS迭代次数）
if exist('tls_iter_total', 'var')
    iter_info.inner_iter_total = tls_iter_total;  % 真正的TLS
else
    % 如果没有记录TLS迭代，使用重加权迭代次数作为备选
    if isfield(iter_info, 'inner_iter') && ~isempty(iter_info.inner_iter)
        iter_info.inner_iter_total = sum(iter_info.inner_iter);
    else
        iter_info.inner_iter_total = 0;
    end
end

% 保存粗差识别信息
iter_info.v_t_final = v_t_final;  % 标准化观测残差
iter_info.R_final = R_final;      % 权重向量
iter_info.rejected_idx = rejected_idx;  % 剔除点索引
iter_info.downweighted_idx = downweighted_idx;  % 降权点索引
iter_info.V = V;  % 观测残差
iter_info.e_t_final = e_t_final;  % 总体残差
iter_info.sigma0_final = sigma0_final;  % 最终sigma0
iter_info.sigma0_estimates = sigma0_estimates;  % 三种估计值
iter_info.sigma0_method = 'robust_median';  % 标记使用的方法

end

%% ========== 辅助函数 ==========

function val = getfield(s, field, default)
% 安全获取结构体字段值
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end

function [x_tls, A_corrected, iter_count] = TLS_solution(A, y, Q_c, x_init)
% TLS解算（基于EIV模型）
% 使用迭代方法求解TLS问题
%
% 注意：Q_c可以是m×m（只考虑y的误差）或m*(n+1)×m*(n+1)（考虑A和y的误差）
% 输出：
%   iter_count - TLS迭代次数（Newton迭代次数）

tol = 1e-8;
max_iter_tls = 50;
x = x_init;
[m, n] = size(A);
k = n + 1;  % 每个观测方程的误差个数（1个y的误差 + n个A的误差）
iter_count = 0;  % 初始化迭代计数

% 检查Q_c的维度
[Q_m, Q_n] = size(Q_c);

if Q_m == m && Q_n == m
    % Q_c是m×m，只考虑了y的误差
    % 需要构建完整的误差协因数矩阵，假设A的误差与y的误差独立
    % 对于EIV模型，构建块对角矩阵：每个观测方程有k个误差
    % 注意：这是完整的TLS实现，与论文一致
    Q_c_full = zeros(m*k, m*k);
    for i = 1:m
        idx_start = (i-1)*k + 1;
        % y的误差协因数来自Q_c
        Q_c_full(idx_start, idx_start) = Q_c(i, i);
        % A的误差协因数，假设与y的误差相同（可以根据实际情况调整）
        for j = 2:k
            Q_c_full(idx_start+j-1, idx_start+j-1) = Q_c(i, i);
        end
    end
    % 使用更稳定的求逆方法，并抑制警告
    cond_threshold = 1e-12;
    if rcond(Q_c_full) < cond_threshold
        warning('off', 'MATLAB:nearlySingularMatrix');
        Q_c_inv = pinv(Q_c_full);
        warning('on', 'MATLAB:nearlySingularMatrix');
    else
        Q_c_inv = inv(Q_c_full);
    end
    Q_c_full_used = true;
elseif Q_m == m*k && Q_n == m*k
    % Q_c已经是完整的误差协因数矩阵
    Q_c_full = Q_c;
    if rcond(Q_c) < 1e-12
        Q_c_inv = pinv(Q_c);
    else
        Q_c_inv = inv(Q_c);
    end
    Q_c_full_used = true;
else
    error('Q_c的维度不正确。应该是m×m或m*(n+1)×m*(n+1)');
end

% 使用完整的EIV模型方法（与论文一致）
if false  % 保持完整TLS实现，不使用简化方法
    % 简化方法：对于Q_c是m×m的情况，使用加权LS迭代
    % 这样可以避免构建巨大的EIV矩阵，提高数值稳定性
    P_c = Q_c_inv;
    
    for iter = 1:max_iter_tls
        iter_count = iter;  % 记录当前迭代次数
        % 计算残差
        v = y - A * x;
        
        % 使用加权LS更新
        % 对于线性回归，如果A的误差很小，可以直接使用LS
        A_Q_inv_A = A' * P_c * A;
        
        if rcond(A_Q_inv_A) < 1e-10
            A_Q_inv_A_inv = pinv(A_Q_inv_A);
        else
            A_Q_inv_A_inv = inv(A_Q_inv_A);
        end
        
        x_new = A_Q_inv_A_inv * (A' * P_c * y);
        
        % 检查收敛
        if norm(x_new - x) < tol
            x = x_new;
            A_corrected = A;  % 简化方法不修正A
            break;
        end
        
        x = x_new;
    end
    
    A_corrected = A;
else
    % 完整EIV模型方法
    for iter = 1:max_iter_tls
        iter_count = iter;  % 记录当前迭代次数
        % 计算残差
        v = y - A * x;
        
        % 构建B矩阵 (EIV模型)
        b = [1; x];  % 列向量
        B = kron(eye(m), b');  % m × m*k
        
        % 计算 P_v = (B * Q_c^(-1) * B')^(-1)
        Pv_inv = B * Q_c_inv * B';  % m × m
        
        % 使用更稳定的求逆方法，提高阈值以避免数值问题
        cond_threshold = 1e-12;
        if rcond(Pv_inv) < cond_threshold
            warning('off', 'MATLAB:nearlySingularMatrix');
            P_v = pinv(Pv_inv);
            warning('on', 'MATLAB:nearlySingularMatrix');
        else
            P_v = inv(Pv_inv);
        end
        
        % 计算误差估计
        vp = P_v * v;  % m × 1
        B_T_vp = B' * vp;  % m*k × 1
        e_hat = Q_c_inv * B_T_vp;  % m*k × 1
        
        % 重塑误差向量
        e_hat_reshaped = reshape(e_hat, k, m)';  % m × k
        
        % 提取系数矩阵的误差
        e_A = e_hat_reshaped(:, 2:end);  % m × n
        e_y = e_hat_reshaped(:, 1);  % m × 1
        
        % 修正系数矩阵和观测向量
        A_corr = A + e_A;
        y_corr = y + e_y;
        
        % 计算新的参数估计（使用修正后的A和y）
        % 提取y对应的权重
        P_c_y = zeros(m, m);
        for i = 1:m
            q_ii = Q_c_full((i-1)*k+1, (i-1)*k+1);
            if q_ii > 1e-12  % 避免除零
                P_c_y(i, i) = 1 / q_ii;
            else
                P_c_y(i, i) = 1e12;  % 给很小的协因数很大的权重
            end
        end
        
        A_Q_inv_A = A_corr' * P_c_y * A_corr;
        % 使用更稳定的求逆方法
        cond_threshold = 1e-12;
        if rcond(A_Q_inv_A) < cond_threshold
            warning('off', 'MATLAB:nearlySingularMatrix');
            A_Q_inv_A_inv = pinv(A_Q_inv_A);
            warning('on', 'MATLAB:nearlySingularMatrix');
        else
            A_Q_inv_A_inv = inv(A_Q_inv_A);
        end
        
        x_new = A_Q_inv_A_inv * (A_corr' * P_c_y * y_corr);
        
        % 检查收敛
        if norm(x_new - x) < tol
            x = x_new;
            A_corrected = A_corr;
            break;
        end
        
        x = x_new;
    end
end

x_tls = x;
if ~exist('A_corrected', 'var')
    A_corrected = A_corr;
end
end

function [e_t, Q_e] = compute_total_residuals(A, y, x, Q_c)
% 计算总残差（公式23）和协因数矩阵（公式24）
%
% 公式23: ~e_t^(j+1) = V - A^T δζ^(j+1)
% 其中V是观测残差，δζ是参数修正量
% 公式24: Q_e^(j+1) = Q_c - A^T (A Q_c^(-1) A^T)^(-1) A^T
%
% 注意：对于EIV模型，总残差是A和y的预测残差的组合

% 抑制矩阵奇异值警告（这些警告已经被pinv处理）
warning('off', 'MATLAB:nearlySingularMatrix');

[m, n] = size(A);
k = n + 1;

% 检查Q_c的维度
[Q_m, Q_n] = size(Q_c);

% 计算观测残差向量 V
V = y - A * x;

% 如果Q_c是m×m，使用简化的LS方法计算总残差
if Q_m == m && Q_n == m
    % 基于LS理论的线性化公式（公式23和24）
    % A是m×n，Q_c是m×m，A'是n×m
    % 抑制inv的警告
    warning('off', 'MATLAB:nearlySingularMatrix');
    Q_c_inv = inv(Q_c);
    warning('on', 'MATLAB:nearlySingularMatrix');
    
    % 计算 A' * Q_c^(-1) * A，结果是n×n
    A_Q_inv_A = A' * Q_c_inv * A;  % (n×m) * (m×m) * (m×n) = n×n
    
    % 使用更稳定的求逆方法，并抑制警告
    cond_threshold = 1e-12;
    if rcond(A_Q_inv_A) < cond_threshold
        warning('off', 'MATLAB:nearlySingularMatrix');
        A_Q_inv_A_inv = pinv(A_Q_inv_A);
        warning('on', 'MATLAB:nearlySingularMatrix');
    else
        warning('off', 'MATLAB:nearlySingularMatrix');
        A_Q_inv_A_inv = inv(A_Q_inv_A);
        warning('on', 'MATLAB:nearlySingularMatrix');
    end
    
    % 计算参数修正量（基于LS理论）
    % delta_zeta = (A' * Q_c^(-1) * A)^(-1) * A' * Q_c^(-1) * V
    % 简化：delta_zeta = A_Q_inv_A_inv * A' * Q_c_inv * V
    delta_zeta = A_Q_inv_A_inv * (A' * Q_c_inv * V);
    
    % 总残差（公式23）
    % e_t = V - A * delta_zeta
    % V是m×1，A是m×n，delta_zeta是n×1，所以A*delta_zeta是m×1 ✓
    e_t = V - A * delta_zeta;
    
    % 协因数矩阵（公式24）
    % 根据LS理论，残差的协因数矩阵是：
    % Q_e = Q_c - A * (A' * Q_c^(-1) * A)^(-1) * A'
    % Q_c是m×m，A是m×n，A_Q_inv_A_inv是n×n，A'是n×m
    % A * A_Q_inv_A_inv * A' 是 (m×n) * (n×n) * (n×m) = m×m ✓
    Q_e = Q_c - A * A_Q_inv_A_inv * A';
    
elseif Q_m == m*k && Q_n == m*k
    % Q_c是完整的误差协因数矩阵，使用EIV模型
    % 提取y对应的协因数矩阵（每个观测方程的第一个误差对应y）
    Q_c_y = zeros(m, m);
    for i = 1:m
        Q_c_y(i, i) = Q_c((i-1)*k+1, (i-1)*k+1);
    end
    
    % 使用简化的LS方法（类似于上面的情况）
    Q_c_y_inv = inv(Q_c_y);
    A_Q_inv_A = A' * Q_c_y_inv * A;  % n×n
    
    if rcond(A_Q_inv_A) < eps
        A_Q_inv_A_inv = pinv(A_Q_inv_A);
    else
        A_Q_inv_A_inv = inv(A_Q_inv_A);
    end
    
    delta_zeta = A_Q_inv_A_inv * (A' * Q_c_y_inv * V);
    e_t = V - A * delta_zeta;
    
    % 协因数矩阵（公式24）
    Q_e = Q_c_y - A * A_Q_inv_A_inv * A';
else
    error('Q_c的维度不正确。应该是m×m或m*(n+1)×m*(n+1)');
end

% 确保Q_e是对称正定的
Q_e = (Q_e + Q_e') / 2;

% 恢复警告
warning('on', 'MATLAB:nearlySingularMatrix');

% 如果Q_e不是正定的，使用修正方法
[V_eig, D_eig] = eig(Q_e);
D_eig = diag(max(diag(D_eig), 1e-10));  % 确保特征值非负
Q_e = V_eig * D_eig * V_eig';
end

function v_t = compute_standardized_residuals(e_t, Q_e, sigma0)
% 计算标准化总残差（公式25）- 第二次标准化
%
% 公式25: ~v_ti^(j+1) = |~e_ti^(j+1)| / (δ_0 * sqrt(Q_e_ti^(j+1)))
%
% 这是第二次标准化过程：
%   第一次标准化（在compute_posterior_std中）：
%     计算 ratio_i = |e_ti| / sqrt(Q_e_ti)，然后取中位数得到δ_0
%   第二次标准化（本函数）：
%     使用δ_0对总体残差进行标准化：v_ti = |e_ti| / (δ_0 * sqrt(Q_e_ti))
%
% 输入参数:
%   e_t    - 总体残差向量（公式23）
%   Q_e    - 总体残差的协因数矩阵（公式24）
%   sigma0 - 单位权后验标准差δ_0（公式26，第一次标准化的结果）

% 提取Q_e的对角线元素
Q_e_diag = diag(Q_e);

% 第二次标准化：使用δ_0对总体残差进行标准化
v_t = abs(e_t) ./ (sigma0 * sqrt(Q_e_diag));

% 避免除零
v_t(Q_e_diag <= 0) = 0;
end

function sigma0 = compute_posterior_std_median(e_t, Q_e)
% 计算单位权后验标准差（公式26）- 使用中位数法（仅用于初值）
%
% 公式26: δ_0 = 1.4286 * med_{i=1}^m (|~e_ti^(j+1)| / sqrt(Q_e_ti^(j+1)))
%
% 这是第一次标准化过程：
%   对每个观测方程计算 ratio_i = |e_ti| / sqrt(Q_e_ti)
%   然后取所有ratio_i的中位数，乘以1.4286得到δ_0
%   这个δ_0将用于第二次标准化（计算标准化总残差）
%
% 输入参数:
%   e_t - 总体残差向量（公式23）
%   Q_e - 总体残差的协因数矩阵（公式24）
%
% 输出参数:
%   sigma0 - 单位权后验标准差δ_0

% 提取Q_e的对角线元素
Q_e_diag = diag(Q_e);

% 第一次标准化：计算 |e_ti| / sqrt(Q_e_ti)
% 这是对总体残差的初步标准化（相对于各自的协因数）
ratio = abs(e_t) ./ sqrt(Q_e_diag);

% 避免除零和无效值
ratio(Q_e_diag <= 0 | ~isfinite(ratio)) = [];

if isempty(ratio)
    sigma0 = 1.0;
else
    % 计算中位数（使用中位数法提高抗差性）
    med_ratio = median(ratio);
    % 公式26：乘以1.4286得到单位权后验标准差
    sigma0 = 1.4286 * med_ratio;
    
    % 避免sigma0过小
    if sigma0 < 1e-10
        sigma0 = 1.0;
    end
end
end

function sigma0 = compute_posterior_std_traditional(e_t, Q_c, R, m, n)
% 使用传统加权残差平方和形式计算单位权中误差
% 公式: σ₀ = sqrt(v' * P * v / r)
%       其中 P = diag(R) * Q_c^(-1) * diag(R) 是等价权矩阵
%             r = m - n 是自由度
%
% 输入参数:
%   e_t - 总体残差向量
%   Q_c - 当前协因数矩阵
%   R   - 权重向量（来自IGGIII权函数）
%   m   - 观测数量
%   n   - 参数数量
%
% 输出参数:
%   sigma0 - 单位权中误差

% 构建等价权矩阵
% P = diag(R) * Q_c^(-1) * diag(R)
R_diag = diag(R);
Q_c_inv = inv(Q_c);
P = R_diag * Q_c_inv * R_diag;

% 计算加权残差平方和
rho = e_t' * P * e_t;

% 计算自由度
r = m - n;
if r <= 0
    r = max(1, m - n);
end

% 计算单位权中误差
sigma0 = sqrt(rho / r);

% 确保非负且合理
if ~isfinite(sigma0) || sigma0 < 1e-10
    sigma0 = 1.0;
end
end

function R = compute_IGGIII_weights(v_t, k0, k1)
% 计算IGGIII等价权函数（公式27）- 基于标准化总残差计算R矩阵
%
% 公式27: 分段函数，根据标准化总残差v_ti的值确定权值R_i
%   R_i = 1,                           if |~v_ti^(j+1)| ≤ k_0  (保守区)
%   R_i = (k_0 / |~v_ti^(j+1)|) * ((k_1 - k_0) / (k_1 - |~v_ti^(j+1)|))^2,  
%         if k_0 < |~v_ti^(j+1)| ≤ k_1  (关注区，降权)
%   R_i = 0,                           if |~v_ti^(j+1)| > k_1  (排除区，权重为0)
%
% 注意：这里与论文不同，将排除区的权重改为0（硬剔除）
%      而不是1e10，以真正排除粗差的影响
%
% 输入参数:
%   v_t - 标准化总残差向量（公式25，第二次标准化的结果）
%   k0  - 保守区阈值（推荐值：2.0-3.0）
%   k1  - 排除区阈值（推荐值：4.5-8.5）
%
% 输出参数:
%   R   - 权值向量，用于计算等价协因数矩阵进行抗差

m = length(v_t);
R = ones(m, 1);

abs_v_t = abs(v_t);

% 条件1: 保守区
idx_conserved = abs_v_t <= k0;
R(idx_conserved) = 1;

% 条件2: 关注区
idx_concern = (abs_v_t > k0) & (abs_v_t <= k1);
if any(idx_concern)
    v_concern = abs_v_t(idx_concern);
    R(idx_concern) = (k0 ./ v_concern) .* ((k1 - k0) ./ (k1 - v_concern)).^2;
end

% 条件3: 排除区 - 改为硬剔除（权重=0）
idx_exclude = abs_v_t > k1;
R(idx_exclude) = 0;  % 直接设为0，而不是1e10
end

function Q_c_bar = compute_equivalent_cofactor_matrix(Q_c, R)
% 计算等价协因数矩阵（公式28）- 抗差核心步骤
%
% 公式28: Q̄_c = √(R Q_c R)
% 其中R是对角矩阵，对角线元素为R_i（来自IGGIII权函数）
%
% 这是抗差的核心步骤：
%   通过R矩阵调整协因数矩阵，使得：
%   - 正常观测（保守区）：R_i=1，保持原权重
%   - 可疑观测（关注区）：0<R_i<1，降低权重
%   - 异常观测（排除区）：R_i=0，完全剔除
%
% 输入参数:
%   Q_c - 当前协因数矩阵
%   R   - 权值向量（来自IGGIII权函数）
%
% 输出参数:
%   Q_c_bar - 等价协因数矩阵，用于下一次TLS迭代

% 处理R=0的情况：将其设为极大值，使得权重≈0
% 对于R=0的点，Q_c_ii → ∞，权重 → 0
R_adjusted = R;
R_adjusted(R == 0) = 1e15;  % 将R=0替换为极大值

% 构建对角矩阵√R
R_diag = diag(sqrt(R_adjusted));

% 计算等价协因数矩阵
% Q̄_c = diag(√R) * Q_c * diag(√R)
Q_c_bar = R_diag * Q_c * R_diag;

% 确保对称性
Q_c_bar = (Q_c_bar + Q_c_bar') / 2;
end

function [x_new, A_new, tls_iter_count] = update_matrices_and_solve(A, y, Q_c, x_current)
% 更新矩阵并求解RTLS解

% 使用当前Q_c求解TLS
[x_new, A_new, tls_iter_count] = TLS_solution(A, y, Q_c, x_current);
end

function sigma0_sq = compute_posterior_variance_factor(e_t, Q_c, m, n, p)
% 计算单位权方差因子后验估计（公式29）
%
% 公式29: σ̂_0^2 = (ṽ_t^(j+1)^T Q_c^(-1) ṽ_t^(j+1)) / (m - 4 - p)

Q_c_inv = inv(Q_c);

% 计算分子
numerator = e_t' * Q_c_inv * e_t;

% 计算分母
denominator = m - 4 - p;

% 避免分母为零或负值
if denominator <= 0
    denominator = max(1, m - n);
end

sigma0_sq = numerator / denominator;

% 确保非负
if sigma0_sq < 0
    sigma0_sq = 0;
end
end

function Q_x = compute_parameter_cofactor_matrix(A, Q_c)
% 计算参数估计的协因数矩阵
% 基于TLS的协因数传播定律
%
% 根据LS理论，参数的协因数矩阵是：
% Q_x = (A' * Q_c^(-1) * A)^(-1)

[m, n] = size(A);

% 检查Q_c的维度
[Q_m, Q_n] = size(Q_c);

if Q_m == m && Q_n == m
    % Q_c是m×m
    Q_c_inv = inv(Q_c);
    % 计算 A' * Q_c^(-1) * A，结果是n×n
    A_Q_inv_A = A' * Q_c_inv * A;  % (n×m) * (m×m) * (m×n) = n×n
    
    if rcond(A_Q_inv_A) < eps
        Q_x = pinv(A_Q_inv_A);
    else
        Q_x = inv(A_Q_inv_A);
    end
else
    % Q_c可能是完整矩阵，提取y对应的部分
    k = n + 1;
    Q_c_y = zeros(m, m);
    for i = 1:m
        Q_c_y(i, i) = Q_c((i-1)*k+1, (i-1)*k+1);
    end
    Q_c_inv = inv(Q_c_y);
    A_Q_inv_A = A' * Q_c_inv * A;
    
    if rcond(A_Q_inv_A) < eps
        Q_x = pinv(A_Q_inv_A);
    else
        Q_x = inv(A_Q_inv_A);
    end
end

% 确保对称正定
Q_x = (Q_x + Q_x') / 2;
end

%% ========== 使用示例 ==========
%{
% 示例1: 基本使用
% 生成测试数据
m = 20;  % 观测数
n = 3;   % 参数个数
A = randn(m, n);
x_true = [1; 2; 3];
y = A * x_true + 0.1 * randn(m, 1);  % 添加噪声

% 添加一些异常值（粗差）
y([5, 10, 15]) = y([5, 10, 15]) + 5 * randn(3, 1);

% 设置协因数矩阵（单位矩阵）
Q_c = eye(m);

% 调用RTLS估计器
options.k0 = 2.5;
options.k1 = 6.0;
options.max_iter = 100;
options.tol = 1e-6;

[x_est, sigma0_sq, Q_x, iter_info] = RTLS_Eqn(A, y, Q_c, [], options);

fprintf('估计的参数:\n');
disp(x_est);
fprintf('真实参数:\n');
disp(x_true);
fprintf('单位权方差因子: %.6f\n', sigma0_sq);
fprintf('外部迭代次数: %d\n', iter_info.outer_iter);
fprintf('收敛状态: %s\n', char(string(iter_info.converged)));

% 示例2: 使用自定义初始值
x0 = [0.5; 1.5; 2.5];
[x_est2, ~, ~, iter_info2] = RTLS_Eqn(A, y, Q_c, x0, options);

% 示例3: 使用非单位协因数矩阵
% 假设观测值有不同的精度
weights = ones(m, 1);
weights([5, 10, 15]) = 0.1;  % 降低异常值的权重
Q_c_weighted = diag(1 ./ weights);

[x_est3, ~, ~, iter_info3] = RTLS_Eqn(A, y, Q_c_weighted, [], options);
%}

