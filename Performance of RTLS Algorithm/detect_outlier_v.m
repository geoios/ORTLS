function [detected, w_tests, v, x_hat, F_critical, results] = detect_outlier_v(A_obs, L_obs, P, alpha)
% DETECT_OUTLIER_V 基于残差v的粗差探测函数（分量压缩法）
%
% 功能：对EIV模型进行参数估计，并基于残差v进行粗差探测（w检验）
% 原理：使用分量压缩法估计参数，计算残差v的协因数阵，进行w检验
%
% 输入参数：
%   A_obs - 设计矩阵 [n x m]，n为观测点数量，m为参数个数
%   L_obs - 观测值向量 [n x 1]
%   P     - 权阵 [3 x n]，格式为 [py; px1; px2]'，其中：
%           py  - L的权（第1行）
%           px1 - A第1列的权（第2行）
%           px2 - A第2列的权（第3行，通常为常数项，权值很大如1e14）
%   alpha - 显著性水平（可选，默认0.1）
%
% 输出参数：
%   detected   - 粗差检测结果 [n x 1] 逻辑向量，true表示检测到粗差
%   w_tests    - w检验统计量 [n x 1]
%   v          - 残差向量 [n x 1]，v = L_obs - A_obs * x_hat
%   x_hat      - 参数估计值 [m x 1]
%   F_critical - F检验临界值（标量）
%   results    - 结构体，包含中间计算结果：
%                .Qv          - 残差的协因数阵
%                .sit0_1       - 单位权方差
%                .Sigma_e      - 残差的方差协方差阵
%                .e_hat        - 误差估计值
%                .iter         - 迭代次数
%
% 使用示例：
%   % 设置权阵
%   n = 20; m = 2;
%   py = 1/sigma_y^2 * ones(n, 1);
%   px1 = 1/sigma_x^2 * ones(n, 1);
%   px2 = 1e14 * ones(n, 1);
%   P = [py, px1, px2]';
%   
%   % 进行粗差探测
%   [detected, w_tests, v, x_hat, F_critical] = detect_outlier_v(A_obs, L_obs, P, 0.1);
%
% 注意：
%   1. 本函数已包含所有依赖的子函数，可独立运行
%   2. 权阵P的格式必须为3×n
%   3. 检测标准：|w| > F_critical 时判定为粗差

% 参数检查
if nargin < 4
    alpha = 0.05;  % 默认显著性水平
end

% 获取维度信息
[n, m] = size(A_obs);

% 检查输入维度
if size(L_obs, 1) ~= n
    error('L_obs的维度必须与A_obs的行数相同');
end

if size(P, 1) ~= 3 || size(P, 2) ~= n
    error('权阵P的维度必须为3×n');
end

% ========== 步骤1：参数估计（使用TLS_XG_newton3）==========
PP = diag(P(:));
[x_hat, e_hat, iter, Pv_inv_final] = TLS_XG_newton3(A_obs, L_obs, PP);

% ========== 步骤2：计算误差传播 ==========
Q_e = inv(PP);
[H, e_A, B, e] = Hessian(A_obs, L_obs, PP, x_hat);
[Sigma_e_L, de_dL, sit0_1] = simplified_error_v(A_obs, L_obs, P, x_hat, Q_e, H, e_A, B, e);
Sigma_e_a_total = extract_dv(A_obs, L_obs, x_hat, P, H, Q_e, sit0_1);
Sigma_e = Sigma_e_a_total + Sigma_e_L;

% ========== 步骤3：计算残差 ==========
v = L_obs - A_obs * x_hat;

% ========== 步骤4：计算残差的协因数阵（不包含σ₀²）==========
Qv = Sigma_e / sit0_1;  % 协因数阵

% ========== 步骤5：w检验 ==========
% 根据公式：w = v / (σ₀ * sqrt(Qv))
% 其中Qv是残差的协因数阵（不包含σ₀²），σ₀是单位权中误差
w_tests = v ./ (sqrt(sit0_1) * sqrt(diag(Qv)));

% ========== 步骤6：计算临界值（F检验）==========
df1 = 1;
df2 = n - m;
F_critical = sqrt(finv(1 - alpha, df1, df2));

% ========== 步骤7：粗差检测 ==========
detected = abs(w_tests) > F_critical;

% ========== 步骤8：组织输出结果 ==========
if nargout > 5
    results = struct();
    results.Qv = Qv;
    results.sit0_1 = sit0_1;
    results.Sigma_e = Sigma_e;
    results.e_hat = e_hat;
    results.iter = iter;
    results.H = H;
    results.e_A = e_A;
    results.B = B;
    results.e = e;
end

end

%% ========== 子函数定义 ==========

function [x_tls, e_hat, iter, Pv_inv_final] = TLS_XG_newton3(A_obs, L_obs, P)
% TLS_XG_newton3 使用牛顿法求解加权总体最小二乘问题
% 输入：
%   A_obs - 设计矩阵 [n x m]
%   L_obs - 观测值向量 [n x 1]
%   P     - 权矩阵 [k*n x k*n]，其中k=m+1
% 输出：
%   x_tls - 参数估计值 [m x 1]
%   e_hat - 误差估计值 [k*n x 1]
%   iter  - 迭代次数
%   Pv_inv_final - 最终的Pv_inv矩阵

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
    Pv_inv_final = [];

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
        
        % 保存最后一次的Pv_inv
        Pv_inv_final = Pv_inv;
        
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
        
        H1 = - A_corr' * P_v * A_corr2;

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
% Hessian 计算Hessian矩阵及相关中间变量
% 输入：
%   A - 设计矩阵 [n x m]
%   L - 观测值向量 [n x 1]
%   P - 权矩阵 [k*n x k*n]
%   x - 参数估计值 [m x 1]
% 输出：
%   H     - Hessian矩阵 [m x m]
%   e_A   - 系数矩阵的误差 [n x m]
%   B     - B矩阵 [k*n x n]
%   e_hat - 误差估计值 [k*n x 1]

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
        
    H1 = - A_corr' * P_v * A_corr2;

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

function [Sigma_e, de_dL, sit0_1] = simplified_error_v(A, L, P, X, Q_e, H, e_A, B, e)
% simplified_error_v 计算误差传播（L方向的误差）
% 输入：
%   A     - 设计矩阵 [n x m]
%   L     - 观测值向量 [n x 1]
%   P     - 权矩阵 [3 x n]
%   X     - 参数估计值 [m x 1]
%   Q_e   - 协因数矩阵
%   H     - Hessian矩阵
%   e_A   - 系数矩阵的误差
%   B     - B矩阵
%   e     - 误差向量（未使用，为兼容性保留）
% 输出：
%   Sigma_e - 误差向量e的方差协方差阵
%   de_dL   - 误差对观测值的偏导矩阵
%   sit0_1  - 单位权方差
        
    [n, m] = size(A);
    k = m + 1;
    
    % 提取权阵信息
    py = P(1, :)';
    Py = diag(py);    

    % 3. 计算 P_v = (B * Q_e * B')^{-1}
    Pv_inv = B * Q_e * B';
    if rcond(Pv_inv) < eps
        P_v = pinv(Pv_inv);
    else
        P_v = inv(Pv_inv);
    end
    
    % 计算单位权方差   
    In = eye(n);
    Gamma = (A + e_A)' * P_v;
    v = L - A * X;
    sit0_1 = (v' * P_v * v) / (n - m);

    if rcond(H) < eps
        dx_dL = -pinv(H) * Gamma;  % m×n
    else
        dx_dL = -H \ Gamma;        % 更高效，等价于 -inv(H) * Gamma
    end

    de_dL = In - A * dx_dL;
    
    % 步骤9：计算误差向量的方差协方差阵
    Sigma_L = sit0_1 * inv(Py);
    Sigma_e = de_dL * Sigma_L * de_dL';
end

function Sigma_e_a_total = extract_dv(A, L, X, P, H, Q_e, sit0_1)
% extract_dv 计算所有设计矩阵元素ai的误差传播
% 输入：
%   A      - 设计矩阵 [n x m]
%   L      - 观测向量 [n x 1]
%   X      - 参数估计 [m x 1]
%   P      - 权矩阵 [3 x n]，格式为[py, px1, px2]'
%   H      - Hessian矩阵
%   Q_e    - 协因数矩阵
%   sit0_1 - 单位权方差
% 输出：
%   Sigma_e_a_total - 总误差协方差矩阵（所有ai贡献之和）

    % 从P中提取设计矩阵元素的权阵
    [n, m] = size(A);
    Pa = cell(1, m);
    for i = 1:m
        Pa{i} = diag(P(i+1, :));
    end

    % 1. 计算基础变量
    v = L - A * X;
    k = m + 1;

    % 2. 计算B矩阵
    I_n = eye(n);
    B = kron(I_n, [1, X']);

    % 3. 计算P_v矩阵
    Pv_inv = B * Q_e * B';
    if rcond(Pv_inv) < eps
        P_v = pinv(Pv_inv);
    else
        P_v = inv(Pv_inv);
    end
    
    % 4. 计算e_hat和e_A
    e_hat = Q_e * B' * P_v * v;
    e_hat_reshaped = reshape(e_hat, k, n)';
    e_A = e_hat_reshaped(:, 2:end);

    % 5. 初始化输出变量
    Sigma_e_a_total = zeros(n, n);
    de_da_all = cell(1, m);
    dx_da_all = cell(1, m);
    
    % 6. 对每个设计矩阵参数ai进行循环计算
    for param_idx = 1:m  % 对每个参数列（A1, A2, ...）
        Sigma_a_i = sit0_1 * inv(Pa{param_idx});
        
        % 初始化dF_da_i矩阵
        dF_da_i = zeros(m, n);
        term1 = zeros(n, n);
        
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
            dF_da_single = (R1' + dET_da) * P_v * v + (A + e_A)' * P_v * R2;
            dF_da_i(:, obs_idx) = dF_da_single;
            
            % 计算term1
            da = zeros(n, m);
            da(obs_idx, param_idx) = 1;
            term1(:, obs_idx) = da * X;
        end
        
        if rcond(H) < eps
            dx_da_i = -pinv(H) * dF_da_i;  % m×n
        else
            dx_da_i = -H \ dF_da_i;        % 更高效，等价于 -inv(H) * de_da_i
        end
        dx_da_all{param_idx} = dx_da_i;

        de_da = -(term1 + A * dx_da_all{param_idx});
        
        % 保存该参数的所有偏导
        de_da_all{param_idx} = de_da;

        % 计算该参数贡献的协方差矩阵
        Sigma_e_i = de_da_all{param_idx} * Sigma_a_i * de_da_all{param_idx}';

        % 累加到总协方差
        Sigma_e_a_total = Sigma_e_a_total + Sigma_e_i;
    end
end
