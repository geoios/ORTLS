function [X, residuals, iter_info] = IRTLS_Mahboub2013(A, L, sigma_A, sigma_y, options)
% IRTLS_Mahboub2013 - Iteratively Reweighted Total Least Squares
% 
% 严格实现Mahboub等人2013年论文中的IRTLS算法
% 论文: "Iteratively reweighted total least squares: a robust estimation 
%        in errors-in-variables models"
%        Survey Review, 2013, Vol 45, No 329
%
% 算法严格按照论文公式(17)-(20)实现
%
% 输入:
%   A        - 系数矩阵 [n x m]，包含误差
%   L        - 观测向量 [n x 1]，包含误差
%   sigma_A  - A的标准差 [n x m]（每个元素的标准差）
%   sigma_y  - L的标准差 [n x 1]（每个观测的标准差）
%   options  - 选项:
%       .weight_type - 'huber'(公式21) 或 'proposed'(公式23)
%       .k          - Huber常数，默认1.5
%       .max_iter   - 最大迭代次数，默认100
%       .tol        - 收敛阈值，默认1e-6
%       .verbose    - 是否显示详细信息，默认true
%
% 输出:
%   X          - 参数估计 [m x 1]
%   residuals  - 残差结构
%   iter_info  - 迭代信息

%% 参数设置
if nargin < 5
    options = struct();
end
if ~isfield(options, 'weight_type'), options.weight_type = 'huber'; end
if ~isfield(options, 'k'), options.k = 1.5; end
if ~isfield(options, 'max_iter'), options.max_iter = 100; end
if ~isfield(options, 'tol'), options.tol = 1e-6; end
if ~isfield(options, 'verbose'), options.verbose = true; end
if ~isfield(options, 'damping'), options.damping = 0; end  % 阻尼因子(0-1)，0=无阻尼

[n, m] = size(A);
y = L;  % 使用论文中的符号

if options.verbose
    fprintf('========== IRTLS算法 (Mahboub 2013) - 严格实现 ==========\n');
    fprintf('观测数n: %d, 参数数m: %d\n', n, m);
    fprintf('权函数: %s, k=%.2f\n\n', options.weight_type, options.k);
end

%% Step 1: Schur分解Q_A和Cholesky分解Q_y (公式16)

% 构建协方差矩阵
Q_y = diag(sigma_y.^2);  % [n x n]
Q_A_diag = (sigma_A(:)).^2;  % [nm x 1]
Q_A = diag(Q_A_diag);  % [nm x nm]

% Cholesky分解 Q_y = G_y^T * G_y
try
    G_y = chol(Q_y);  % 上三角矩阵
catch
    % 如果不是正定，添加小的正则化项
    Q_y = Q_y + eye(n) * 1e-10;
    G_y = chol(Q_y);
end

% Schur分解 Q_A = S * Lambda * S^T (论文公式16)
% 对于对角矩阵，Schur分解就是其本身
[S, Lambda] = schur(Q_A);
% 对于对角协方差矩阵，S可以取为单位阵

if options.verbose
    fprintf('Step 1: Schur分解和Cholesky分解完成\n');
    fprintf('  Q_y范围: [%.2e, %.2e]\n', min(diag(Q_y)), max(diag(Q_y)));
    fprintf('  Q_A范围: [%.2e, %.2e]\n', min(Q_A_diag), max(Q_A_diag));
end

%% 初始化：使用TLS作为初始估计（对噪声数据更稳定）
% 选项1：加权最小二乘（原始）
% P_y_init = inv(Q_y);
% X_current = (A' * P_y_init * A) \ (A' * P_y_init * y);

% 选项2：TLS初始值（改进 - 对高噪声数据更稳定）
C = [A, -y];
[~, ~, V] = svd(C, 0);
v = V(:, end);
X_current = -v(1:m) / v(m+1);

% 如果TLS失败，回退到OLS
if any(isnan(X_current)) || any(isinf(X_current))
    P_y_init = inv(Q_y);
    X_current = (A' * P_y_init * A) \ (A' * P_y_init * y);
end

if options.verbose
    fprintf('\n初始估计x̂(0): ');
    fprintf('%.6f ', X_current);
    fprintf('\n\n');
end

%% Step 2-4: 迭代 (公式17-20)
iter_count = 0;
converged = false;
sigma0_history = zeros(options.max_iter, 1);
convergence_history = zeros(options.max_iter, 1);

% 初始权重矩阵（单位阵）
W_y = eye(n);
W_A = eye(n*m);

% 内层迭代参数
if ~isfield(options, 'max_inner_iter'), options.max_inner_iter = 50; end
if ~isfield(options, 'inner_tol'), options.inner_tol = 1e-8; end

% 累积内层迭代次数
inner_iter_total = 0;

while ~converged && iter_count < options.max_iter
    iter_count = iter_count + 1;
    X_prev_outer = X_current;  % 外层迭代前的参数
    
    %% 内层循环：在给定权重W_y和W_A下，迭代求解加权TLS问题（公式17-20）
    inner_iter_count = 0;
    X_inner = X_current;  % 内层迭代的初始值
    
    for inner_iter = 1:options.max_inner_iter
        inner_iter_count = inner_iter;
        X_inner_prev = X_inner;
        
        %% 公式(17): 计算 R_1^(i)
        % R_1^(i) = [G_y^T * W_y^(-1) * G_y + (I_n ⊗ x̂^(i-1))^T * S * W_A^(-1/2) * W_A^(-1/2) * S^T * (I_n ⊗ x̂^(i-1))]^(-1)
        % 注意：这里使用 (I_n ⊗ x̂) 而非论文符号中可能的 (x̂ ⊗ I_m)，因为EIV模型中vec(A)是[nm x 1]
        
        % 计算 Kronecker乘积: I_n ⊗ x̂ (正确的维度应该是[nm x n])
        I_n = eye(n);
        I_m = eye(m);
        In_kron_X = kron(I_n, X_inner);  % [nm x n]
        
        % 计算第一项
        W_y_inv = inv(W_y);
        term1 = G_y' * W_y_inv * G_y;
        
        % 计算第二项
        W_A_inv_sqrt = sqrt(inv(W_A));  % W_A^(-1/2)
        temp = S' * In_kron_X;  % S^T * (I_n ⊗ x̂), [nm x n]
        term2 = temp' * W_A_inv_sqrt * W_A_inv_sqrt * temp;  % [n x n]
        
        % R_1^(i) - 这是[n x n]矩阵
        try
            R_1 = inv(term1 + term2);
        catch
            R_1 = pinv(term1 + term2);
        end
        
        %% 公式(18): 计算 λ̂^(i)
        % λ̂^(i) = R_1^(i) * (y - A*x̂^(i-1))
        residual_y = y - A * X_inner;
        lambda_hat = R_1 * residual_y;
        
        %% 公式(19): 计算 R_2^(i)
        % R_2^(i) = (特殊矩阵) * S * W_A^(-1/2) * W_A^(-1/2) * S^T * (I_n ⊗ x̂^(i-1)) * R_1^(i)
        % 
        % 注意：这里需要一个 [m x nm] 的矩阵来提取每个参数的贡献
        % 标准的 kron(I_m, x^T) 维度是 [m x m]，不符合需求
        % 因此手动构造一个特殊的选择矩阵 [m x nm]：
        %   - 每行 i 对应参数 i
        %   - 在每个观测 j 的对应位置 (j-1)*m+i 处放置 x(i)
        % 这个矩阵从 vec(A) 中提取并加权每个参数的贡献
        Im_kron_XT = zeros(m, n*m);
        for i = 1:m
            for j = 1:n
                % 第i个参数在第j个观测的位置是 (j-1)*m + i
                Im_kron_XT(i, (j-1)*m + i) = X_inner(i);
            end
        end
        
        % 计算R_2
        temp_R2_1 = S' * In_kron_X * R_1;  % [nm x n] * [n x n] = [nm x n]
        temp_R2_2 = W_A_inv_sqrt * W_A_inv_sqrt * temp_R2_1;  % [nm x nm] * [nm x n] = [nm x n]
        temp_R2_3 = S * temp_R2_2;  % [nm x nm] * [nm x n] = [nm x n]
        R_2 = Im_kron_XT * temp_R2_3;  % [m x nm] * [nm x n] = [m x n]
        
        %% 公式(20): 更新参数 x̂^(i)
        % x̂^(i) = (A^T * R_1^(i) * A + R_2^(i) * A)^(-1) * (A^T * R_1^(i) * y + R_2^(i) * y)
        % 注意：R_2是[m x n]，A是[n x m]，所以R_2*A是[m x m]
        try
            term_ARA = A' * R_1 * A;  % [m x n] * [n x n] * [n x m] = [m x m]
            term_R2A = R_2 * A;  % [m x n] * [n x m] = [m x m]
            coeff_matrix = term_ARA + term_R2A;  % [m x m]
            rhs_vector = A' * R_1 * y + R_2 * y;    % [m x 1]
            X_inner_new = coeff_matrix \ rhs_vector;
        catch ME
            if options.verbose
                fprintf('\n错误发生在公式(20): %s\n', ME.message);
            end
            coeff_matrix = A' * R_1 * A + R_2 * A;
            rhs_vector = A' * R_1 * y + R_2 * y;
            X_inner_new = pinv(coeff_matrix) * rhs_vector;
        end
        
        % 更新内层迭代参数
        X_inner = X_inner_new;
        
        % 检查内层循环收敛
        inner_param_change = norm(X_inner - X_inner_prev) / (norm(X_inner_prev) + 1e-10);
        if inner_param_change < options.inner_tol
            break;
        end
    end
    
    % 累积内层迭代次数
    inner_iter_total = inner_iter_total + inner_iter_count;
    
    % 内层循环结束后，更新外层参数
    X_current = X_inner;
    
    % 应用阻尼更新（防止振荡）- 基于外层迭代
    if options.damping > 0 && iter_count > 1
        X_current = (1 - options.damping) * X_inner + options.damping * X_prev_outer;
    end
    
    %% 计算残差 (公式8和9) - 用于更新权重（使用最终的内层解）
    % ẽ_y = Q_y^T * λ̂ (公式8)
    e_y_tilde = Q_y' * lambda_hat;
    
    % ẽ_A = vec(Ẽ_A) = -Q_A * (I_n ⊗ x̂) * λ̂ (公式9)
    In_kron_X_final = kron(I_n, X_current);  % [nm x n]
    e_A_tilde = -Q_A * In_kron_X_final * lambda_hat;  % [nm x 1]
    
    %% 估计方差因子 σ_0^2 (公式22的简化版)
    % σ_0^2 = (ẽ_y^T * P_y * ẽ_y + ẽ_A^T * P_A * ẽ_A) / (n - m)
    P_y = inv(Q_y);
    P_A = inv(Q_A);
    
    sigma0_sq = (e_y_tilde' * P_y * e_y_tilde + e_A_tilde' * P_A * e_A_tilde) / (n - m);
    sigma0 = sqrt(max(sigma0_sq, 1e-10));
    sigma0_history(iter_count) = sigma0;
    
    %% Step 2: 更新权重矩阵 (公式21或23)
    
    % 标准化残差（添加数值稳定性保护）
    std_e_y = abs(e_y_tilde) ./ (sigma_y * sigma0 + 1e-10);
    std_e_A = abs(e_A_tilde) ./ (sigma_A(:) * sigma0 + 1e-10);
    
    % 根据权函数计算新权重
    switch lower(options.weight_type)
        case 'huber'
            % 公式(21): Huber权函数
            w_y = huber_weight(std_e_y, options.k);
            w_A = huber_weight(std_e_A, options.k);
            
        case 'proposed'
            % 公式(23): 改进权函数
            w_y = proposed_weight(std_e_y, options.k);
            w_A = proposed_weight(std_e_A, options.k);
            
        otherwise
            error('未知权函数: %s', options.weight_type);
    end
    
    W_y = diag(w_y);
    W_A = diag(w_A);
    
    %% 检查外层循环收敛
    param_change = norm(X_current - X_prev_outer) / (norm(X_prev_outer) + 1e-10);
    convergence_history(iter_count) = param_change;
    
    if options.verbose && (mod(iter_count, 5) == 0 || param_change < options.tol)
        fprintf('  迭代 %3d: σ₀=%.6f, Δx=%.2e, x̂=[', iter_count, sigma0, param_change);
        fprintf('%.6f ', X_current);
        fprintf(']\n');
    end
    
    % 改进的收敛判断：至少迭代3次，并且参数变化和sigma0变化都小于阈值
    if iter_count >= 3
        if param_change < options.tol
            % 额外检查：sigma0是否也稳定
            if iter_count > 1
                sigma0_change = abs(sigma0 - sigma0_history(iter_count-1)) / (sigma0_history(iter_count-1) + 1e-10);
                if sigma0_change < options.tol
                    converged = true;
                end
            else
                converged = true;
            end
        end
    end
end

%% 输出结果
X = X_current;

% 最终残差
residuals = struct();
residuals.v = y - A * X;  % 观测残差
residuals.e_y_tilde = e_y_tilde;  % 公式8
residuals.e_A_tilde = e_A_tilde;  % 公式9
residuals.lambda_hat = lambda_hat;  % Lagrange乘子

iter_info = struct();
iter_info.iterations = iter_count;  % 外层迭代次数
iter_info.inner_iter_total = inner_iter_total;  % 内层迭代总数
iter_info.converged = converged;
iter_info.sigma0 = sigma0_history(1:iter_count);
iter_info.convergence = convergence_history(1:iter_count);
iter_info.weight_type = options.weight_type;
iter_info.final_weights_y = diag(W_y);
iter_info.final_weights_A = diag(W_A);

if options.verbose
    fprintf('\n========================================\n');
    fprintf('✓ 算法%s (迭代%d次)\n', ...
        iif(converged, '收敛', '达到最大迭代次数'), iter_count);
    fprintf('最终σ₀: %.6f\n', sigma0);
    fprintf('最终参数估计x̂: ');
    fprintf('%.8f ', X);
    fprintf('\n');
    fprintf('观测残差v的范数: %.6f\n', norm(residuals.v));
    fprintf('========================================\n');
end

end

%% 辅助函数

function w = huber_weight(std_residual, k)
% Huber权函数 (公式21)
% [W]_jj = 1           if |e_j| ≤ σ_0*k
%        = σ_0*k/|e_j|  if |e_j| ≥ σ_0*k
w = ones(size(std_residual));
mask = std_residual > k;
w(mask) = k ./ std_residual(mask);
end

function w = proposed_weight(std_residual, k)
% 改进权函数 (公式23)
% [W]_jj = 1              if |e_j| ≤ σ_0*k
%        = (σ_0*k/|e_j|)²  if |e_j| ≥ σ_0*k
w = ones(size(std_residual));
mask = std_residual > k;
w(mask) = (k ./ std_residual(mask)).^2;
end

function result = iif(condition, true_val, false_val)
% 简单的三元运算符
if condition
    result = true_val;
else
    result = false_val;
end
end
