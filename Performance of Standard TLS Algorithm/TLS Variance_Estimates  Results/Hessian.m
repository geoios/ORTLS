function [H, e_A, B, e_hat] = Hessian(A, L, P, x)



    % 获取问题维度
    [n, m] = size(A);
    k = m + 1;
    
    % 只计算一次Q_e和I_n，不预提取块
    Q_e = inv(P);
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
    e_hat = -Q_e * B_T_vp;
        
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