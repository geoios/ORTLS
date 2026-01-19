 function [x_tls, P_v, iter] = TLS_XG_newton3(A_obs, L_obs, P)
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