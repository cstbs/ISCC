% =========================================================================
% 混合场通感一体化隐蔽通信系统仿真 (IEEE Transactions 级排版精装版)
% 最终定稿版：移除热力图，顺延图表编号 (共4个核心性能图)
% =========================================================================
clear; clc; close all;

%% 0. 全局图形与排版设置 (开启 LaTeX 渲染)
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);

%% 1. 系统核心物理参数设置
fc = 28e9;                  % 载波频率 28 GHz
c = 3e8;                    % 光速
lambda = c / fc;            % 波长 (~10.7 mm)
d_a = lambda / 2;           % 天线间距
N = 128;                    % 默认基站天线数
D = (N-1) * d_a;            % 阵列物理孔径
d_ray = 2 * D^2 / lambda;   % 瑞利距离 (N=128时约为86.4m)
P_max_dBm = 30;             % 最大发射功率 30 dBm
P_max = 10^((P_max_dBm - 30)/10); 
sigma2_dBm = -90;           % 噪声功率 -90 dBm 
sigma2 = 10^((sigma2_dBm - 30)/10); 
kappa = 0.1;                % 隐蔽门限参数

% 节点几何拓扑 (Willie 设定在 15m 深度近场以最大化展现混合场优势)
theta_W = 0;   d_W = 15;    % Willie 默认位置 (0°, 15m) -> 近场
theta_B = 0;   d_B = 200;   % Bob 默认位置 (0°, 200m)   -> 远场
theta_T = 30;  d_T = 200;   % Target 默认位置 (30°, 200m) -> 远场

% 信道大尺度衰落系数
beta_0 = (lambda / (4 * pi))^2; 
beta_B = beta_0 / (d_B^2);
n_idx = (-(N-1)/2 : (N-1)/2).'; % 对称阵列天线索引

%% 2. 基础信道矩阵与门限预计算 (供全局共享使用)
a_B_1 = get_a_far(theta_B, n_idx);
a_W_1 = get_a_near(theta_W, d_W, n_idx, lambda, d_a);
Phi_T_1 = get_Phi_T(theta_T, n_idx);

HB_bar = a_B_1 * a_B_1';
HW_bar = a_W_1 * a_W_1';
beta_W_1 = beta_0 / (d_W^2);
noise_margin_1 = kappa * sigma2 / beta_W_1; 
gamma_max_1 = P_max * max(eig(Phi_T_1));
gamma_th_1 = 0.1 * gamma_max_1; 

%% ========================================================================
%% 图2：可达隐蔽速率 vs. 窃听者距离 d_W (跨越瑞利距离)
%% ========================================================================
fprintf('正在求解 图2 (距离 d_W 的演进)...\n');
d_W_vec = linspace(10, 150, 15); 
Rate_proposed_d = zeros(size(d_W_vec));
Rate_farfield_d = zeros(size(d_W_vec));
Rate_upper_d    = zeros(size(d_W_vec));

a_W_far = get_a_far(theta_W, n_idx);
HW_far_bar = a_W_far * a_W_far';

for i = 1:length(d_W_vec)
    d_w_cur = d_W_vec(i);
    a_W_cur = get_a_near(theta_W, d_w_cur, n_idx, lambda, d_a);
    HW_cur_bar = a_W_cur * a_W_cur';
    noise_margin_cur = kappa * sigma2 / (beta_0 / (d_w_cur^2));
    
    [R_obj] = solve_proposed(HB_bar, HW_cur_bar, Phi_T_1, P_max, kappa, gamma_th_1, N, noise_margin_cur);
    Rate_proposed_d(i) = log2(1 + (beta_B * R_obj) / sigma2);
    
    [R_far] = solve_proposed(HB_bar, HW_far_bar, Phi_T_1, P_max, kappa, gamma_th_1, N, noise_margin_cur);
    Rate_farfield_d(i) = log2(1 + (beta_B * R_far) / sigma2);
    
    [R_up]  = solve_upper(HB_bar, Phi_T_1, P_max, gamma_th_1, N);
    Rate_upper_d(i) = log2(1 + (beta_B * R_up) / sigma2);
end

figure(2);
set(gcf, 'Position', [150, 150, 600, 450], 'Color', 'w');
plot(d_W_vec, Rate_upper_d, '--ks', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(d_W_vec, Rate_proposed_d, '-ro', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r'); 
plot(d_W_vec, Rate_farfield_d, '-b^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
plot([d_ray, d_ray], [0, max(Rate_upper_d)+1], '--k', 'LineWidth', 1.5);
text(d_ray+3, 1, sprintf('Rayleigh Distance \n$\\approx %.1f$ m', d_ray), 'Interpreter', 'latex', 'FontSize', 11);
grid on; ylim([0, max(Rate_upper_d)+1]);
xlabel('Distance of Willie $d_W$ (m)'); 
ylabel('Covert Rate (bps/Hz)');
legend('Performance Upper Bound', 'Proposed MF-ISCC Scheme', 'Baseline 1: Pure Far-Field Design', 'Location', 'NorthEast');
box on;

%% ========================================================================
%% 图3：可达隐蔽速率 vs. 天线数量 N
%% ========================================================================
fprintf('正在求解 图3 (天线数量 N 与 Scaling Law，耐心等待)...\n');
N_vec = 32:16:192; 
Rate_proposed_N  = zeros(size(N_vec));
Rate_farfield_N  = zeros(size(N_vec));
Rate_upper_N     = zeros(size(N_vec));

for i = 1:length(N_vec)
    N_cur = N_vec(i);
    n_idx_cur = (-(N_cur-1)/2 : (N_cur-1)/2).'; 
    
    a_B_N = get_a_far(theta_B, n_idx_cur);
    a_W_N = get_a_near(theta_W, d_W, n_idx_cur, lambda, d_a);
    a_W_far_N = get_a_far(theta_W, n_idx_cur);
    Phi_T_N = get_Phi_T(theta_T, n_idx_cur);
    
    HB_N_bar = a_B_N * a_B_N';  
    HW_N_bar = a_W_N * a_W_N';
    HW_far_N_bar = a_W_far_N * a_W_far_N';
    
    gamma_th_N = 0.1 * P_max * max(eig(Phi_T_N)); 
    
    [R_prop] = solve_proposed(HB_N_bar, HW_N_bar, Phi_T_N, P_max, kappa, gamma_th_N, N_cur, noise_margin_1);
    Rate_proposed_N(i) = log2(1 + (beta_B * R_prop) / sigma2);
    
    [R_far] = solve_proposed(HB_N_bar, HW_far_N_bar, Phi_T_N, P_max, kappa, gamma_th_N, N_cur, noise_margin_1);
    Rate_farfield_N(i) = log2(1 + (beta_B * R_far) / sigma2);
    
    [R_up] = solve_upper(HB_N_bar, Phi_T_N, P_max, gamma_th_N, N_cur);
    Rate_upper_N(i) = log2(1 + (beta_B * R_up) / sigma2);
end

figure(3);
set(gcf, 'Position', [300, 300, 600, 450], 'Color', 'w');
plot(N_vec, Rate_upper_N, '--ks', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(N_vec, Rate_proposed_N, '-ro', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(N_vec, Rate_farfield_N, '-b^', 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
grid on;
xlabel('Number of Transmit Antennas $N$'); 
ylabel('Covert Rate (bps/Hz)');
legend('Performance Upper Bound', 'Proposed MF-ISCC Scheme', 'Baseline 1: Pure Far-Field Design', 'Location', 'NorthWest');
box on;
disp('图表生成完毕！');

%% ========================================================================
%% 图4：可达隐蔽速率 vs. 总发射功率 P_max
%% ========================================================================
fprintf('正在求解 图4 (功率 P_max 的博弈)...\n');
P_vec_dBm = 20:2.5:40;
P_vec_W = 10.^((P_vec_dBm - 30)/10);
Rate_proposed_P  = zeros(size(P_vec_W));
Rate_baselineAN_P = zeros(size(P_vec_W));
Rate_upper_P     = zeros(size(P_vec_W));

for i = 1:length(P_vec_W)
    P_cur = P_vec_W(i);
    gamma_th_cur = 0.1 * P_cur * max(eig(Phi_T_1)); 
    
    [R_prop] = solve_proposed(HB_bar, HW_bar, Phi_T_1, P_cur, kappa, gamma_th_cur, N, noise_margin_1);
    Rate_proposed_P(i) = log2(1 + (beta_B * R_prop) / sigma2);
    
    [R_baseAN] = solve_baseline_AN(HB_bar, HW_bar, Phi_T_1, P_cur, kappa, gamma_th_cur, N, theta_T, theta_B, n_idx, noise_margin_1);
    Rate_baselineAN_P(i) = log2(1 + (beta_B * R_baseAN) / sigma2);
    
    [R_upper] = solve_upper(HB_bar, Phi_T_1, P_cur, gamma_th_cur, N);
    Rate_upper_P(i) = log2(1 + (beta_B * R_upper) / sigma2);
end

figure(4);
set(gcf, 'Position', [200, 200, 600, 450], 'Color', 'w');
plot(P_vec_dBm, Rate_upper_P, '--ks', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(P_vec_dBm, Rate_proposed_P, '-ro', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(P_vec_dBm, Rate_baselineAN_P, ':gd', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'g');
grid on;
xlabel('Maximum Transmit Power $P_{\max}$ (dBm)'); 
ylabel('Covert Rate (bps/Hz)');
legend('Performance Upper Bound', 'Proposed MF-ISCC Scheme', 'Baseline 2: AN-Aided ISAC Scheme', 'Location', 'NorthWest');
box on;

%% ========================================================================
%% 图5：可达隐蔽速率 vs. 感知精度约束
%% ========================================================================
fprintf('正在求解 图5 (感知精度 CRB 阈值)...\n');
gamma_ratio_vec = linspace(0.05, 0.5, 10);
Rate_proposed_gamma  = zeros(size(gamma_ratio_vec));
Rate_baselineAN_gamma = zeros(size(gamma_ratio_vec));
Rate_upper_gamma     = zeros(size(gamma_ratio_vec));

for i = 1:length(gamma_ratio_vec)
    gamma_th_cur = gamma_ratio_vec(i) * gamma_max_1;
    
    [R_prop] = solve_proposed(HB_bar, HW_bar, Phi_T_1, P_max, kappa, gamma_th_cur, N, noise_margin_1);
    Rate_proposed_gamma(i) = log2(1 + (beta_B * R_prop) / sigma2);
    
    [R_baseAN] = solve_baseline_AN(HB_bar, HW_bar, Phi_T_1, P_max, kappa, gamma_th_cur, N, theta_T, theta_B, n_idx, noise_margin_1);
    Rate_baselineAN_gamma(i) = log2(1 + (beta_B * R_baseAN) / sigma2);
    
    [R_up] = solve_upper(HB_bar, Phi_T_1, P_max, gamma_th_cur, N);
    Rate_upper_gamma(i) = log2(1 + (beta_B * R_up) / sigma2);
end

figure(5);
set(gcf, 'Position', [250, 250, 600, 450], 'Color', 'w');
plot(gamma_ratio_vec, Rate_upper_gamma, '--ks', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(gamma_ratio_vec, Rate_proposed_gamma, '-ro', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot(gamma_ratio_vec, Rate_baselineAN_gamma, ':gd', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'g');
grid on;
xlabel('Normalized Sensing CRB Threshold $\gamma_{\rm th} / \gamma_{\rm max}$'); 
ylabel('Covert Rate (bps/Hz)');
legend('Performance Upper Bound', 'Proposed MF-ISCC Scheme', 'Baseline 2: AN-Aided ISAC Scheme', 'Location', 'NorthEast');
box on;


%% ========================================================================
%% 核心求解函数区 
%% ========================================================================
function a = get_a_far(theta_deg, n_idx)
    theta_rad = deg2rad(theta_deg);
    a = exp(1j * pi * n_idx * sin(theta_rad)); 
end

function a = get_a_near(theta_deg, d, n_idx, lambda, d_a)
    theta_rad = deg2rad(theta_deg);
    r_n = sqrt(d^2 + (n_idx * d_a).^2 - 2 * d .* n_idx * d_a * sin(theta_rad));
    a = exp(-1j * (2*pi/lambda) * (r_n - d)); 
end

function Phi = get_Phi_T(theta_deg, n_idx)
    theta_rad = deg2rad(theta_deg);
    a_T = exp(1j * pi * n_idx * sin(theta_rad));
    da_dtheta = 1j * pi * n_idx * cos(theta_rad) .* a_T;
    H_dot = da_dtheta * a_T' + a_T * da_dtheta';
    Phi = H_dot' * H_dot;
end

function [R_obj, wc_opt, Ws_opt] = solve_proposed(HB_bar, HW_bar, Phi_T, P_max, kappa, gamma_th, N, noise_margin)
    cvx_begin sdp quiet
        variable Wc(N,N) hermitian
        variable Ws(N,N) hermitian
        maximize( real(trace(HB_bar * Wc)) )
        subject to
            real(trace(HB_bar * Ws)) <= 1e-8;  
            real(trace(HW_bar * Wc)) <= kappa * real(trace(HW_bar * Ws)) + noise_margin;
            real(trace(Phi_T * Ws)) >= gamma_th;
            real(trace(Wc) + trace(Ws)) <= P_max;
            Wc >= 0;  Ws >= 0;
    cvx_end
    
    if ~contains(cvx_status, 'Solved') && ~contains(cvx_status, 'Inaccurate')
        R_obj = 0; wc_opt = zeros(N,1); Ws_opt = zeros(N,N);
    else
        R_obj = real(trace(HB_bar * Wc));
        [V_c, D_c] = eig(Wc); [~, idx_c] = max(diag(D_c)); 
        wc_opt = V_c(:, idx_c) * sqrt(abs(D_c(idx_c, idx_c)));
        Ws_opt = Ws; 
    end
end

function [R_obj] = solve_baseline_AN(HB_bar, HW_bar, Phi_T, P_max, kappa, gamma_th, N, theta_T, theta_B, n_idx, noise_margin)
    a_T = get_a_far(theta_T, n_idx);
    ws_dir = a_T / norm(a_T); 
    Phi_T_gain = real(ws_dir' * Phi_T * ws_dir);
    P_s = gamma_th / Phi_T_gain; 
    
    if P_s > P_max
        R_obj = 0; 
        return; 
    end
    
    Ws_fixed = P_s * (ws_dir * ws_dir');
    leakage_s = real(trace(HW_bar * Ws_fixed)); 
    
    a_B = get_a_far(theta_B, n_idx);
    wc_dir = a_B / norm(a_B);
    Wc_base = wc_dir * wc_dir';
    gain_B = real(trace(HB_bar * Wc_base));         
    leakage_c_factor = real(trace(HW_bar * Wc_base)); 
    
    cvx_begin sdp quiet
        variable p_c nonnegative
        variable Wan(N,N) hermitian
        
        maximize( p_c * gain_B )
        subject to
            real(trace(HB_bar * Wan)) <= 1e-8;
            p_c * leakage_c_factor <= kappa * (real(trace(HW_bar * Wan)) + leakage_s) + noise_margin;
            p_c + real(trace(Wan)) <= P_max - P_s;
            Wan >= 0;
    cvx_end
    
    if ~contains(cvx_status, 'Solved') && ~contains(cvx_status, 'Inaccurate')
        R_obj = 0; 
    else
        R_obj = p_c * gain_B; 
    end
end

function [R_obj] = solve_upper(HB_bar, Phi_T, P_max, gamma_th, N)
    cvx_begin sdp quiet
        variable Wc(N,N) hermitian
        variable Ws(N,N) hermitian
        maximize( real(trace(HB_bar * Wc)) )
        subject to
            real(trace(HB_bar * Ws)) <= 1e-8;
            real(trace(Phi_T * Ws)) >= gamma_th;
            real(trace(Wc) + trace(Ws)) <= P_max;
            Wc >= 0; Ws >= 0;
    cvx_end
    
    if ~contains(cvx_status, 'Solved') && ~contains(cvx_status, 'Inaccurate')
        R_obj = 0; 
    else
        R_obj = real(trace(HB_bar * Wc)); 
    end
end