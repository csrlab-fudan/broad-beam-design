%% beamforming matrices
N1 = 36;
addpath('gsc', 'golay', 'bmwss', 'ebmwss', 'opt', 'mllmi');
N2 = 64;
u0 = 0;
v0 = 0;
targetTheta = pi/2;
targetPhi = 0;
alpha = 1/4;
beta = 1/4;

[W1Golay, W2Golay] = gen_golay_hier(N1, N2, 1/alpha, 1/beta, u0, v0);

WBmwss = gen_bmwss_matrix(N1, N2, alpha, beta, u0, v0);
[WEBmwss, wEBMWSS1, wEBMWSS2] = gen_ebmwss_matrix(N1, N2, alpha, beta, u0, v0);
[WGsc, wGsc1, wGsc2] = gen_gsc_matrix(N1, N2, alpha, beta, u0, v0, "GSC");
[WGc, wGc1, wGc2] = gen_gsc_matrix(N1, N2, alpha, beta, u0, v0, "GC");
k = 1.2; % stopband leakage ratio over initialization
factor = 4; % FFT over sampling
init = 'gsc';
tic
[W1Opt, W2Opt] = gen_opt_matrix(N1, N2, alpha, beta, u0, v0, k, factor, init);
toc
time_opt = toc;
W1Opt = W1Opt/norm(W1Opt, "fro");
W2Opt = W2Opt/norm(W2Opt, "fro");
W1Opt = quantize(W1Opt, 64);
W2Opt = quantize(W2Opt, 64);

% tic
% WMLLMI = gen_mllmi_matrix(N1, N2, alpha, beta, u0, v0);
% toc
% time_mllmi = toc;
% disp(['Time of Opt: ', num2str(time_opt)]);
% disp(['Time of MLLMI: ', num2str(time_mllmi)]);
% WMLLMI = quantize(WMLLMI, 64);
load("wllmi_36_64.mat");
W_for_phase = {{wEBMWSS1(:), wEBMWSS2(:)}, {wGc1(:), wGc2(:)}, {wGsc1(:), wGsc2(:)},...
    {WMLLMI(:)}, {[W1Golay(:); W2Golay(:)]}, {[W1Opt(:); W2Opt(:)]}};
W = {{WEBmwss}, {WGc}, {WGsc}, {WMLLMI}, {W1Golay, W2Golay}, {W1Opt, W2Opt}};
% W = {{WEBmwss}, {WGc}, {WGsc}, {W1Golay, W2Golay}, {W1Opt, W2Opt}};
titles = {'EBMWSS ', 'GC', 'GSC', 'MLLMI', 'GASI', 'Opt'};
captions = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
%%  beampatterns and phases
u = 0; v = 0;
angleLim = [pi/3, 2*pi/3, -pi/6, pi/6];
beamSampleHorizonNum = 180;
beamSampleVerticalNum = 180;
beamSampleLength = beamSampleHorizonNum*beamSampleVerticalNum;
pattern = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
pattern_real = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
mode = 'uniform'; % the mode of choosing discrete angles
[beamThetaVec, beamPhiVec] = gen_angle_vec(beamSampleHorizonNum, beamSampleVerticalNum, mode, angleLim);

[~, indexTargetTheta] = min(abs(targetTheta-beamThetaVec));
[~, indexTargetPhi] = min(abs(targetPhi-beamPhiVec));
maskTheta = zeros(size(beamThetaVec));
maskPhi = zeros(size(beamPhiVec));
% stepTheta = round(1/2*abs(acos(u-alpha)-acos(u+alpha))/(2*pi)*length(beamThetaVec));
% stepPhi = round(1/2*abs(asin(v-beta)-asin(v+beta))/pi*length(beamPhiVec));
idxTheta = beamThetaVec <= acos(u-alpha) & beamThetaVec >=acos(u+alpha);
idxPhi = beamPhiVec <= asin(v+beta) & beamPhiVec >= asin(v-beta);
maskTheta(idxTheta) = 1/alpha/beta;
maskPhi(idxPhi) = 1/alpha/beta;
for i = 1:beamSampleHorizonNum
  maskTheta(i) = maskTheta(i)*gen_real_pattern(beamThetaVec(i), targetPhi);
end  
for i = 1:beamSampleVerticalNum
  maskPhi(i) = maskPhi(i)*gen_real_pattern(targetTheta, beamPhiVec(i));
end    

[~, right] = max(beamThetaVec >= acos(u-alpha));
[~, left] = max(beamThetaVec >= acos(u+alpha));
[~, top] = max(beamPhiVec >= asin(v-beta));
[~, down] = max(beamPhiVec >= asin(v+beta));
boundx = [left, left, right, right, left];
boundy = [down, top, top, down, down];

close all
f1 = figure;
t1 = tiledlayout('flow');
f2 = figure;
t2 = tiledlayout('flow');
f3 = figure;
t3 = tiledlayout('flow');
f4 = figure;
t4 = tiledlayout('flow');
f5 = figure;
t5 = tiledlayout('flow');
for ii = 1:length(W)
    w = W{ii};
    for i = 1:beamSampleHorizonNum
        temp = cos(beamThetaVec(i));
        for j = 1:beamSampleVerticalNum
            v = sin(beamPhiVec(j));
            u = temp*cos(beamPhiVec(j));
            F  = exp(-1j*pi*(u*(0:N1-1).' + (v*(0:N2-1))));
            if length(w)==2
                pattern(i, j) = 0.5*(abs(F(:).'*w{1}(:))^2 + abs(F(:).'*w{2}(:))^2);
            else
                pattern(i, j) = abs(F(:).'*w{1}(:))^2;
            end
            pattern_real(i, j) = pattern(i, j) * gen_real_pattern(beamThetaVec(i), beamPhiVec(j));
        end
    end

    figure(f1)
    nexttile
    x = 180/pi*beamThetaVec;
    y = pattern_real(:, indexTargetPhi).';
    plot(x, y, LineWidth=1.0);
    hold on
    plot(x, maskTheta, '--', LineWidth=1.0);
    title(captions{ii});


    figure(f2)
    nexttile
    plot(180/pi*beamPhiVec, pattern_real(indexTargetTheta, :).', LineWidth=1.0);
    hold on
    plot(180/pi*beamPhiVec, maskPhi, '--', LineWidth=1.0);
    title(captions{ii})

    draw_2d_pattern(pattern, f3, boundx, boundy, angleLim, captions{ii})
    draw_2d_pattern(pattern_real, f4, boundx, boundy, angleLim, captions{ii})

    figure(f5)
    nexttile
    plot_phase(W_for_phase{ii});
    title(captions{ii})
end
figure(f1)
xlabel(t1, '$\theta / \circ$', 'Interpreter','latex');
ylabel(t1, 'Normalized power / watt')
t1.TileSpacing = 'compact';
t1.Padding = 'compact';

figure(f2)
xlabel(t2, '$\varphi / \circ$', 'Interpreter','latex');
ylabel(t2, 'Normalized power / watt')
t2.TileSpacing = 'compact';
t2.Padding = 'compact';

figure(f3)
ylabel(t3, '$\varphi / \circ$', 'Interpreter','latex')
xlabel(t3, '$\theta / \circ$', 'Interpreter','latex')
t3.TileSpacing = 'compact';
t3.Padding = 'compact';

figure(f4)
ylabel(t4, '$\varphi / \circ$', 'Interpreter','latex')
xlabel(t4, '$\theta / \circ$', 'Interpreter','latex')
t4.TileSpacing = 'compact';
t4.Padding = 'compact';

%% ber
montNum = 1e6;
frameLen = 1e3;
chan = 'awgn';
if isequal(chan, 'awgn')
    SNR = -10:2:2;
else
    SNR = 0:4:30;
end
ber = zeros(length(W), length(SNR));
for i = 1:length(SNR)
    snr = db2pow(SNR(i));
    for j = 1:length(W)
        w = W{j};
        errorNum = 0;
        parfor k = 1:montNum
            bits = randi([0, 1], 1, frameLen);
            bpsk = 2*bits-1;
            errorNum= errorNum + tranceiver(bpsk, w, snr, alpha, beta, u0, v0, chan);
        end
        ber(j, i) = errorNum/(montNum*frameLen); 
    end
end

figure
marks = {'-^', '-v', '-o', '-*', '-s', '-d'};
for i = 1:size(ber, 1)
    semilogy(SNR, ber(i, :), marks{i}, 'LineWidth', 1.5);
    hold on
end
set(gca, 'LooseInset', [0,0,0,0]);
set(gca, 'FontSize', 12)
legend(titles, fontsize=12)
xlabel('SNR (dB)', FontSize=12)
ylabel('BER', FontSize=12)

%% functions
function plot_phase(A)
if length(A)==1
    a = A{1};
else
    a1 = A{1};
    a2 = A{2};
    if find_resolution(a1)>find_resolution(a2)
        a = a1;
    else
        a = a2;
    end
end

polarscatter(angle(a), ones(length(a), 1), 20);
rlim([0, 1])
ax = gca;
ax.RTickLabel = {};
ax.RGrid = 'off';
% ax.ThetaGrid = 'off';
set(gca, "FontSize", 9)
% set(gca,'LooseInset',get(gca,'TightInset'))
end

function resolution = find_resolution(a)
phases = mod(angle(a), 2*pi)/(2*pi);
phases = abs(phases - phases.');
resolution = 1/min(phases(phases>1e-6));
end

function draw_2d_pattern(pattern, f, boundx, boundy, angleLim, t)
figure(f)
nexttile
imagesc(pattern.')
xtickslabels = linspace(angleLim(1)*180/pi, angleLim(2)*180/pi, 4);
xticks = linspace(1, size(pattern, 2), numel(xtickslabels));
set(gca, 'XTick', xticks, 'XTickLabel', xtickslabels)
ytickslabels = linspace(angleLim(3)*180/pi, angleLim(4)*180/pi, 4);
yticks = linspace(1, size(pattern, 1), numel(ytickslabels));
set(gca, 'YTick', yticks, 'YTickLabel', ytickslabels)
hold on
plot(boundx, boundy, '--r', LineWidth=1.0)
colorbar
title(t)
end

function [theta_vec,phi_vec] = gen_angle_vec(horizon_sample_num, vertical_sample_num, mode, angleLim)
if nargin < 4
    angleLim = [0, 2*pi, -pi/2, pi/2];
end
theta0 = angleLim(1); theta1 = angleLim(2);
phi0 = angleLim(3); phi1 = angleLim(4);
assert(phi1==-phi0, 'Asymmetric elevation!');
if isequal(mode, 'random')
    theta_vec = theta0 + (theta1-theta0)*rand(1, horizon_sample_num);
    phi_vec =phi0 + (phi1-phi0)*rand(1, vertical_sample_num);
%     phi_vec = asin(sin(phi1)*(1-2*rand(1, vertical_sample_num)));
    theta_vec = sort(theta_vec);
    phi_vec = sort(phi_vec);
else
    dTheta = (theta1-theta0)/(horizon_sample_num-1);
    dPhi = (phi1-phi0)/(vertical_sample_num-1);
    theta_vec = theta0:dTheta:theta1;
    phi_vec = phi0:dPhi:phi1;
end
end

function A = gen_real_pattern(theta, phi)
theta = theta/pi*180;
phi = phi/pi*180;
Atheta = -min(12*((theta-90)/65)^2, 30);
Aphi = -min(12*(phi/65)^2, 30);
A = Atheta + Aphi;
A = db2pow(A);
end

function errorNum = tranceiver(bpsk, w, snr, alpha, beta, u0, v0, chan)
u1 = gen_rand_pos(u0, alpha);
v1 = gen_rand_pos(v0, beta);
[N1, N2] = size(w{1});
H = exp(-1j*pi*u1*(0:N1-1)).' * exp(-1j*pi*(v1*(0:N2-1)));
if isequal(chan, 'rayleigh')
    H = 1/sqrt(2)*(randn+1j*randn) * H;
end
N = length(bpsk);
bits = (bpsk>0);
noise = 1/sqrt(snr) * 1/sqrt(2)*(randn(1, N)+1j*randn(1, N));

if length(w) == 1
    y = abs(w{1}(:)'*H(:))* bpsk + noise;
else
    y = sqrt(0.5*(abs(w{1}(:)'*H(:))^2 + abs(w{2}(:)'*H(:))^2)) * bpsk + noise;
end
errorNum = sum(abs(bits-(real(y)>0)));
end

function u1 = gen_rand_pos(u0, alpha)
u1 = (2*rand-1)*alpha*0.9 + u0;
end
