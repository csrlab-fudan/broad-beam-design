N1 = 96;
addpath('./psicd', './opt', 'gsc', 'HierBF', 'bmwss', 'ebmwss');
N2 = 120;
u = 0;
v = 0;
targetTheta = pi/2;
targetPhi = 0;
alpha = 1/4;
beta = 1/4;

[W1Golay, W2Golay] = gen_golay_hier(N1, N2, 1/alpha, 1/beta, u, v);

WBmwss = gen_bmwss_matrix(N1, N2, alpha, beta, u, v);
WEBmwss = gen_ebmwss_matrix(N1, N2, alpha, beta, u, v);

WGsc = gen_gsc_matrix(N1, N2, alpha, beta, u, v);

rfNum = 2; digitNum = ceil(log2(max(N1, N2)));
WPsicd = gen_ps_icd_matrix(N1, N2, alpha, beta, u, v, rfNum, digitNum);

k = 1; % stopband leakage ratio over initialization
factor = 4; % FFT over sampling
init = 'gsc';
% k = alpha*beta; % passband weight
[W1Opt, W2Opt] = gen_opt_matrix(N1, N2, alpha, beta, u, v, k, factor, init);
W1Opt = quantize(W1Opt, max(N1, N2));
W2Opt = quantize(W2Opt, max(N1, N2));

W = {{WBmwss}, {WEBmwss}, {WPsicd}, {W1Golay, W2Golay}, {WGsc}, {W1Opt, W2Opt}};

%% 
u = 0; v = 0;
angleLim = [pi/3, 2*pi/3, -pi/6, pi/6];
beamSampleHorizonNum = 180;
beamSampleVerticalNum = 180;
beamSampleLength = beamSampleHorizonNum*beamSampleVerticalNum;
broadbeampattern = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
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

close all
titles = {'BMWSS [13] ', 'EBMWSS [14]', 'PSICD [20]', 'HierBF', 'GSC', 'Numerical Opt', 'Mask'};
f1 = figure;
t1 = tiledlayout('flow');
f2 = figure;
t2 = tiledlayout('flow');
for ii = 1:length(W)
    w = W{ii};
    if length(w) == 2
        for i = 1:beamSampleHorizonNum
            temp = cos(beamThetaVec(i));
            w1 = w{1}; w2 = w{2};
            for j = 1:beamSampleVerticalNum
                u = sin(beamPhiVec(j));
                v = temp*cos(beamPhiVec(j));
                F  = exp(-1j*pi*(u*(0:N1-1).' + (v*(0:N2-1))));
                W1 = w1; W2 = w2;
                broadbeampattern(i, j) = 0.5*(abs(F(:).'*W1(:))^2 + abs(F(:).'*W2(:))^2);
            end
        end
    else
        for i = 1:beamSampleHorizonNum
            temp = cos(beamThetaVec(i));
            w1 = w{1};
            for j = 1:beamSampleVerticalNum
                u = sin(beamPhiVec(j));
                v = temp*cos(beamPhiVec(j));
                F  = exp(-1j*pi*(u*(0:N1-1).' + (v*(0:N2-1))));
                W1 = w1;
                broadbeampattern(i, j) = abs(F(:).'*W1(:))^2;
            end
        end
    end
%     figure
%     feature('DefaultCharacterSet','UTF-8');
%     [phi, theta] = meshgrid(beamPhiVec, beamThetaVec);
%     [x,y,z] = sph2cart(theta,phi,broadbeampattern);
%     p = mesh(x, y, z, broadbeampattern);
% %     view(180, 90)
%     axis equal
%     xlim([0, inf])
%     ylim([0, inf])
%     ax = gca;
%     set(gca,'XTick',[])
%     set(gca,'YTick',[])
%     set(gca,'ZTick',[])
%     set(gca, 'LooseInset', [0,0,0,0]);
%     set(gcf, "Position", [0, 0, 600, 600]);
%     hidden off;
%     light('position',[0 0 -5]);
%     hold on

    %  cross-section in azimuth plane

    figure(f1)
%     [x, y] = pol2cart(beamThetaVec, broadbeampattern(:, indexTargetPhi).');
%     plot(x, y)
%     axis equal
%     box off; axis off
%     set(gca, 'LooseInset', [0,0,0,0]);
%     set(gcf, "Position", [0, 0, 900, 300]);
%     view(90, 90)
    nexttile
    plot(180/pi*beamThetaVec, broadbeampattern(:, indexTargetPhi).', LineWidth=1.0);
    hold on
    plot(180/pi*beamThetaVec, maskTheta, '--', LineWidth=1.0);
    title(titles{ii})




    % cross-section in elevation plane
    figure(f2)
%     [x, y] = pol2cart(beamPhiVec, broadbeampattern(indexTargetTheta, :));
%     plot(x, y)
%     axis equal
%     box off; axis off
%     set(gca, 'LooseInset', [0,0,0,0]);
%     set(gcf, "Position", [0, 0, 900, 300]);
    nexttile
    plot(180/pi*beamPhiVec, broadbeampattern(indexTargetTheta, :), LineWidth=1.0);
    hold on
    plot(180/pi*beamPhiVec, maskPhi, '--', LineWidth=1.0);
    title(titles(ii))
end
figure(f1)
t1.TileSpacing = 'compact';
t1.Padding = 'compact';
xlabel(t1, '$\theta / \circ$', 'Interpreter','latex');
ylabel(t1, 'Power')
figure(f2)
t2.TileSpacing = 'compact';
t2.Padding = 'compact';
xlabel(t2, '$\varphi$', 'Interpreter','latex');
ylabel(t2, 'Power')

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
