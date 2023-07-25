N1 = 96;
addpath('./psicd', './opt', 'gsc', 'HierBF', 'bmwss');
N2 = 120;
u = 0;
v = 0;
alpha = 1/4;
beta = 1/4;

[W1Golay, W2Golay] = gen_golay_hier(N1, N2, 1/alpha, 1/beta, u, v);

WBmwss = gen_bmwss_matrix(N1, N2, alpha, beta, u, v);
WEBmwss = gen_ebmwss_matrix(N1, N2, alpha, beta, u, v);

WGsc = gen_gsc_matrix(N1, N2, alpha, beta, u, v);

rfNum = 2; digitNum = ceil(log2(max(N1, N2)));
WPsicd = gen_ps_icd_matrix(N1, N2, alpha, beta, u, v, rfNum, digitNum);
power = abs(WPsicd(:).^2) / mean(abs(WPsicd(:).^2));
power = pow2db(power);
close all
figure
[f1,x1]=ecdf(power);
plot(x1,f1,'LineWidth',1.5)
hold on
plot([0, 0], [0, 1], '--', 'LineWidth',1.5)
legend("PSICD", "Others")
xlabel("Power / dB")
ylabel("CDF")
set(gca, "FontSize", 12)
set(gca,'LooseInset',get(gca,'TightInset'))


k = 1;
factor = 4; % FFT over sampling
init = 'gsc';
% k = alpha*beta; % passband weight
[W1Opt, W2Opt] = gen_opt_matrix(N1, N2, alpha, beta, u, v, k, factor, init);
W1Opt = quantize(W1Opt, max(N1, N2));
W2Opt = quantize(W2Opt, max(N1, N2));

W = {{WBmwss}, {WEBmwss}, {WPsicd}, {W1Golay, W2Golay}, {WGsc}, {W1Opt, W2Opt}};


montNum = 1e4;
frameLen = 1e4;
chan = 'awgn';
if isequal(chan, 'awgn')
    SNR = -12:2:4;
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
            errorNum= errorNum + tranceiver(bpsk, w, snr, alpha, beta, u, v, chan);
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
legend('BMWSS [13] ', 'EBMWSS [14]', 'PSICD [19]', 'HierBF', 'GSC', 'Numerical Opt', fontsize=12)
xlabel('SNR (dB)', FontSize=12)
ylabel('BER', FontSize=12)

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
u1 = (2*rand-1)*alpha/1 + u0;
end


