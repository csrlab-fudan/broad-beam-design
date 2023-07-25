function [f, error] = ps_icd_optimize(v, rfNum, digitNum, iterNum)
step = 2*pi/(2^digitNum);
phaseSet = -pi:step:pi-step;
N = length(v);
f = zeros(size(v));

if rfNum==1
    for i = 1:N
        phase = quantize(angle(v(i)), phaseSet);
        f(i) = 1/sqrt(N) * exp(1j*phase);
    end
    error = [];
    return
end

Frf = exp(1j*2*pi*rand(N, rfNum));
i = 1;
error = zeros(iterNum, 1);
while i <= iterNum
    Fbb = (Frf'*Frf)\(Frf'*v);
    Frf = frf_optimize(Frf, Fbb, v, phaseSet);
    error(i) = norm(Frf*Fbb-v);
    i = i+1;
end

f = Frf*Fbb/norm(Frf*Fbb);
end

function phase = quantize(theta, phaseSet)
[~, idx] = min(abs(theta-phaseSet));
phase = phaseSet(idx);
end

function Frf = frf_optimize(Frf, Fbb, v, phaseSet)
[N, rfNum] = size(Frf);
for i = 1:N
    if rfNum==2
        Frf(i, :) = two_phase_opt(v(i), Fbb, phaseSet);
    else
        Frf(i, :) = more_phase_opt(v(i), Frf(i, :), Fbb, phaseSet);
    end
end
end

function frf = two_phase_opt(v, fbb, phaseSet)
alpha = abs(v);
beta = angle(v);
zeta1 = abs(fbb(1));
zeta2 = abs(fbb(2));
psi1 = angle(fbb(1));
psi2 = angle(fbb(2));
theta1 = beta-psi1+acos((alpha^2+(zeta1^2-zeta2^2))/(2*zeta1*alpha));
theta2 = beta-psi2-acos((alpha^2-(zeta1^2-zeta2^2))/(2*zeta2*alpha));
theta1 = quantize(theta1, phaseSet);
theta2 = quantize(theta2, phaseSet);
frf = exp(1j*[theta1, theta2]);
end

function frf = more_phase_opt(v, frf, fbb, phaseSet)
phase = angle(frf);
rfNum = length(frf);
t = 1;
phaseOld = phase;
while true
    p = mod(t-1, rfNum-2)+3;
    phaseOld(p) = phase(p);
    phase(p) = find_best_phase(v, frf, fbb, phaseSet, p);
    frf(p) = exp(1j*phase(p));
    t = t+1;
    if t > rfNum-2
        if isequal(phaseOld(3:end), phase(3:end))
            break
        end
    end
end
end

function phase = find_best_phase(v, frf, fbb, phaseSet, p)
N = length(phaseSet);
loss = zeros(N, 1);
for i = 1:N
    frf(p) = exp(1j*phaseSet(i));
    temp = v-frf(3:end)*fbb(3:end);
    frf(1:2) = two_phase_opt(temp, fbb(1:2), phaseSet);
    loss(i) = abs(v-frf*fbb);
end

[~, idx] = min(loss);
phase = phaseSet(idx);
end
        

        
