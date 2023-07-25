function w = quantize(w, q)
[M, N] = size(w);
step = 2*pi/(q);
phaseSet = -pi:step:pi-step;
for i = 1:M
    for j= 1:N
        w(i, j) = exp(1j*round_phase(angle(w(i, j)), phaseSet));
    end
end
w = w/norm(w, 'fro');
end

function phase = round_phase(theta, phaseSet)
[~, idx] = min(abs(theta-phaseSet));
phase = phaseSet(idx);
end