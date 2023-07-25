function seq = gen_golay_sequences(len)
try
    seq = golay(len, 'binary');
    return;
catch
    try
        seq = golay(len, 'quaternary');
        return;
    catch
    end
end

[expo, num] = expo_dec(len, 2);
assert(num<=245, 'Infeasible length!');
struct = load(sprintf('./atoms/cm_cc_%d.mat', num));
atom = struct.W;

for i = 1:expo
    atom = golay_double(atom);
end
seq = atom';
end

function seq = golay_double(atom)
seq(:, 1) = [atom(:, 1); atom(:, 2)];
seq(:, 2) = [atom(:, 1); -atom(:, 2)];
end

        
        