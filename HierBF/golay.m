function seq = golay(len, type)
e2 = [1, 1; 1, -1];
e10 =[1 -1 -1 1 -1 1 -1 -1 -1 1;
    1 -1 -1 -1 -1 -1 -1 1 1 -1];
e26 = [-1,1,-1,-1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1;
    -1,1,-1,-1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1];
e3 = [1, 1, -1; 1, 1i, 1];
e5 = [1i, 1i, 1, -1, 1;
    1i, 1, 1, 1i, -1];
e11 = [1, 1i, -1, 1, -1, 1i, -1i, -1, 1i, 1i, 1;
    1, 1, -1i, -1i, -1i, 1, 1, 1i, -1, 1, -1];
e13 = [1, 1, 1, 1i, -1, 1, 1, -1i, 1, -1, 1, -1i, 1i;
    1, 1i, -1, -1, -1, 1i, -1, 1, 1, -1i, -1, 1, -1i];

expo = len_dec(len, type);
if length(expo) == 7
    [a, uc, ue, b, c, d, e] = expo{:};
    cSeq = cell_build({e3, e5, e11, e13}, [b, c-uc, d, e-ue]);
else
    [a, uc, ue] = expo{:};
    cSeq = cell(0);
end
rSeq = cell_build({e2, e10, e26}, [a, uc, ue]);
N1 = length(rSeq);
N2 = length(cSeq);

% compositions of quaternary sequences
if N2 > 0
    seq = cSeq{1};
    for i = 2:N2
        s1 = rSeq{i-1};
        s2 = cSeq{i};
        seq = golay_glue(s1(1, :), s1(2, :), seq(1, :), seq(2, :), s2(1, :), s2(2, :));
    end
else
    seq = [1; 1];
end

% compositions of binary sequences and quaternary sequences
if N2 == 0
    m = N1;
else
    m = N1-N2+1;
end
for i = 1:m
    s1 = rSeq{i+N1-m};
    seq = golay_glue(s1(1,:), s1(2, :), 1, 1, seq(1, :), seq(2, :));
end
end

function c = cell_build(e, num)
c = cell(1, sum(num));
for i = 1:length(e)
    for j = 1:num(i)
        c{j+sum(num(1:i-1))} = e{i};
    end
end
c = c(1, randperm(sum(num)));
end
    
    

    
