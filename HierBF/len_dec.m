function expo = len_dec(num, type)
init = num;
switch type
    case 'binary'
        base = [26, 10, 2];
    case 'quaternary'
        base = [2, 3, 5, 11, 13];
    otherwise
        error('Binary or quaternary!');
end

for i = 1:length(base)
    [expo{i}, num] = expo_dec(num, base(i));
end

if num ~= 1
    error(sprintf('%d is undecomposable', init));
else
    if strcmp('binary', type)
        expo = flip(expo);
        return
    else
        [a, b, c, d, e] = expo{:};
        uSet = [];
        for u = 0:a
            if u <= c+e && b+c+d+e <= a+u+1
                uSet = [uSet, u];
            end
        end
        if ~isempty(uSet)
            % 随机选取一个u
            uRand = uSet(randperm(length(uSet)));
            u = uRand(1);
            % u中随机选一个符合条件的uc和ue, 表示10和26的个数
            uc = randi([max(u-e, 0), min(u, c)]);
            ue = u - uc;
            a = a-u;
            expo = {a, uc, ue, b, c, d, e};
        else
            error(sprintf('Exponentials are infeasible'));
        end
    end
end

function [expo, num] = expo_dec(num, base)
expo = 0;
while mod(num, base) == 0
    num = num / base;
    expo = expo + 1;
end

            
            