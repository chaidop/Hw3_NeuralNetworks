function [dist] = SumErrors(p, w)
    sum = 0;
    for i = 1:length(p)
        pow_t = (p(i) - w).^2;
        disp(pow_t);
        sum = sum + pow_t;
    end
    dist = sqrt(sum);