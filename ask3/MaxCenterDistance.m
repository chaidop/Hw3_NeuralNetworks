function [max] = MaxCenterDistance(w, S)
    dist = 0;
    max = 0;
    for k = 1:S-1
        for i = k:S
            pow_t = (w(k) - w(k+i)).^2;
            disp(pow_t);
            dist = pow_t;
            dist = sqrt(dist);
            if(max<dist)
                max = dist;
            end
        end
    end
    