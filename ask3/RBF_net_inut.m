function [n] = RBF_net_inut(p, w, b, i, S)
%p: inputs se sthlodianysnma
%w: weights se sthlodianysnma
%b: biases se sthlodianysnma
%i: number of input data points
%n: net output, opou kathe sthlh i einai to net output gia to data point i
disp("======= IN RBF NET INPUT");
%disp(p);
%disp(w);
%disp(b);
n = [];
    pp = [];
    for j = 1:S
        pp = [pp p];% ftiaxnw tis diastaseis
    end
    pp = pp';
    disp(pp);
    disp(b);
    pow_t = (pp - w).^2;
    disp(pow_t);
    n_t = sqrt(pow_t).*b;
    disp(n_t);
    n = [n n_t];
disp("THE N1 IS ");
disp(n);