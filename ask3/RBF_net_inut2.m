function [n] = RBF_net_inut2(p, w, b, i)
%p: inputs se sthlodianysnma
%w: weights se sthlodianysnma
%b: biases se sthlodianysnma
%i: number of input data points
%n: net output, opou kathe sthlh i einai to net output gia to data point i
disp(p);
disp(w);
disp(b);
n = [];
for k = 1:i
    disp(k);
    pp = [p(k) p(k)]';% ftiaxnw tis diastaseis
    disp(pp);
    pow_t = (pp - w).^2;
    disp(pow_t);
    n_t = sqrt(pow_t).*b;
    disp(n_t);
    n = [n n_t];
end  
disp("THE N1 IS ");
disp(n);