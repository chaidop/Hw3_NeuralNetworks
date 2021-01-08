function [n] = ADALINE_net_input(p, w, b, i, S)
%p: inputs se sthlodianysnma
%w: weights se sthlodianysnma
%b: biases se sthlodianysnma
%i: number of input data points
%n: net output, opou kathe sthlh i einai to net output gia to data point i
disp("ADALINEEEEEEEEEEEEEEEE");
disp(p);
disp(w);
disp(b);
n = [];
for i = 1:length(p)
    pp = [];
    for j = 1:S
        pp = [pp p(i)];% ftiaxnw tis diastaseis
    end
    
    disp(pp);
    n = [n w'*pp' + b];
    disp(n);
   % pause;
end
disp("THE N2 IS ");
disp(n);