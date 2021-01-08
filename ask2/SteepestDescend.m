%% Uses Steepest Descend algorithm to calculate new weights and biases
%%for all the input data points
%%every row i has the new weights for the input data point p(i)
function  [wnew, bnew]= SteepestDescend(p, t, a, a1, w, b, n)
%p: o pinakas pou periexei se kathe sthlh ta shmeia eisodou
%t: o pinakas-grammh pou periexei se kathe i-osth thesh tou to target tou
%shmeiou pou einai sthn i-osth sthlh tou p
%a: learning rate
%a1: output of radbas
%w: pinakas varwn
%b: bias
%n: arithmos shmeiwn input 
    for i=1:n
    e=t(i)- a1(i);
    disp("t");
    disp(t(i));
    disp("a1");
    disp(a1);
    disp("e");
    disp(e);
    %%kanontas prakseis, to dF = de^2 = 2e*(de/dw) gia ta varh kai
    %%2e*(de/db)gia ta biases
    %%ara xnew = xold - a*2e*z, where x is the weight vector of length m, where the
    %%last element is the bias (bnew). Also, z is a vector of length 2,
    %%where its first element is the de/dw and the second is de/db
    %%We can also write: wnew = wold -as*2e*z(1)
    %%                   bnew = bold -as*2e*z(2)
    %%Apo prakseis, oi ananewseis tvn varwn einai:
    pp = [p(i) p(i)]';% ftiaxnw tis diastaseis
    pow_t = (pp - w).^2;
    n_t = sqrt(pow_t).*b;  
    disp((exp(-n_t)./n_t));
    disp(-exp(n_t).*sqrt(pow_t));
    disp((exp(-n_t)./n_t).*(pp - w).*b);
    z = [(exp(-n_t)./n_t).*(pp - w).*b; -exp(n_t).*sqrt(pow_t)];
    disp('wold =');
    disp(w);
    wnew = w;
    bnew = b;
    disp('2*a.*e) =');
    disp(2*a.*e);
    wnew = wnew + (2*a.*e).*z(1);
    disp('wnew =');
    disp(wnew);
    disp('bnew =');
    disp(bnew);
    bnew = bnew+ 2*a.*e.*z(2);
    w = wnew;
    b = bnew;
    end
