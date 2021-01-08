%% Needs backpropagation??

%15 neurons, points sample inputs
clear;
no_inputs = 10;
p = linspace(-2, 2, no_inputs)';
S = 2;%neurons of hidden layer (or centers). We experimented with 2, 4 and 8
S2 = 1;%neurons of output linear layer (ADALINE)
trainedFin = [];
predictedFin = [];
iter = 1;
grid on;
j = 1;
a = -0.5;
b = 0.5;
w1 = linspace(-2, 2, S)';
%w1 = [0.5 0.5 0.2];
b1 = SumErrors(w1, S); %dialeksh 16, selida 8
%b1 = [0.2 0.5 -0.1];
%w2 = linspace(-2, 2, S)';
%w2 = [0.5 0.3 0.5];
%b2 = linspace(-2, 2, S2)';
maxiter = 100;

%get targets
f2=@(p) 1+sin(p.*(pi/8));
t=feval(f2,p);
a_rate = 0.06;

%%%%%%%%%%%%% FIND W2 AND B2 WITH LLS (ANALYTIC METHOD)%%%%%%%%%%%%%%%
p_t = p';
w1_t = w1';
b1_t = b1';
t_t = t';

r = 0;%r = 0 gia to A, B ,C erwthma kai r = 4 gia to D
I = eye(S +1);

n1 = RBF_net_inut2(p_t', w1_t', b1_t', no_inputs);
%to n einai sthn i tou sthlh to net input gia to shmeio eisodou p[i]
%ara periexei to net input gia ola ta shmeia eisodou
a1 = radbas2(n1, no_inputs);

%%% Analoga me to S, an exw S=4, tote exw sym w12 w22 w32 w42
syms w12 w22 b2;
%dianysma me tis agnwstes parametrous mou(to bias to kserw)
x = [w12 w22 b2];
%pragmatiko dianysma eisodou sto layer 2
%einai oles oi eisodoi pou ftanoun sto ADALINE layer, mazi me to bias 1
z = [];
for i = 1:no_inputs
    disp(i);
    z_t = [a1(:,i)' 1];
    disp(z_t);
    z = [z ; z_t];
    disp(z);
end
U = z;
disp("U' = ");
disp(U');

disp(x);
%a2 = x*z';

%%%%%%%%%%% LLS %%%%%%%%%%%%%%%%%%%%
%pairnw ola ta shmeia mazi
%c = t^Tt
c = (t_t*t_t');
disp("C = ");
disp(c);
%d = -2*U^T*t
d = -2*U'*t_t';
disp("d = ");
disp(d);
%A = 2[U^TU + rI]
disp(t);
disp(U'*U);
disp(U');
disp(U);
disp(U'*U + r.*I);
A = 2.*(U'*U + r.*I);
disp("A = ");
disp(A);
d = d';
x = x';
f  =  d*x;
ff = (1/2).*x'*A*x;
disp(f);
disp(ff);
%performance function
F = c + d*x + (1/2).*x'*A*x;
disp("Mean Square Error is:");
disp(F);


%%% Twra ksekina h evresh varwn (tou x vector) me th meleth ths F
%tha vroume tis idiotimes tou 2R(Hessianou)
flag = 0;
eigenvalues = eig(A);
for i = 1:length(eigenvalues)
    v = eigenvalues(i);
    disp("Eigenvalues");
    disp(v);
    if v<=0 
        flag = 1;%it has a weak minimm or no minimum on vector d
    end
end
if flag == 0;
    disp("There is a STRONG minimum!");
else
    disp("There is a WEAK minimum :(");
end

%Locating the stationary point
R= U'*U + r.*I;
disp(R);
xstar = inv(R)*(U'*t_t');

w2 = xstar(1:length(xstar)-1);
b2 = xstar(length(xstar));

disp("w2 and b2 are: ");
disp(w2);
disp(b2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thelei allagh analoga to S
w12 = w2(1); 
w22 = w2(2); 
b2 = b2;
x = [w12 w22 b2]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = c + d*x + (1/2).*x'*A*x;
disp("Sum Square Error after LLS for w2, b2 is:");
disp(F);
pause;

%%%%%%%%%%%%%%%%%%%%% FINISHED WITH LLS %%%%%%%%%%%%%%%%%%%%%%%%%%



while((j <= maxiter))
    for k =1:no_inputs
        %% training 1st RBF layer
        n1 = RBF_net_inut(p(k), w1, b1, no_inputs, S);
        %to n einai sthn i tou sthlh to net input gia to shmeio eisodou p[i]
        %ara periexei to net input gia ola ta shmeia eisodou
        a1 = radbas(n1, no_inputs);
        %[w1new , b1new]= SteepestDescend(p, t, a_rate, a1, w1, b1, no_inputs);

        %trainedFin = [trainedFin a1];

        %% training 2nd ADALINE layer
        n2 = ADALINE_net_input(p(k), w2, b2, no_inputs, S);
        a2 = purelin(n2);
        [e_all, w2 , b2] = LMS_func1(p(k), t(k), a_rate,a2, w2, b2, no_inputs, S);
        disp("w2 ========== ");
        disp(w2);
        disp("b2 =========== ");
        disp(b2);
        %%Backpropagation with steepest descend


        j = j +1;
    end
end


%% test output for training data
    %% training 1st RBF layer
    n1 = RBF_net_inut(p, w1, b1, no_inputs, S);
    %to n einai sthn i tou sthlh to net input gia to shmeio eisodou p[i]
    %ara periexei to net input gia ola ta shmeia eisodou
    a1 = radbas(n1, no_inputs);
  
    %trainedFin = [trainedFin a1];
    
    %% training 2nd ADALINE layer
    n2 = ADALINE_net_input(p, w2, b2, no_inputs, S);
    a2 = purelin(n2);
    disp("a2 = ");
    disp(a2);
    e = t - a2';
    disp("e all = ");
    disp(e);
    %% SUM SQUARED ERROR
    F = e.^2;
    tt = 0;
    for i = 1:length(e)
        tt = F(i) + tt;
    end
    disp("SSE = ");
    disp(tt); 
    for i = 1:no_inputs
        disp("======================");
        spf = sprintf('The REAL is %.2f',t(i));
        sspf = sprintf(' and TRAINED is %.2f.', a2(i));
        disp(spf);
        disp(sspf);
    end
    %pause;
hold on;
stem(p,t,'--or');
stem(p,a2,'-.*c');
legend('real', ' trained');
grid on;


