%% Needs backpropagation??

%15 neurons, points sample inputs
clear;
S = 8;%neurons of hidden layer (or centers). We experimented with 2, 4 and 8
S2 = 1;%neurons of output linear layer (ADALINE)
trainedFin = [];
predictedFin = [];
no_inputs = 10;
iter = 1;
grid on;
j = 1;
a = -0.5;
b = 0.5;
%w1 = linspace(-0.5, 0.5, 3)';
w1 = (b-a).*rand(S,1) + a;
%w1 = [0.5 0.5 0.2];
%b1 = linspace(-0.5, 0.5, 3)';
b1 = (b-a).*rand(S,1) + a;
%b1 = [0.2 0.5 -0.1];
%w2 = linspace(-0.5, 0.5, 3)';
w2 = (b-a).*rand(S,1) + a;
%w2 = [0.5 0.3 0.5];
%b2 = linspace(-0.5, 0.5, 3)';
b2 = (b-a).*rand(S2,1) + a;
maxiter = 100;

rng(0,'twister');
%Create a vector of 10 random values using the rand function to draw the 
%values from a uniform distribution in the open interval, (-2, 2).
a = -2;
b = 2;
p = (b-a).*rand(10,1) + a;
%get targets
f2=@(p) 1+sin(p.*(pi/8));
t=feval(f2,p);
a_rate = 0.06;

%{
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
%}

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

