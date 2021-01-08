%problem11 question B
clear all;
clf;

a = 0.01;%amax is around 0.2113
w1 = [-1 1];
b1 = [0.5 0.5];
b2 = 0;
% w2 will be found through LLS

%inputs, plus the ones-vector for the bias
p = [1 0 -1];
% ta targets. To t[i] einai to target gia thn eisodo p[i]
t = [-1 0 1];
no_inputs = length(p);
no_weights = length(w1);
r = 0;%r = 0 gia to A, B ,C erwthma kai r = 4 gia to D
I = eye(no_inputs);

n1 = RBF_net_inut(p', w1', b1', no_inputs);
%to n einai sthn i tou sthlh to net input gia to shmeio eisodou p[i]
%ara periexei to net input gia ola ta shmeia eisodou
a1 = radbas(n1, no_inputs);


syms w12 w22;
%dianysma me tis agnwstes parametrous mou(to bias to kserw)
x = [w12 w22 0];
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
c = (t*t');
disp("C = ");
disp(c);
%d = -2*U^T*t
d = -2*U'*t';
disp("d = ");
disp(d);
%A = 2[U^TU + rI]
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
xstar = inv(R)*(U'*t');
disp(inv(R));
disp("And the best weghts are: ");
disp(xstar);
ze = zeros(3);
plot([10*(-xstar(2)), 10*xstar(2)], [10*(xstar(1)),10*(-xstar(1))],"b");
grid on;
hold on;
