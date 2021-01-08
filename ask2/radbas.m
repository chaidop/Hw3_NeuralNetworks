function [a] = radbas(n, i)
%n: arithmos shmeiwn input, opou kathe sthlh einai gia ena input point
%a: h eksodos tou RBF layer opou kathe eshtlh einai 
%i: arithos input points
a = [];
    disp("---->>>> radbass ");
    %disp(k);
    a_t = exp(-(n(:,1).^2));
    disp(a_t);
    a = [a a_t];
%disp("a=");
disp(a);