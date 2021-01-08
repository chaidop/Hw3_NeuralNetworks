function [a] = radbas(n, i)
%n: arithmos shmeiwn input, opou kathe sthlh einai gia ena input point
%a: h eksodos tou RBF layer opou kathe eshtlh einai 
%i: arithos input points
a = [];
for k = 1:i
    disp("---->>>> k ");
    disp(k);
    a_t = exp(-(n(:,k).^2));
    disp(a_t);
    a = [a a_t];
end
disp("a=");
 disp(a);