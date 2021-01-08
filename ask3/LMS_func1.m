function [e, w, b] = LMS_func1(p,t, a,a2, w, b,n, S)
%p: o pinakas pou periexei se kathe sthlh ta shmeia eisodou
%t: o pinakas-grammh pou periexei se kathe i-osth thesh tou to target tou
%shmeiou pou einai sthn i-osth sthlh tou p
%a: learning rate
%w: pinakas varwn
%b: bias
%n: arithmos shmeiwn input 
e = 0;
disp("LMSSSSSSSSSSSSSSSSSS");

    disp(w);
    disp(p);
    pp = [];
    for j = 1:S
        pp = [pp p];% ftiaxnw tis diastaseis
    end
    pp = pp';
    disp(pp);
    disp(b);
    disp(t);
    disp(a2);
    e=t- a2;
    disp("error = ");
    disp(e);
    
    w = w + (2*a*e)*p;
    %disp('w =');
    %disp(w);
    %disp('b =');
    %disp(b);
    b = b+ 2*a*e;
