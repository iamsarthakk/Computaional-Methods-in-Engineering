digits(7);
syms f1(x);
f1(x) = input('f`(x): ');
g = (input('initial guess: '));
e = input('error: ');
iter = 0;
err = 100;
i = input('Enter the max iterations: ');
while ( (iter < i) && (err > e))
    Fi = vpa(f(g));
    fI = vpa(f1(g));
    temp = vpa(Fi/fI);
    x1 = vpa(g - temp);
    if(x1 == g)
        err = 0; 
        iter = iter + 1;
        break;
    end
    err = abs(vpa((x1 - g)/x1));
    z(iter + 1) = iter + 1;
    E(iter + 1) = err;
    iter = iter + 1;
    g = x1;
end
fprintf('Root is %0.6f in %d iterations',g,iter);