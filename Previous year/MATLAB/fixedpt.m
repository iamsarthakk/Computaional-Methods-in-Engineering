syms phi(x);
phi(x) = input('Enter phi of x: ');
n = double(input('initial value: '));
err = 100;
iter = 0;
e = input('minimum relative error: ');
z = double.empty(i,0);
E = double.empty(i,0);
while((err > e)&&(iter < i))
    if(imag(double(phi(n)))==0)
        y = phi(n);
    end
    if(imag(double(phi(n)))~=0)
        y = -1*abs(double(phi(n)));
    end
    if(y == n)
        err = 0;
        iter = iter + 1;
        break;
    end
    err = double(abs((y - n)/n));
    z(iter + 1) = iter + 1;
    E(iter + 1) = err;
    n = y;
    iter = iter + 1;
end


