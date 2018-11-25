digits(7)
syms f(x);
f(x) = input('f(x): ');
g = input('three guesses(in square brackets): ');
i = input('max iterations: ');
e = input('error: ');
z = int16.empty(i,0);
E = double.empty(i,0);
x1 = g(1);
x2 = g(2);
x3 = g(3);
iter = 0;
err = 100;
x = x1:0.01:x3;
y = f(x);
while((err>e)&&(iter<i))
    temp1 = vpa(f(x3) - f(x2));
    temp2 = vpa(x3 - x2);
    temp3 = vpa(f(x2) - f(x1));
    temp4 = vpa(x2 - x1);
    temp5 = vpa(temp1/temp2);
    temp6 = vpa(temp3/temp4);
    temp7 = vpa(x3 - x1);
    a = vpa((temp5 - temp6)/temp7);
    b = vpa(a*temp2 + temp5);
    c = f(x3);
    temp8 = vpa((b*b - 4*a*c)^(1/2));
    temp9 = abs(vpa(b - temp8));
    temp10 = abs(vpa(b + temp8));
    if(b>0)
        delx = abs(vpa(-2*c/temp10));
    end
    if(b<=0)
        delx = vpa(2*c/temp9);
    end
    err = abs(delx/x3);
    z(iter + 1) = iter + 1;
    E(iter + 1) = err;
    x1 = x2;
    x2 = x3; 
    x3 = x3 + delx;
    iter = iter + 1;
end
disp(x3);
plot(x,y);
plot(E,iter);
