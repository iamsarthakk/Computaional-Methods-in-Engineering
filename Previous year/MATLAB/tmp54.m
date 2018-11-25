digits(7)
a = input('Enter the coefficients(in square brackets) in increasing power of x: ');
err = input('enter error: ');
alpha = input('enter alpha0 and alpha1(in square brackets): ');
n = numel(a);
roots = double.empty(n,0);
cnt = 1;
p = 0;
c_p = double.empty(n,0);
while(p<n)
    c_p(p+1) = a(n-p);
    p = p + 1;
end
syms f(x)
f(x) = poly2sym(c_p);
while(n>3)
dela1 = 100;
dela2 = 100;
while((abs(dela1) > err)||(abs(dela2) > err))
d = double.empty(n,0);
d(n) = a(n);
d(n-1) = a(n-1) + alpha(2)*d(n);
j = n - 2;
while(j>0)
    d(j) = a(j) + alpha(2)*d(j+1) + alpha(1)*d(j+2);
    j = j - 1;
end
del = double.empty(n-1,0);
del(n-1) = d(n);
del(n-2) = d(n-1) + alpha(2)*del(n-1);
k = n-3;
while(k>0)
    del(k) = d(k+1) + alpha(2)*del(k+1) + alpha(1)*del(k+2);
    k = k - 1;
end
temp1 = vpa(d(2)*del(2) - d(1)*del(3));
temp2 = vpa(del(1)*del(3) - del(2)*del(2));
temp3 = vpa(del(2)*d(1) - del(1)*d(2));
dela2 = vpa(temp1/temp2);
dela1 = vpa(temp3/temp2);
alpha(1) = dela1 + alpha(1);
alpha(2) = dela2 + alpha(2);
end
temp4 = vpa(alpha(2)*alpha(2) + 4*alpha(1));
temp5 = vpa((temp4)^(1/2));
temp6 = vpa(0.5*(alpha(2) - temp5));
temp7 = vpa(0.5*(alpha(2) + temp5));
roots(cnt) = temp6;
cnt = cnt + 1;
roots(cnt) = temp7;
cnt = cnt + 1;
d(n) = a(n);
d(n-1) = a(n-1) + alpha(2)*d(n);
j = n - 2;
while(j>0)
    d(j) = a(j) + alpha(2)*d(j+1) + alpha(1)*d(j+2);
    j = j - 1;
end
b = n;
while(b>2)
    a(b-2) = d(b);
    b = b - 1;
end
n = n - 2;
end
if(n == 3)
    temp14 = vpa(a(2)*a(2) - 4*a(1)*a(3));
    temp15 = vpa((temp14)^(1/2));
    temp16 = vpa(0.5*(-1*a(2) + temp15)/a(3));
    temp17 = vpa(0.5*(-1*a(2) - temp15)/a(3));
    roots(cnt) = temp16;
    cnt = cnt + 1;
    roots(cnt) = temp17;
    cnt = cnt + 1;
end
if(n == 2)
    temp8 = -1*a(1)/a(2);
    roots(cnt) = temp8;
    cnt = cnt + 1;
end
disp(roots);


