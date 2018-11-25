digits(7);
l = r(1);
u = r(2);
e = input('error: ');
iter = 0;
err = 100;
while((err > e)&&(iter < i))
    fl = vpa(f(l));
    fu = vpa(f(u));
    temp1 = vpa(fu*vpa(u - l));
    temp2 = vpa(fu -fl);
    temp3 = vpa(fu*(u-l)/(fu-fl));
    m2 = vpa(u - temp3);
    err = vpa(abs((u - m2)/m2));
    z(iter + 1) = iter + 1;
    E(iter + 1) = err;
    l = u;
    u = m2;
end
fprintf('Root is %0.6f in %d iterations',m2,iter);
