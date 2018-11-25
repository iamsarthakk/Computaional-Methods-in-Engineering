syms f(x);
prompt = 'Enter your function' ;
f(x) = input(prompt);
r = input('enter the range of root: ');
i = input('max number of iterations: ');
n = input(['Choose method.',... 
    '\n1. Bisection Method.',... 
    '\n2. False Position Method',... 
    '\n3. Fixed Point Method',... 
    '\n4. Newton-Raphson Method',... 
    '\n5. Secant Method: ']);
z = int16.empty(i,0);
E = double.empty(i,0);
    switch n
        case 1 
                prompt3 = 'enter max approx error' ;
                e = input(prompt3) ;
                err = 100;
                iter = 0;
                l = r(1);
                u = r(2);
                fl = f(l);
                fu = f(u);
                while((err > e)&&(iter < i))

                        m = (l + u)*0.5;
                        fm = f(m);
                        if( fl * fm < 0)
                                err = (u - m)/m;
                                u = m;
                                fu = f(m);
                        end
                        if( fu * fm < 0)
                                err = (m - l)/m;
                                l = m;
                                fl = f(m);
                        end
                        if( fm == 0)
                            err = 0;
                        end
                        E(iter + 1) = err;
                        z(iter + 1) = iter + 1;
                        iter = iter + 1;
                end     
                fprintf('Root is %0.6f in %d iterations',m,iter);
                
        case 2  
                prompt3 = 'enter max approx error' ;
                e = input(prompt3) ;
                err = 100;
                iter = 0;
                l = r(1);
                u = r(2);
                fl = f(l);
                fu = f(u);
                while((err > e)&&(iter < i))
                    m = u - (fu*(u - l)/(fu - fl));
                    fm = f(m);
                        if( fl * fm < 0)
                                err = (u - m)/m;
                                u = m;
                                fu = f(m);
                        end

                        if( fu * fm < 0)
                                err = (m - l)/m;
                                l = m;
                                fl = f(m);
                        end
                        if( fm == 0)
                            err = 0;
                        end
                        z(iter + 1) = iter + 1;
                        E(iter + 1) = err;
                        iter = iter + 1;
                end
                fprintf('Root is %0.6f in %d iterations',m,iter);
        case 3
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
                fprintf('Root is %0.6f in %d iterations',n,iter);
        case 4
                digits(7);
                syms f1(x);
                f1(x) = input('f`(x): ');
                g = (input('initial guess: '));
                e = input('error: ');
                iter = 0;
                err = 100;
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
        case 5
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
                    iter = iter + 1;
                end
                fprintf('Root is %0.6f in %d iterations',m2,iter);
        otherwise 
            disp('error');
    end
x = r(1):0.01:r(2);
y = f(x);
figure
subplot(2,1,1)
plot(x,y)
title('f(x) vs x')

subplot(2,1,2)
plot(z,E)
title('error vs iteration number')
               