function find_roots
c=input('Is equation a polynomial? Y/N ','s');
if c=='N'
    fprintf('Choose one of the method:\n');
    y=input('1.Bisection\n2.False-position\n3.Fixed-Point\n4.Newton-Raphson\n5.Secant\n');
    
    switch y
        case 1,syms x
               z = input(' Enter The Equation ( Like sin(x) ): ','s');
               f = str2func(['@(x) ', z]);
               a = input('Enter two starting points (in brackets): ');
               fprintf('Enter the stopping criteria: \n');
               e1 = input('Enter the max. relative approximate error in successive iterations: ');
               e2 = input('How much closer to zero can function go: ');
               itr = input('Enter maximum number of iterations: ');
               if(f(a(1))*f(a(2))<0)
                   for i=1:itr
                       m=(a(1)+a(2))/2;
                      % E(itr) = abs(f(m));
                      
                       if abs(f(m))<e2
                           a1=m;
                           fprintf("%f",a1);
                           fprintf(' Termination due to convergence of function value');
                           break;
                       end
                       
                       if f(a(1))*f(m)<0
                           a(2)=m;
                       else
                           a(1)=m;
                       end
                       E(i) = a(2)-a(1);
                       iter(i) = i;
                        
                       if abs(a(2)-a(1))<e1
                           a1=m;
                           fprintf('%f',a1);
                           fprintf(' Termination due to convergence of interval');
                           break;   
                       end
                       z=i+1;
                       u=m;
                   end
                   if z==itr+1
                       fprintf('%f',u);
                       fprintf(' Termination due to limited number of iterations');
                   end
               fplot(f,[-50 50]);grid on;title('The Function');ylabel('f(x)');xlabel('x');
               figure;
               plot(iter,E);title('Error vs No. of iterations'); grid on; ylabel('Relative Error');xlabel('Number of iterations');
               else
                   fprintf('Wrong choice of initial points');
               end
               
        case 2,syms x
               z = input(' Enter The Equation ( Like sin(x) ): ','s');
               f = str2func(['@(x) ', z]);
               a = input('Enter two starting points (in bracket): ');
               fprintf('Enter the stopping criteria: \n');
               e1 = input('Enter the max. relative approximate error in successive iterations: ');
               e2 = input('How much closer to zero can function go: ');
               itr = input('Enter maximum number of iterations '); 
               if(f(a(1))*f(a(2))<0)
                   for i=1:itr
                       m=a(1)-(a(2)-a(1))/(f(a(2))-f(a(1)))*f(a(1));
                       if abs(f(m))<e2
                           a1=m;
                           fprintf("%f",a1);
                           fprintf(' Termination due to convergence of function value');
                           break;
                       end
                       
                       E(i) = a(2)-a(1);
                       iter(i) = i;
                       if f(a(1))*f(m)<0
                           a(2)=m;
                       else
                           a(1)=m;
                       end
                       
                       if abs(a(2)-a(1))<e1
                           a1=m;
                           fprintf('%f',a1);
                           fprintf(' Termination due to convergence of interval');
                           break;   
                       end
                       z=i+1;
                       u=m;
                   end
                   if z==itr+1
                       fprintf('%f',u);
                       fprintf(' Termination due to limited number of iterations');
                   end
               fplot(f);grid on;title('The Function');ylabel('f(x)');xlabel('x');
               figure;plot(iter,E);grid on;title('Error vs No. of iterations');ylabel('Relative Error');xlabel('Number of iterations');
               else
                   fprintf('Wrong choice of initial points');
               end
               
        case 3,syms x
               z = input('Enter function g(x) such that x=g(x): ','s');
               g = str2func(['@(x) ', z]);
               x1 = input('Enter initial guess: ');
               fprintf('Enter the stopping criteria: \n');
               e1 = input('Enter the max. relative approximate errors in successive iterations: ');
               e2 = input('How close much closer to zero can x-g(x) go: '); 
               itr = input('Enter maximum number of iterations '); 
               for i=1:itr
                   if abs(x1-g(x1))<e2
                       fprintf('%f',x1);
                       fprintf(' Termination due to convergence of function value');
                       break;
                   end
                   x2=g(x1);
                   E(i) = (x2-x1)/x1;
                   iter(i) = i;
                   
                   if abs((x2-x1)/x1)<e1
                       fprintf('%f',x2);
                       fprintf(' Termination due to convergence in relative approximate error');
                       break;
                   else
                       x1=x2;
                   end
                   if i==itr
                       fprintf('%f',x2);
                       fprintf('Termination due to limited number of iterations');
                       break;
                   end
               end
               
               fplot(g);grid on;title('Plot of y=g(x) and y=x');
               hold on
               fplot(x);legend('g(x)','x');
               hold off
               figure;
               plot(iter,E);grid on;title('Error vs No. of iterations');ylabel('Error');xlabel('Number of iterations');
               
        case 4,syms x
               z = input('Enter function f(x): ','s');
               f = str2func(['@(x) ', z]);
               z1 = input('Enter function derivative of f(x): ','s');
               S = vectorize(char(z)); % For Signs Like %,^,&,.
               f1 = str2func(['@(x) ', z1]);
               x1 = input('Enter starting point: ');
               fprintf('Enter the stopping criteria: \n');
               itr = input('Enter maximum number of iterations ');
               e1 = input('Enter the max. relative approximate errors in successive iterations: ');
               e2 = input('How close can funtion go to zero: ');
               for i=1:itr
                   
                   x2=x1-f(x1)/f1(x1);
                   E(i) = (x2-x1)/x1;
                   iter(i) = i;
                   
                   if abs(f(x1))<e2
                       fprintf('%f',x1);
                       fprintf(' Termination due to convergence of function value');
                       break;
                   end
                   if abs((x2-x1)/x1)<e1
                       fprintf('%f',x2);
                       fprintf(' Termination due to convergence of relative approximate error')
                       break;
                   else
                       x1=x2;
                   end
                   if i==itr
                       fprintf('%f',x2)
                       fprintf(' Termination due to limited number of iterations');
                       break;
                   end
               end
               fplot(f);grid on;title('The Function');ylabel('f(x)');xlabel('x');
               figure;
               plot(iter,E);grid on;title('Error vs No. of iterations');ylabel('Error');xlabel('Number of iterations');
               
        case 5,syms x
               z = input(' Enter The Equation ( Like sin(x) ): ','s');
               f = str2func(['@(x) ', z]);
               a = input('Enter two starting points (in bracket): ');
               fprintf('Enter the stopping criteria: \n');
               itr = input('Enter maximum number of iterations: ');
               e1 = input('Enter the max. relative approximate errors in successive iterations: ');
               e2 = input('How much close can function go to the zero: ');
               for i=1:itr
                   m=a(2)-f(a(2))*(a(2)-a(1))/(f(a(2))-f(a(1)));
                   a(1)=a(2);
                   a(2)=m;
                   E(i) = (a(2)-a(1))/a(1);
                   iter(i) = i;
                   
                   if abs(f(a(2)))<e2
                       fprintf('%f',a(2));
                       fprintf(' Termination due to convergence of function value');
                       break;
                   end
                   
                   if abs((a(2)-a(1))/a(1))<e1
                       fprintf('%f',a(2));
                       fprintf(' Termination due to convergence of relative approximate error');
                       break;
                   end
                   if i==itr
                   fprintf('%f',a(2));
                   fprintf(' Termination due to limited number of iterations');
                   end
               end      
              
               fplot(f);grid on;title('The Function');ylabel('f(x)');xlabel('x');
               figure;
               plot(iter,E);grid on;title('Relative Error vs No. of iterations');ylabel('Relative Error');xlabel('Number of iterations');
    end
    
elseif c=='Y'
    
    n = input('Enter the order of polynomial: ');
    a = input('Enter the vector of the polynomial coefficients (in brackets) in decreasing order of power of x: ');
    fprintf('Choose the method:\n');
    y = input('1.Muller\n2.Bairstow\n');
    syms f(x);
    f(x)=0;
    for i=1:n+1
        f(x)=f(x)+a(i)*(x^(n-i+1));
    end
    
    if y==1
        g = input('Enter three starting points (in bracket): ');
        fprintf('Enter the stopping criteria: \n');
        i = input('Enter maximum number of iterations ');
        e = input('Enter the max. relative approximate errors in successive iterations: ');
        e1 = input('How much close to zero can function value go: ');
        z = int16.empty(i,0);
        E = double.empty(i,0);
        x1 = g(1);
        x2 = g(2);
        x3 = g(3);
        iter = 0;
        err = 100;
        x = x1:0.01:x3;
        y = f(x);
        flag = 0;
        while((err>e)&&(iter<i))
            if abs(f(x3))<e1
                flag=1;
                break;
            end
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
        fprintf('%f',x3);
        if flag==1
            fprintf(' Termination due to convergence of function value');
        elseif i==iter
            fprintf('Termination due to limited number of iterations');
        else
            fprintf(' Termination due to convergence of relative approximate error');
        end 
         fplot(f,[-50 50]);grid on;title('The Function');ylabel('f(x)');xlabel('x');
         figure;
         plot(z,E);grid on;title('Absolute relative Error vs No. of iterations');ylabel('Error');xlabel('Number of iterations');
        
    elseif y==2
      s3=input('Enter the max. Number of iterations: ');
      s2=input('Enter the max. relative approximate errors in successive iterations: ');
      s1=input('How close can function go to zero: ');
      k=1;
      temp=length(a);
      c=double.empty(temp,0);
      for i=1:temp
          c(i) = a(temp+1-i);
      end
      for i=1:temp
          a(i) = c(i);
      end
      
        while n>1
            r=input('Enter r: ');
            s=input('Enter s: ');
            d=zeros(n+1,1);
            d1=zeros(n+1,1);
            d(1)=100;
            e1=[300];
            e2=[300];
            t=1;
            while (abs(d(1))>s2||abs(d(2))>s2)&&t<s3&&(abs(e1(t))>s1||abs(e2(t))>s1)
            
                d(n+1)=a(n+1);
                d(n)=a(n)+s*d(n+1);
                for i=n-1:-1:1
                    d(i)=a(i)+s*d(i+1)+r*d(i+2);
                end
        
                d1(n)=d(n+1);
                for j=n-1:-1:1
                    d1(j)=d(j+1)+r*d1(j+2)+s*d1(j+1);
                end
                r1=(d1(1)*d(2)-d1(2)*d(1))/(d1(2)^2-d1(3)*d1(1));
                s1=(d1(2)*d(2)-d1(3)*d(1))/(d1(1)*d1(3)-d1(2)^2);
            
                r=r+r1;
                s=s+s1;
                e1=[e1,r1/r*100];
                e2=[e2,s1/s*100];
                t=t+1;
            end
            it=1:t-1;
          
            e1=e1(2:t);
            e2=e2(2:t);
            if t>s3
               fprintf('Termination due to limited number of iterations');
            elseif abs(e1(i))<=s1||abs(e1(i))<=s1
                fprintf('Termination due to convergence of function value');
            elseif (abs(d(1))<=s2||abs(d(2))<=s2)
                fprintf('Termination due to convergence of relative approximate error');
            end
            figure;
            plot(it,e1);grid on;title('Relative Error vs No. of iterations');xlabel('Number of iterations');ylabel('Relative Error');
            figure;
            plot(it,e2);grid on;title('Relative Error vs No. of iterations');xlabel('Number of iterations');ylabel('Relative Error');
            
            if k==1
                root1=0.5*(s+(s^2+4*r)^0.5);
                root2=0.5*(s-(s^2+4*r)^0.5);
                root=[root1,root2];
                disp(root);
            else
                root1=0.5*(s+(s^2+4*r)^0.5);
                root2=0.5*(s-(s^2+4*r)^0.5);
                root=[root,root1,root2];
                disp(root);
            end
           
            d(n+1)=a(n+1);
            d(n)=a(n)+s*d(n+1);
            for i=n-1:-1:1
                d(i)=a(i)+s*d(i+1)+r*d(i+2);
            end
            a=d(3:n+1);
            n=n-2;
            k=k+1;
            if n==1
            root=[root,-a(1)/a(2)];
            disp(root);
            end
            
        end
    end 
end
end
