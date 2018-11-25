n=input('Enter Degree');
s3=input('Number of iterations');
s2=input('Enter error');
s1=input('Closeness to zero');
a=input('Enter coefficiets');
k=1;
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
            end
            disp('Quadratic is x^2-s*x-r: ' );
            
           
            
        end