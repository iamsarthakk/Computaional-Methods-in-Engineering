function funct_approx
 [filename,filepath] = uigetfile('.txt','File Selector');
 s = fullfile(filepath,filename);
 fid = fopen(s , 'r');
 X = fscanf(fid, '%f', inf);
 fclose(fid);
 n = size(X);
 x = [];
 y = [];
 for i=1:n/2
     x = [x,X(2*i-1,1)];
     y = [y,X(2*i,1)];
 end
 n = n/2;
 
 choice1 = input('What do you want to do?\nA. Fit a least quare polynomial\nB. Fit a Lagrange Interpolation Polynomial\nC. Fit Cubic splines\n','s');
 
 switch choice1
     case 'A' 
         m = input('What should be the order of polynomial? ');
         c = input('choose whether the intercept should be zero or not (Y or N) ','s');
         C=zeros(m,m+1);
         for i=1:m
             for j=1:m+1
                 if j>m
                     for k=1:n
                         C(i,j)=C(i,j)+(x(1,k)^(i-1))*(y(1,k));
                     end
                 else
                     for k=1:n
                         C(i,j)=C(i,j)+(x(1,k)^(i-1))*(x(1,k)^(j-1));
                     end
                 end
             end
             ans = inv(C(1:m,1:m))*C(:,m+1);
            
         end
         disp(ans);
         plot(x,y,'g*');
         
     case 'B'
         A=zeros(n,n);
         A(:,1)=y;
         for i=2:n
             for j=1:n-i+1
                 x1(1,j+i-1);
                 x1(1,1);
                 A(j,i)=(A(j,i-1)-A(j+1,i-1))/(-x(1,j+i-1)+x1(1,j));
             end
         end
         c=A(1,:);
         
     case 'C', choice2 = input("Which type of Spline?\n1. Natural Spline\n2. Not-a-knot\n3. Periodic\n4. Clamped Spline ");
         switch choice2
             case 1,NaturalCubicSpline(x,y,t);
             case 2,NotAKnotCubicSpline(x, y, t );
             case 3,PeriodicCubicSpline(x, y, t );
             case 4,ClampedCubicSpline( x, y, t ,s0 ,s1 );
         end
 end
end