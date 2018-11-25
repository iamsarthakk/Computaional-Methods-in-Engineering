f1=fopen('A5Q7.txt','r');
A=fscanf(f1,'%f');
[n,m]=size(A)
x1=[];
y1=[];
for i=1:n/2
    x1=[x1,A(2*i-1,1)];
    y1=[y1,A(2*i,1)];
end
n=n/2;
fprintf('choose from following option\n');
fprintf('A. Fit a least square polynomial\n')
fprintf('B. Fit a Lagrange Interpolation Polynomial\n');
fprintf('C. Fit Cubic splines\n');
Z=input('input options\n','s')
if Z=='A'
    m=input('chose order of the polynomial\n');
    b=input('weather the intercept is zero\n','s');
    C=zeros(m,m+1);
    for i=1:m
        for j=1:m+1
            if j>m
                for k=1:n
                    C(i,j)=C(i,j)+(x1(1,k)^(i-1))*(y1(1,k))
                end
            else
                for k=1:n
                    C(i,j)=C(i,j)+(x1(1,k)^(i-1))*(x1(1,k)^(j-1))
                end
            end
        end
        c=inv(C(1:m,1:m))*C(:,m+1);
    end
elseif Z=='B'
    A=zeros(n,n);
    A(:,1)=y1;
    for i=2:n
        for j=1:n-i+1
            x1(1,j+i-1)
            x1(1,1)
            A(j,i)=(A(j,i-1)-A(j+1,i-1))/(-x1(1,j+i-1)+x1(1,j))
        end
    end
    c=A(1,:);
else
    fprintf('Select a given option\n');
    fprintf('1. Natural Spline\n2. Not a knot');
    fprintf('3. Periodic\n4. Clamped Spline');
    k=input('Enter your option');
    switch(k)
        case(1)
            
    
end
    
            
        
    
    