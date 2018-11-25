function[] = regression()
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s, 'r');
I5 = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I5(1) - 1;
z = double.empty(n+1,0);
y = double.empty(n+1,0);
j = 2;
i = 3;
while(j<=(2*n+2))
    z(j/2) = I5(j);
    j = j + 2;
end
delta = z(1);
delta1 = z(1);
for del = 1:n+1
    if(z(del)>=delta)
        delta = z(del);
    end
    if(z(del)<=delta1)
        delta1 = z(del);
    end
end
sum5 = 0;
while(i<=(2*n + 3))
    y((i-1)/2) = I5(i);
    sum5 = double(sum5 + y((i-1)/2));
    i = i + 2;
end
mean = double(sum5/(n+1));
spread = 0;
for k3 = 1:n+1
    spread = double(spread + (y(k3) - mean)^2);
end
d = input('enter the degree of fitting polynomial: ');
coeff = double.empty(d + 1,0);
val = double.empty(2*d+1,0);
for l = 1:2*d+1
    sum = 0;
    for I = 1:n+1
        sum = sum + z(I)^(l-1);
    end
    val(l) = sum;
end
val1 = double.empty(d+1,0);
for k = 1:d+1
    sum = 0;
    for J = 1:n+1
        sum = sum + z(J)^(k-1)*y(J);
    end
    val1(k) = sum;
end
A = zeros(n,n,'double');
for I1 = 1:d+1
    for J1 = 1:d+1
        A(I1,J1) = val(I1 + J1 - 1);
    end
end
%disp(A);
%disp(val1);
coeff = LU_crout(A,val1,d+1);
%disp(coeff);
syms x;
syms f(x);
syms sum2;
sum2 = 0;
for tmp = 1:d+1
    sum2 = sum2 + coeff(tmp)*x^(tmp-1);
end
f(x) = sum2;
dev = 0;
for j3 = 1:n+1
    dev = dev + double((y(j3) - f(z(j3)))^2);
end
cod = double((spread - dev)/spread);
x = delta1:(0.001):delta;
y1 = f(x);
figure 
plot(x,y1,'b')
hold on
for tmp2 = 1:n+1
    plot(z(tmp2),y(tmp2),'r*')
    hold on
end
[filename,filepath] = uigetfile('.txt','File Selector');
s1= fullfile(filepath,filename);
fileid1 = fopen(s1 , 'w');
fprintf(fileid1,'Regression \r\n');
formatspec2 = '%d degree polynomial \r\n';
fprintf(fileid1,formatspec2,d);
formatspec1 = '%0.3f\r\n';
fprintf(fileid1,formatspec1,coeff);
fprintf(fileid1,'COD = %0.3f\r\n',cod);
end
function[c] = LU_crout(a,y,n)
A = zeros(n,n+1,'double');
for b = 1:n
    for r = 1:n
        A(b,r) = a(b,r);
    end
end
for w = 1:n
    A(w,n+1) = y(w);
end
l = zeros(n,n,'double');
u = zeros(n,n,'double');
i = 1;
while(i<=n)
    j = 1;
    while(j<=n)
        if(i==j)
            u(i,j) = 1;
        end
        if(i<j)
            l(i,j) = 0;
        end
        if(i>j)
            u(i,j) = 0;
        end
        j = j + 1;
    end
    i = i + 1;
end
k = 1;
while(k<=n);
     I = k;
    while(I<=n)
        m = 1;
        while(m<k)
            A(I,k) = A(I,k) - l(I,m)*u(m,k);
            m = m + 1;
        end
        l(I,k) = A(I,k);
        I = I + 1;
    end
     J = k + 1;
    while(J<=n)
        m = 1;
        while(m<k)
            A(k,J) = A(k,J) - l(k,m)*u(m,J);
            m = m + 1;
        end
        u(k,J) = A(k,J)/l(k,k);
        J = J+1;
    end
    k = k + 1;
end
z(1) = double(A(1,n+1)/l(1,1));
j = 2;
while(j<=n)
    y = 1;
    while(y<j)
        A(j,n+1) = A(j,n+1) - z(y)*l(j,y);
        y = y + 1;
    end
    z(j) = double(A(j,n+1)/l(j,j));
    j = j + 1;
end
x(n) = double(z(n)/u(n,n));
j = n-1;
while(j>0)
y = n;
while(y>j)
    z(j) = z(j) - x(y)*u(j,y);
    y = y - 1;
end
x(j) = double(z(j)/u(j,j));
j = j - 1;
end
c = trans_pose(x,1,n);
end
function[c] = trans_pose(a,m,n)
c = zeros(n,m,'double');
for i = 1:m
    for j = 1:n
            c(j,i) = a(i,j);
    end
end
end