function sspline()
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s, 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I(1) - 1;
z = double.empty(n+1,0);
y = double.empty(n+1,0);
j = 2;
j1 = 3;
while(j<=(2*n+2))
    z(j/2) = I(j);
    j = j + 2;
end
while(j1<=(2*n+3))
   y((j1-1)/2) = I(j1); 
   j1 = j1 + 2;
end
v = double.empty(n+1,0);
v(1) = 0;
v(n+1) = 0;
u = double.empty(n,0);
dd = double.empty(n,0);
h = double.empty(n,0);
for k = 1:n
    dd(k) = (y(k+1) - y(k))/(z(k+1) - z(k));
    h(k) = (z(k+1) - z(k));
end
A = zeros(n-1,n-1,'double');
for j3 = 2:n
    A(j3-1,j3-1) = 2*(h(j3-1) + h(j3));
    if(j3>2)
    A(j3-1,j3-2) = h(j3-1);
    end
    if(j3<n)
    A(j3-1,j3) = h(j3);
    end
end
B = double.empty(n-1,0);
for j4 = 2:n
    B(j4-1) = 6*dd(j4) - 6*dd(j4-1);
end
tmp2 = LU_crout(A,B,n-1);
for k1 = 2:n
    v(k1) = tmp2(k1-1);
end
for k2 = 1:n
    u(k2) = dd(k2) - (h(k2)*(v(k2+1) + 2*v(k2)))/6;
end
syms x;
syms f(x);
syms sum;
sum = 0;
for k3 = 1:n
    sum = sum + (y(k3) + (x - z(k3))*u(k3) + (x - z(k3))^2*(v(k3)/2 )+ ((x - z(k3))^3*(v(k3+1) - v(k3))/(6*h(k3))))*(heaviside(x - z(k3)) - heaviside(x - z(k3 + 1)));
end
    f(x) = sum;
    x = z(1):0.001:z(n+1);
    y1 = f(x);
    figure
    plot(x,y1);
    
hold on
for tmp5 = 1:n+1
    plot(z(tmp5),y(tmp5),'ro');
    hold on
end
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
    