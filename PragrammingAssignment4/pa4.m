syms p(x);
syms q(x);
syms r(x);
syms s(x);
p(x)=input('Input function p(x): ');
q(x)=input('Input function q(x): ');
r(x)=input('Input function r(x): ');
s(x)=input('Input function s(x): ');
h=input('Input grid size h: ');
a1=input('Left limit: ');
a2=input('Right limit: ');
op=input('Choose any one from following\n1 2nd order backward difference\n2 Central difference with ghost node\n');
n=(a2-a1)/h;
y1=input('Derivative: ');
l=zeros(n,1);
d=zeros(n,1);
u=zeros(n,1);
b=zeros(n,1);
a11=input('Boundary Value:');
for j=1:n
    l(j)=p(a1+j*h)/h^2-q(a1+j*h)/(2*h);
    d(j)=-2*p(a1+j*h)/h^2+r(a1+j*h);
    u(j)=p(a1+j*h)/h^2+q(a1+j*h)/(2*h);
    b(j)=s(a1+j*h);
end
b(1)=b(1)-l(1)*a11;
l(1)=0;
if op==1
    b(n-1)=b(n-1)-2*u(n-1)/3*y1*h;
    l(1)=0;
    l(n-1)=l(n-1)-u(n-1)/3;
    d(n-1)=d(n-1)+4/3*u(n-1);
    u(n-1)=0;
    %fun1(l(1:n-1),d(1:n-1),u(1:n-1),b(1:n-1));
    l=l(1:n-1);
    d=d(1:n-1);
    u=u(1:n-1);
    b=b(1:n-1);
    [n1,qp]=size(d);
al=[d(1)];
be=[b(1)];
for j=2:n1
    al1=d(j)-l(j)*u(j-1)/al(j-1);
    be1=b(j)-l(j)*be(j-1)/al(j-1);
    al=[al,al1];
    be=[be,be1];
end
x=[be(n1)/al(n1)];
for j=n1-1:-1:1
    x1=(be(j)-u(j)*x(n1-j))/al(j);
    x=[x,x1];
end
y=[x(n1)];
for j=2:n1
    y=[y,x(n1-j+1)];
end
y=[y,(2*h*y1+4*y(n-1)-y(n-2))/3]
%fclose(a);
a=fopen('sol1.txt','wt');
fprintf(a,'%f\n',y);
fclose(a);
else
    b(n)=b(n)-2*y1*h*u(n);
    l(n)=l(n)+u(n);
    %fun1(l,d,u,b);
    u(n)=0;
    [n1,qp]=size(d);
al=[d(1)];
be=[b(1)];
for j=2:n1
    al1=d(j)-l(j)*u(j-1)/al(j-1);
    be1=b(j)-l(j)*be(j-1)/al(j-1);
    al=[al,al1];
    be=[be,be1];
end
x=[be(n1)/al(n1)];
for j=n1-1:-1:1
    x1=(be(j)-u(j)*x(n1-j))/al(j);
    x=[x,x1];
end
y=[x(n1)];
for j=2:n1
    y=[y,x(n1-j+1)];
end
y
%fclose(a);
a=fopen('sol2.txt','wt');
fprintf(a,'%f\n',y);
fclose(a);
end

    
    

	