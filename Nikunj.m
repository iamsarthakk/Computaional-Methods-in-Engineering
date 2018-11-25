function fun5(c)
a = fopen(c,'rt');
n=fgets(a);
n=str2num(n);
b=fgets(a);
b=split(b);
b=str2double(b);
b=[b(1:n)'];
for j=2:n
    d=fgets(a);
    d=split(d);
    d=str2double(d);
    b=[b;d(1:n)'];
end
d=fgets(a);
d=split(d);
d=str2double(d);
e=b;
t1=fopen('log1.txt','wt');
for j=1:n-1
    p1=find(abs(e(j:n,j:n))==max(max(abs(e(j:n,j:n)))),1);
    p=mod(p1,n-j+1);
    if p==0
        p=p1/(n-j+1);
    end
    if p~=0
        q=(p1+n-j+1-p)/(n-j+1);
    else
        q=p1/(n-j+1);
    end
    p=p+j-1;
    q=q+j-1;
    if p~=j && q~=j
        t=e(j,:);
        e(j,:)=e(p,:);
        e(p,:)=t;
        t=b(j,:);
        b(j,:)=b(p,:);
        b(p,:)=t;
        fprintf(t1,"%s %d %s %d\n","Exchange row",j,"with row",p);
        t=e(:,j);
        e(:,j)=e(:,q);
        e(:,q)=t;
        t=b(:,j);
        b(:,j)=b(:,q);
        b(:,q)=t;
        fprintf(t1,"%s %d %s %d\n","Exchange column",j,"with column",q);
    end
       for k=j+1:n
          fac=e(k,j)/e(j,j);
          e(k,:)=e(k,:)-fac*e(j,:);
       end
end
l=zeros(n);
u=zeros(n);
for k=1:n
    u(k,k)=1;
    for i=k:n
        r=b(i,k);
        for m=1:k-1
            r=r-l(i,m)*u(m,k);
        end
        l(i,k)=r;
    end
    for j=k+1:n
        r=b(k,j);
        for m=1:k-1
            r=r-l(k,m)*u(m,j);
        end
        u(k,j)=r/l(k,k);
    end
end
fprintf(t1,"%s","Matrix L\n");
for j=1:n
    fprintf(t1,"%f ",l(j,:));
    fprintf(t1,"\n");
end
fprintf(t1,"%s","Matrix U\n");
for j=1:n
    fprintf(t1,"%f ",u(j,:));
    fprintf(t1,"\n");
end
fclose(t1);
fclose(a);