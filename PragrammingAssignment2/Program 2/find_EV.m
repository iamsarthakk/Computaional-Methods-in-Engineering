function find_EV
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fid = fopen(s,'r');
X = fscanf(fid, '%f', inf);
fclose(fid);
n = X(1);
m = size(X);
error = X(m);
A = zeros(n,n,'double');
for i=1:n
    for j=1:n
        A(i,j) = X((i-1)*(n-1) + (i+j));
    end
end
 
choice = input('Wants only largest eigen value(type L) or all eigen values(type A)?\n','s'); 

if choice == 'L'
     z = double.empty(n,0);
     z(1) = 1;
     for i=2:n
         z(i) = 0;
     end
     e = 100;
     y = transpose(z);
     z = A*y;
     lambda = transpose(y)*z;
     while(e>error)
         div = sqrt(transpose(z)*z);
         y = z./div;
         z = A*y;
         lambda1 = transpose(y)*z;
        
         e = abs((lambda-lambda1)/lambda*100);
         lambda = lambda1;
     end
    
        fileid1 = fopen("Power-Method.txt" , 'wt');
        fprintf(fileid1,'Largest Eigen Value: \n');
        fprintf(fileid1,'lambda = ');
        formatspec1 = '%0.5f\n';
        fprintf(fileid1,formatspec1,lambda);
        fclose(fileid1);
    
elseif choice == 'A'
    Q = zeros(n,n,'double');
    
    e = 100;
    iter = 1;
    while e>error
    %div = 0;
    a = A(:,1);
    a = a./sqrt(transpose(a)*a);
    Q(:,1) = a;
    
    for i=2:n
        z=A(:,i);
        for k=2:i
            z = z-transpose(A(:,i))*Q(:,k-1).*Q(:,k-1);
        end
        Q(:,i) = z./sqrt(transpose(z)*z);
    end
    
    for i=1:n
        for j=1:n
            if i>j
                R(i,j) = 0;
            else
                R(i,j) = transpose(Q(:,i))*A(:,j);
            end
        end
    end
     
     for i=1:n
        for j=1:n
            if i==j
                lambda(i) = R(i,j);
            end
        end
     end
     if iter>1
    for i=1:n
        for j=1:n
            if i==j
                er(i) = abs((lambda(i)-l(i))/l(i)*100);
            end
        end
    end
    e = max(er);
     end
    A = R*Q;
    l = lambda;
    iter=iter+1;
    
    end
   
        fileid1 = fopen("Q-R_Decomp.txt" , 'wt');
        fprintf(fileid1,'All Eigen Values: \n');
        fprintf(fileid1,'lambda =  \n');
        formatspec1 = '%0.5f\n';
        fprintf(fileid1,formatspec1,lambda);
        fclose(fileid1);
  
else
    fprintf('Wrong Choice');
end
end