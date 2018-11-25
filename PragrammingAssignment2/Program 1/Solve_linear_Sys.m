function Solve_linear_Sys
c = input('Choose one option:\nA. Solve a System of Equation\nB. Perform a LU decomposition\nC. Perform a Matrix Inversion\n','s');

switch c
    case 'A', choice = input('Is the system is tri-diagonal? (Y/N)\n','s');
        if choice=='Y'
            [filename,filepath] = uigetfile('.txt','File Selector');
            s = fullfile(filepath,filename);
            fid = fopen(s , 'r');
            X = fscanf(fid, '%f', inf);
            fclose(fid);
            n = X(1);
            m = size(X);
            x=0;
            for i=2:n:m
                l(mod(i-1,n)+x) = X(i);
                d(mod(i-1,n)+x) = X(i+1);
                u(mod(i-1,n)+x) = X(i+2);
                b(mod(i-1,n)+x) = X(i+3);
                x=x+1;
            end
            alpha(1) = d(1);
            beta(1) = b(1);
            for i=2:n
                alpha(i) = d(i)-l(i)*u(i-1)/alpha(i-1);
                beta(i) = b(i)-l(i)*beta(i-1)/alpha(i-1);
            end
            x(n) = beta(n)/alpha(n);
            for i=n-1:-1:1
                x(i) = (beta(i)-u(i)*x(i+1))/alpha(i);
            end
          
        fileid1 = fopen("Thomas_Output.txt" , 'wt');
        fprintf(fileid1,'Thomas Algorithm Solutions: \n');
        fprintf(fileid1,'X =\n');
        formatspec1 = '%0.5f\n';
        fprintf(fileid1,formatspec1,x);
        fclose(fileid1);
        
        elseif choice == 'N'
            [filename,filepath] = uigetfile('.txt','File Selector');
            s = fullfile(filepath,filename);
            fid = fopen(s , 'r');
            X = fscanf(fid, '%f', inf);
            fclose(fid);
            n = X(1);
            A = zeros(n,n+1,'double');
           for i=1:n
               for j=1:n+1
                   A(i,j) = X((i-1)*(n-1) + (i+j)+i-1);
               end
           end
          
           r = 1;
        while(r<=n)
            tmp = A(r,r);
            tmp1 = r + 1;
            tmpr = r;
            z = r+1;
            while(z<=n)
                if(abs(A(z,r))>=abs(tmp))
                    tmp = A(z,r);
                    tmpr = z;
                end
                z = z+1;
            end
            if(tmpr~=r)
                temp = A(r,:);
                A(r,:) = A(tmpr,:);
                A(tmpr,:)= temp;
            end
            
            while(tmp1<=n)
                if(A(tmp1,r)~=0)
                    A(tmp1,:) = A(tmp1,:) - A(tmp1,r)*A(r,:)/tmp;
                end
                tmp1 = tmp1 + 1;
            end
            
            r = r + 1;
        end
        x(n) = A(n,n+1)/A(n,n);
        j = n-1;
        while(j>0)
        y = n;
        while(y>j)
            A(j,n+1) = A(j,n+1) - x(y)*A(j,y);
            y = y - 1;
        end
        x(j) = A(j,n+1)/A(j,j);
        j = j - 1;
        end
    
        fileid1 = fopen("Gauss_Elimination_Output.txt" , 'wt');
        fprintf(fileid1,'Solutions using Gauss Elimination: \n');
        fprintf(fileid1,'X = \n');
        formatspec1 = '%0.5f\n';
        fprintf(fileid1,formatspec1,x);
        fclose(fileid1);
        end
        
    case 'B', choice = input('Is the matrix is symmetric and positive definite? [Y/N]\n','s');
        if choice=='Y'
            [filename,filepath] = uigetfile('.txt','File Selector');
            s = fullfile(filepath,filename);
            fid = fopen(s , 'r');
            X = fscanf(fid, '%f', inf);
            fclose(fid);
            n = X(1);
            C = zeros(n,n+1,'double');
            A = zeros(n,n,'double');
            l = zeros(n,n,'double');
           for i=1:n
               for j=1:n+1
                      C(i,j) = X((i-1)*(n-1) + (i+j)+i-1);
               end
           end
           
           for i=1:n
               for j=1:n
                   A(i,j) = C(i,j);
               end
           end
               
           fout=fopen('Cholesky_output.txt','wt');
           b = A;
           for j=1:n
               p1 = find(abs(b(j:n,j:n)) == max(max(abs(b(j:n,j:n)))),1);
               p = mod(p1,n-j+1);
           end
           if p~=0
               q = (p1+n-j+1-p)/(n-j+1);
           else
               q = p1/(n-j+1);
           end
           p = p+j-1;
           q = q+j-1;
           if p~=j && q~=j
               t = b(j,:);
               b(j,:) = b(p,:);
               b(p,:) = t;
               t = A(j,:);
               A(j,:) = A(p,:);
               A(p,:) = t;
               fprintf(fout,"%s %d %s %d\n","Exchange row",j,"with row",p);
               t = b(:,j);
               b(:,q) = t;
               t = A(:,j);
               A(:,j) = A(:,q);
               A(:,q) = t;
               fprintf(fout,"%s %d %s %d\n","Exchange column",j,"with column",q);
           end
           for k=j+1:n
               fac = b(k,j)/b(j,j);
               b(k,:) = b(k,:)-fac*b(j,:);
           end
       
               
           x=0;
           for i=1:n
               for j=1:n
                   if i==j
                      for k=1:j-1
                          x=x+l(j,k)*l(j,k);
                      end
                      l(i,j) = sqrt(A(i,j)-x);
                      x = 0;
                   elseif i>j
                       for k=1:j-1
                           x = x+l(i,k)*l(j,k);
                       end
                       l(i,j) = (A(i,j) - x)/l(j,j);
                       x = 0;
                   end
                   
               end       
           end
           
           fprintf(fout,"Matrix L:\n");
           for i=1:n
           fprintf(fout,"%f ",l(i,:));
           fprintf(fout,"\n");
           end
           
           fclose(fout);
            
        elseif choice == 'N'
            [filename,filepath] = uigetfile('.txt','File Selector');
            s = fullfile(filepath,filename);
            fid = fopen(s , 'r');
            X = fscanf(fid, '%f', inf);
            fclose(fid);
            n = X(1);
            C = zeros(n,n,'double');
            A = zeros(n,n,'double');
            b = zeros(n,n,'double');
            l = zeros(n,n,'double');
            u = zeros(n,n,'double');
           for i=1:n
               for j=1:n+1
                   C(i,j) = X((i-1)*(n-1) + (i+j)+i-1);
               end
              
           end
          
           for i=1:n
               for j=1:n
                   A(i,j) = C(i,j);
               end
           end
           b = A;
           for j=1:n
               p1 = find(abs(b(j:n,j:n)) == max(max(abs(b(j:n,j:n)))),1);
               p = mod(p1,n-j+1);
           end
           if p~=0
               q = (p1+n-j+1-p)/(n-j+1);
           else
               q = p1/(n-j+1);
           end
           p = p+j-1;
           q = q+j-1;
           if p~=j && q~=j
               t = b(j,:);
               b(j,:) = b(p,:);
               b(p,:) = t;
               t = A(j,:);
               A(j,:) = A(p,:);
               A(p,:) = t;
               fprintf(fout,"%s %d %s %d\n","Exchange row",j,"with row",p);
               t = b(:,j);
               b(:,q) = t;
               t = A(:,j);
               A(:,j) = A(:,q);
               A(:,q) = t;
               fprintf(fout,"%s %d %s %d\n","Exchange column",j,"with column",q);
           end
           for k=j+1:n
               fac = b(k,j)/b(j,j);
               b(k,:) = b(k,:)-fac*b(j,:);
           end
            fprintf('Choose the method of decomposition:\n');
            c1 = input('1.Doolittle\n2.Crout\n');
            if c1==1
                for i=1:n
                    for j=1:n
                       if i == j
                            l(i,j) = 1;
                        end
                     end
                end
                x = 0;
                for i=1:n
                    for j=1:n
                        if i<=j
                            for k=1:i-1
                                x = x+l(i,k)*u(k,j);
                            end
                            u(i,j) = A(i,j) - x;
                            x = 0;
                        else
                            for k=1:j-1
                                x = x+l(i,k)*u(k,j);
                            end
                            l(i,j) = (A(i,j)-x)/u(j,j);
                            x = 0;
                        end
                    end
                   
                end 
                fout=fopen('Doolittle.txt','wt');
                fprintf(fout,"Matrix L:\n");
                for j=1:n
                    fprintf(fout,"%f ",l(j,:));
                    fprintf(fout,"\n");
                end
                fprintf(fout,"Matrix U:\n");
                for j=1:n
                fprintf(fout,"%f ",u(j,:));
                fprintf(fout,"\n");
                end
                
               
                fclose(fout);
        
                
            elseif c1 == 2
                 for i=1:n
                    for j=1:n
                       if i == j
                            u(i,j) = 1;
                        end
                     end
                end
                x = 0;
                for i=1:n
                    for j=1:n
                        if j<=i
                            for k=1:j-1
                                x = x+l(i,k)*u(k,j);
                            end
                            l(i,j) = A(i,j) - x;
                            x = 0;
                        else
                            for k=1:i-1
                                x = x+l(i,k)*u(k,j);
                            end
                            u(i,j) = (A(i,j)-x)/l(i,i);
                            x = 0;
                        end
                    end
                   
                end  
                 
            end
                fout = fopen("Crout.txt",'wt');
                fprintf(fout,"Matrix L:\n");
                for j=1:n
                    fprintf(fout,"%f ",l(j,:));
                    fprintf(fout,"\n");
                end
                fprintf(fout,"Matrix U:\n");
                for j=1:n
                fprintf(fout,"%f ",u(j,:));
                fprintf(fout,"\n");
                end
                
                fclose(fout);
                
        end
        
    case 'C',[filename,filepath] = uigetfile('.txt','File Selector');
            s = fullfile(filepath,filename);
            fid = fopen(s , 'r');
            X = fscanf(fid, '%f', inf);
            fclose(fid);
            n = X(1);
            A = zeros(n,2*n,'double');
           for i=1:n
               for j=1:n
                   A(i,j) = X((i-1)*(n-1) + (i+j));
               end
           end
           for i=1:n
               for j=n+1:2*n
                   if j-i==n
                       A(i,j) = 1;
                   end
               end
           end
         
           for k=1:n
               A(k,:) = A(k,:)./A(k,k);
               
               for i=1:n
                   if i~=k
                       A(i,:)=A(i,:)-A(i,k)*A(k,:);
                   end
               end
           end
           for i=1:n
               for j=n+1:2*n
                   Inv(i,j-n) = A(i,j);
               end
           end
       
        fileid1 = fopen('Gauss-Jorden.txt', 'wt');
        fprintf(fileid1,'Inverse of matrix using Gauss Jorden Elimination: \n');
        for i=1:n
        fprintf(fileid1,"%f ",Inv(i,:));
        fprintf(fileid1,"\n");
        end
        disp(A);
        disp(Inv);
        fclose(fileid1);
       
end
end