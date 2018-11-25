[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s, 'r');
tline = fgetl(fileid);
ip = cell(4,1);
tmp = 1;
while ischar(tline)
    ip(tmp) = {tline};
    tline = fgetl(fileid);
    tmp = tmp + 1;
end
fclose(fileid);
str = char(ip(1));
f = inline(str);
str1 = char(ip(2));
str2 = strsplit(str1,',');
a = str2double(str2(1));
b = str2double(str2(2));
str3 = char(ip(3));
err = str2double(str3);
disp(err);
Z = input(['Choose method.',...
    '\n1. Romberg Method',... 
    '\n2. Gauss Legendre Quadrature: ']); 
switch Z
    case 1
         h = (b - a);
         i = 1;
         e = 100;   
         sum(i,i) = h*(f(a) + f(b))/2;
         while(abs(e) >= err)
                 h = h/2;
                 i = i+1;
                 tmp = 0;
                 tmp1 = a+h;
                 while(tmp1<=b-h)
                    tmp = tmp + f(tmp1);
                    tmp1 = tmp1 + h;
                 end
                 sum(i,1) = ((f(a) + f(b))/2 + tmp)*h;
                 e = (sum(i,1) - sum(i-1,i-1))/sum(i,1)*100;
                 if(abs(e)>err)
                    for tmp2 = 1:i-1
                        i2hk = sum(i-1,tmp2);
                        ihk = sum(i,tmp2);
                        ihkp2 = (2^(2*tmp2)*ihk - i2hk)/(2^(2*tmp2) - 1);
                        sum(i,tmp2+1) = ihkp2;
                        e = ((sum(i,tmp2+1) - sum(i,tmp2))/sum(i,tmp2+1))*100;
                        if(abs(e)>err)
                            continue;
                        elseif(abs(e)<=err)
                            break;
                        end
                    end
                 end
         end
    dim = size(sum);
    row = dim(1);
    col = dim(2);
    fprintf('%0.8f',sum(row-1,col));
    disp(e);
    case 2
        tmp1 = 1;
        z = vpasolve(legendreP(tmp1,x) == 0);
        disp(z);
end