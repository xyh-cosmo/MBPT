d=load('xyh_1020_tk.dat');
N=length(d);
fp=fopen('xyh.dat','w');
for i=1:N
    fprintf(fp,'%10.5E %10.5E %10.5E %10.5E %10.5E\n',d(i,1)/0.71,-d(i,4),-d(i,3),0,-d(i,8));
end
fclose(fp);
