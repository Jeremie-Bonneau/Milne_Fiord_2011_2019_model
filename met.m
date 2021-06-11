%% Data     
t0=datenum(2011,05,01);
tf=datenum(2019,08,01);
fields={'Time', 'temp_1m', 'temp_2m', 'pressure', 'RH', 'radiation', 'wind', 'gust','wind_dir', 'battery','snow'};


[m,time]=met_info(fields{2},t0,tf,'hour');
data=nan(length(m'),11);
data(:,1)=time';
data(:,2)=m';

for i=3:11
    [m,~]=met_info(fields{i},t0,tf,'hour');
    data(:,i)=m';
end
time=dat(data(:,1));

%%
 rootfolder = ['C:\Users\jerem\OneDrive\Documents\UBC\Archive'];
 outputdir = [rootfolder];
 fid = fopen([outputdir '/' 'Purple_Valley_met_2011_2019'],'w'); 
fprintf(fid,'time;temp_1m;temp_2m;pressure;RH;radiation;wind;gust;wind_dir;battery;snow\n');


for i=1:length(data(:,1))
    fprintf(fid,'%8s;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f\n', ...
        time(i),data(i,2),data(i,3),data(i,4),data(i,5),data(i,6),data(i,7),data(i,8),data(i,9),data(i,10),data(i,11));    
end
%format = ['%s %8.3f %8.1f %8.5f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n'];
%dlmwrite([outputdir '/' 'Purple_Valley_met_2011_2019'], data,'delimiter','\t','precision','%8.4f', '-append');
fclose(fid);