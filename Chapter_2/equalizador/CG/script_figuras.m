files = dir('*.fig');
nfigs = length(files); 
 
files = {files.name};
 
dir_sav='last';
mkdir(dir_sav);
 
for i=1:nfigs
    fprintf(['Abrindo: ' files{i} '\n'])
    open(files{i});
    caxis([-2.2 1]); 
    saveas(gcf,[dir_sav '/' files{i}]);
    close(gcf);
end
 
cd(dir_sav);
 

addpath('..');
for i=1:nfigs
    fprintf(['Formatando: ' files{i} '\n'])
    open(files{i});
    formatFig_mod(files{i}(1:end-4)); 
    close(gcf);
end
 
 

 

system('mv *.eps jon04/OneDrive/Tese/Codes/Cap 2/equalizador/CG/figs')
%system('mv *.eps /nfs/home/marcele.kuhfuss/Dropbox')
