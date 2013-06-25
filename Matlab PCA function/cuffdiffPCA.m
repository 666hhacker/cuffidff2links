function [ ] = cuffdiffPCA( fname )
%Parse a PCA input file and pass the data to pcafigure() function

z = importdata(fname);
header = z.textdata(1,2:end);
pcafigurex(z.data,header);
end

function [  ] = pcafigurex( mat,colLabel )
%Perform Principal Component Analysis on the input data and generate a
%3D scatter plot as result.
%   Input should be matrix and column labels. For a gene expression matrix,
%rows are genes and columns are experiments. Experiments belonging to 
%the same comparison group should have the same labels, which give their 
%identical color in the final figure. 
%
%usage:
%     pcafigure(matrix,columnlabel);
%
%The most possible error for running this function is OUT OF MEMORY. When
%this happens, try to reduce the size of your matrix
%
%Author:
% Qiaonan Duan, Ma'ayan lab, Icahn School of Medicine at Mount Sinai
% *Adapted from Simon Gordonov's PCA script.

expressions = mat;
colorLabels = colLabel;
cmap = [[240 163 255];[0 117 220]; [153 63 0] ;[255 80 5]; [255 204 153] ;[116 10 255] ;[43 206 72] ;[0 0 0]; [128 128 128];[148 255 181];[143 124 0];[157 204 0];[194 0 136];[0 51 128];[255 164 5];[255 168 187];[66 102 0];[255 0 16];[94 241 242];[0 153 143];[224 255 102];[116 10 255];[153 0 0];[255 255 128];[255 255 0];[255 80 5]]./255;

   
[colorLabels, idx] = sort(colorLabels);
expressions = expressions(:,idx);

expressions = expressions';
[pc, score, pcvars] = princomp(expressions);

x = zscore(score(:,1));
y = zscore(score(:,2));
z = zscore(score(:,3));

figure
gscat = gscatter(x,y,colorLabels',cmap,'.',18,'on');
uniqueLabel  = unique(colorLabels);
for i = 1:numel(uniqueLabel)
    eachLabelIndex = ismember(colorLabels,uniqueLabel(i));
    set(gscat(i),'ZData', z(eachLabelIndex));
end
legendHandle = legend('location','bestoutside');
set(legendHandle,'fontsize',11);
set(gca,'Xgrid','on','ygrid','on','fontsize',18);
box on;
xlabel(['PC-1 (',num2str(round(pcvars(1)/sum(pcvars)*100)),'%)'],'fontsize',18);
ylabel(['PC-2 (',num2str(round(pcvars(2)/sum(pcvars)*100)),'%)'],'fontsize',18);
zlabel(['PC-3 (',num2str(round(pcvars(3)/sum(pcvars)*100)),'%)'],'fontsize',18);

set(gcf,'color',[1 1 1]);

end

