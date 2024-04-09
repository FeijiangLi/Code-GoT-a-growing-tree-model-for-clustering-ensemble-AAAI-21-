load Iris.mat
gt=data(:,end);
H=51;
k=length(unique(gt));
data_feature=data(:,1:end-1);
% data_feature=predata(data_feature);
[clusterings] =creat_clusters_fixk_kmeans(data_feature,H);
result=GoT(clusterings,k );
[ac,ARI,NMI]=evaluate2(result,gt,k)