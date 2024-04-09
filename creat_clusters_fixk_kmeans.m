function [clusters] =creat_clusters_fixk_kmeans(data,H)
[n,~]=size(data);
clusters=zeros(n,H);
maxk=floor(sqrt(n));
k=min(maxk,50);
for i=1:H
%     intial_center=randperm(n,k)';
%     intial_center=data(intial_center,:);
%     clusters(:,i)=kmeans(data,k,'emptyaction','singleton','start',intial_center);

%     clusters(:,i)=kmeans(data,k,'emptyaction','singleton','distance','cosine');
   clusters(:,i)=kmeans(data,k,'emptyaction','singleton');


%         clusters(:,i)=SpectralCluster(data,k);

end

end