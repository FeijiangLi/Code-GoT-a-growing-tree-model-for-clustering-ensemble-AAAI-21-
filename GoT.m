function [cl] = GoT(clusters,k)
[sim] = simnumber(clusters);
[dens]=cal_density2(clusters);
sim(sim==0)=0.0001;
dis=-log(sim);
sdis=sparse(dis);
% dis=graphallshortestpaths(sdis);
G=graph(sdis);
dis=distances(G);
[lden] = local_dens(dis,dens);
[~,initial]=sort(lden,'descend');
initial=initial(1:k);
sim=-dis;
[cl] = link(sim,initial,k);




function [sim] = simnumber(clusters) 
[N,M] = size(clusters); % no. of clustering
newE = zeros(N,M);
ucl = unique(clusters(:,1)); % all clusters in i-th clustering
if (max(clusters(:,1) ~= length(ucl)))
    for j = 1:length(ucl)
        newE(clusters(:,1) == ucl(j),1) = j; % re-labelling
    end
end
for i = 2:M
    ucl = unique(clusters(:,i)); % all clusters in i-th clustering
    prevCl = length(unique(newE(:,[1:i-1])));
    for j = 1:length(ucl)
        newE(clusters(:,i) == ucl(j),i) = prevCl + j; % re-labelling
    end
end
no_allcl = max(max(newE));
[n,m]=size(newE);
max_E=max(newE);
min_E=min(newE);
relabel=zeros(n,no_allcl);
for i=1:m
    cl=newE(:,i);
    for j=min_E(i):max_E(i)
        relocat=cl==j;
        relabel(relocat,j)=1;
    end
end
sim=relabel*relabel'./M;



function [dens]=cal_density2(clusters)
[n,l]=size(clusters);
clb = BI_clusters(clusters);
[~,cl]=size(clb);
dens=zeros(1,n);
for i=1:n
    zero_clb=zeros(n,cl);
    cl_now=clb(i,:)==1;
    zero_clb(:,cl_now)=clb(:,cl_now);
    dens(i)=sum(sum(zero_clb)./n)/l;
end




function clb = BI_clusters(cl)
[n,m]=size(cl);
[newE, no_allcl] = relabelCl(cl);
clb=zeros(n,no_allcl);
for i=1:m
    now_cl=newE(:,i);
    for j=min(now_cl):max(now_cl)
        clb(:,j)=now_cl==j;
    end
end


function [lden] = local_dens(dis,dens)
n=length(dens);
lden=zeros(1,n);
maxdis=max(max(dis))+1;

lo_maxdens=find(dens==max(dens));
lden(lo_maxdens(1))=maxdis+1;
if length(lo_maxdens)>1
    dens(lo_maxdens(2:end))=dens(lo_maxdens(2:end))-0.0001;
end

for i=[1:lo_maxdens(1)-1, lo_maxdens(1)+1:n]
    now_dis=dis(i,:);
    now_dis(i)=maxdis+1;
    lo_deni=dens==dens(i);
    
    if sum(lo_deni)>1
        dens(lo_deni)=dens(lo_deni)-0.000001;
        dens(i)=dens(i)+0.000001;
    end
    
    lo= dens>dens(i);
    no_dis=min(now_dis(lo));
    lden(i)=no_dis(1);
end

function [cl] = link(sim,initial,k)
n=size(sim,1);
cl=zeros(n,1);
for i=1:k
    cl(initial(i))=i;
end
while sum(cl==0)>0
    liu=find(cl==0);
    labeled=n-length(liu);
    sim_in=zeros(length(liu),length(labeled));
    for i=1:k
        lo_label_i=cl==i;
        liusimi=sim(liu,lo_label_i);
        sim_in(:,i)=max(liusimi,[],2);
    end
    re_sim_in=sort(sim_in,2,'descend');
    cha_sim_in=abs(re_sim_in(:,1)-re_sim_in(:,2));
    
    if sum(sum(cha_sim_in))==0
        now_label_lo=1:1:length(cha_sim_in);
        cl(liu(now_label_lo))=k+1;
    else
        cha_sim_in=(cha_sim_in-min(cha_sim_in))./(max(cha_sim_in)-min(cha_sim_in));
        th=graythresh(cha_sim_in);
        now_label_lo=find(cha_sim_in>=th);
        if isempty(now_label_lo)
            now_label_lo=1:1:length(cha_sim_in);
        end
        [~,label]=max(sim_in(now_label_lo,:),[],2);
        cl(liu(now_label_lo))=label;
    end
end


function [newE, no_allcl] = relabelCl(E) 

[N,M] = size(E); 
newE = zeros(N,M);

ucl = unique(E(:,1)); 
if (max(E(:,1) ~= length(ucl)))
    for j = 1:length(ucl)
        newE(E(:,1) == ucl(j),1) = j; 
    end
end


for i = 2:M
    ucl = unique(E(:,i)); 
    prevCl = length(unique(newE(:,[1:i-1])));
    for j = 1:length(ucl)
        newE(E(:,i) == ucl(j),i) = prevCl + j; 
    end
end

no_allcl = max(max(newE));