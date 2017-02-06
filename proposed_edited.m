clc;
clear all;
close all;
 
% Simulation Parameters
 
%Field dimensions in meters
xm = 200;
ym = 200;
 
%initial x and y coordinates of the Sink
sink.x = 0;
sink.y = ym*0.5;
 
%Number of Nodes in the field
n = 600;
 
%Optimal Election Probability of a node
%to become cluster head
p_opt = 0.1;
 
L = 4000;   % packet length
 
%Energy Model
Eo = 0.5;                      %Initial Energy
E_elec=50*0.000000001;         % energy consumed by radio electronics in transmit/receive mode(J/bit)
E_fs=10*0.000000000001;        %energy consumed by the power amplifier on the free space model(J/bit/m2)
E_mp=0.0013*0.000000000001;    %energy consumed by the power amplifier on the multi path model(J/bit/m4)
E_DA=5*0.000000001;            % energy consumed for data aggregation(J/bit/signal)
 
INFINITY = 999999999999999;
 
%maximum number of rounds
rmax=15000;
 
%threshold distance
do= 70;
 
m = 0.75;
a = 2;
m0 = 0.525;
b = 2.5;
m1 = 0.225;
u = 3;
 
c=0.02;
z=0.71;
T_absolute= z*Eo;
 
cr=5;                      % compression ratio
 
S = struct;
S.xd = zeros(n,1);
S.yd = zeros(n,1);
S.E = zeros(n,1);
S.type = zeros(n,1);
S.G = zeros(n,1);
S.CH = zeros(n,1);
S.D = zeros(n,1);
S.T = zeros(n,1);
S.min_dis_cluster= zeros(n,1);
S.dis_to_cluster= zeros(n,1);
 
C = struct;
C.id= zeros(25,1);
C.xd= zeros(25,1);
C.yd= zeros(25,1);
C.E= zeros(25,1);
C.member= zeros(25,1);
 
R = struct;
R.id= zeros(100,1);
R.xd= zeros(100,1);
R.yd= zeros(100,1);
 
p = zeros(n,1);
E = zeros(n,1);
DEAD = zeros(rmax,1);
RES = zeros(rmax,1);
 
% Creation of four-level heterogeneous network 
% Normal Nodes
for i=1:1:240
    S(i).xd = xm*rand(1,1);
    S(i).yd = ym*rand(1,1);
    S(i).E = Eo;
    S(i).type = 'N';
    S(i).G = 1;
    S(i).CH = 0;
    S(i).D = 0;
    S(i).R = 0;
    S(i).min_dis_cluster=0;
    S(i).dis_to_cluster=0;
    figure(3);
    plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 3, 'MarkerFaceColor', 'g');
    hold on;
end
% Advanced Nodes
for i=241:1:420
    S(i).xd = xm*rand(1,1);
    S(i).yd = ym*rand(1,1);
    S(i).E = Eo*(1+a);
    S(i).type = 'A';
    S(i).G = 1;
    S(i).CH = 0;
    S(i).D = 0;
    S(i).R = 0;
    S(i).min_dis_cluster=0;
    S(i).dis_to_cluster=0;
    figure(3);
    plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 3, 'MarkerFaceColor', 'b');
    hold on;
end
% Super Nodes
for i=421:1:546
    S(i).xd = xm*rand(1,1);
    S(i).yd = ym*rand(1,1);
    S(i).E = Eo*(1+b);
    S(i).type = 'S';
    S(i).G = 1;
    S(i).CH = 0;
    S(i).D = 0;
    S(i).R = 0;
    S(i).min_dis_cluster=0;
    S(i).dis_to_cluster=0;
    figure(3);
    plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
    hold on;
end
% Ultra-Super Nodes
for i=547:1:600
    S(i).xd = xm*rand(1,1);
    S(i).yd = ym*rand(1,1);
    S(i).E = Eo*(1+u);
    S(i).type = 'U';
    S(i).G = 0;
    S(i).CH = 0;
    S(i).D = 0;
    S(i).R = 0;
    S(i).min_dis_cluster=0;
    S(i).dis_to_cluster=0;
    figure(3);
    plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y');
    hold on;
end
 
%E_total= (n*m1*Eo*(1+u))+(n*m0*Eo*(1+b))+(n*m*Eo*(1+a))+(n*(1-m1-m0-m)*Eo);
E_total= 0;
for i=1:1:n
    E_total=E_total+S(i).E;
end
 
% Selection of Rendezvous Nodes
rn=1;
    for i=1:1:n
        if (S(i).type=='S')
        Rx= rand;
        if ((ym*(1-Rx)/2)<=S(i).yd && (ym*(1+Rx)/2)>=S(i).yd)
            if S(i).E > 0
            S(i).R= 1;
            R(rn).xd= S(i).xd;
            R(rn).yd= S(i).yd;
            R(rn).id = i;
            rn=rn+1;
            end
        end
        end
    end
 
for r=1:1:rmax
    r
    RES(r)=0;
    for i=1:1:n
        if S(i).E>0
        RES(r)= RES(r)+ S(i).E;
        end
    end
    dead=0; 
   for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            S(i).D=1;
            dead=dead+1;
        end
    end
    if (dead == n)
        break;
    end
    
    DEAD(r+1)=dead; 
    
    for i=1:1:n
        S(i).CH=0;
    end
    
    d_toBS = 0.765*(xm/2);
    k_opt = (sqrt(n/(2*pi)))*(sqrt(E_fs/E_mp))*(xm/(d_toBS^2));
    d_toCH = xm/(sqrt(2*pi*k_opt));
    E_round = L*((2*n*E_elec) + (n*E_DA) + (k_opt*E_mp*(d_toBS^4)) + (n*E_fs*(d_toCH^2)));
    R = E_total/E_round;
    
    N= round(n/(k_opt))+1;
    M= round(N/cr)+1;
    q= mod(r,M+2);
    
    % Average Energy of the network in this round
    Ei(r)= 0;
    for i=1:1:n
        if S(i).E>0
        Ei(r)=Ei(r)+S(i).E;
        end
    end
    E(r)= Ei(r)/n;
    
    % Calculation of election probability of nodes to become CH
   
    for i=1:1:n
        if S(i).E <= T_absolute
            p(i)= (c* p_opt *(1+u)*S(i).E)/((1+m*(a+m0*(-a+b+m1*(-b+u))))*E(r));
            %S(i).G= 1;
        else
            if S(i).type == 'N'
                p(i)= (p_opt *S(i).E)/((1+m*(a+m0*(-a+b+m1*(-b+u))))*E(r));
                %S(i).G= 1;
            end
            if S(i).type == 'A'
                p(i)= (p_opt*(1+a)*S(i).E)/((1+m*(a+m0*(-a+b+m1*(-b+u))))*E(r));
                %S(i).G= 1;
            end
            if S(i).type == 'S'
                p(i)= (p_opt*(1+b)*S(i).E)/((1+m*(a+m0*(-a+b+m1*(-b+u))))*E(r));
                %S(i).G= 1;
            end
            if S(i).type == 'U'
                p(i)= (p_opt*(1+u)*S(i).E)/((1+m*(a+m0*(-a+b+m1*(-b+u))))*E(r));
                %S(i).G= 1;
            end
        end
    end
    
    % Operation for epoch
    for i=1:1:n
        if( (mod(r, round(1/p(i))))==0)                             % checking if node has been CH in last 1/p rounds
            S(i).G=0;
            %S(i).cl=0;
        end
    end
    
    for i=1:1:n
        if S(i).G ==0                                                % i.e. if the node belongs to G
            S(i).T = p(i)/(1-(p(i)*(mod(r,round(1/p(i))))));             % Threshold Probability
        else
            S(i).T = 0;
        end
    end
    
    % Election of CLUSTER HEADS (CHs)
    if q==1
    cluster=1;
    for xcord=0:xm/6:((5*xm)/6)
        for ycord=0:ym/4:((3*ym)/4)
            for i=1:1:n
                if (S(i).xd>= xcord && S(i).xd<(xcord+(xm/6)) && S(i).yd>= ycord && S(i).yd<(ycord+(ym/4)))
                    Rx = rand;
                    %disp('a')
                    if S(i).T >= 0.2
                        if S(i).R~=1
                            if S(i).E>0
                                %disp('b')
                                S(i).CH = 1;                               % node becomes CH
                                C(cluster).id= cluster;
                                C(cluster).xd= S(i).xd;
                                C(cluster).yd= S(i).yd;
                                C(cluster).E= S(i).E;
                                S(i).G = round(1/p(i))-1;
                                break;
                            end
                        end
                    end
                end
            end
            cluster=cluster+1;
        end
    end
    
    cluster1=cluster-1;
    
   for xcord=0:xm/6:((5*xm)/6)
        for ycord=0:ym/4:((3*ym)/4)
            countCH=0;
            for i=1:1:n
                if (S(i).xd>= xcord && S(i).xd<(xcord+(xm/6)) && S(i).yd>= ycord && S(i).yd<(ycord+(ym/4)))
                    if S(i).CH==1
                        countCH=countCH+1;
                        break;
                    end
                end
            end
            if countCH==0
                maxenergy=0;
                for i=1:1:n
                    if (S(i).xd>= xcord && S(i).xd<(xcord+(xm/6)) && S(i).yd>= ycord && S(i).yd<(ycord+(ym/4)))
                        if S(i).E > maxenergy
                            maxenergy=S(i).E;
                        end
                    end
                end
                for i=1:1:n
                    if (S(i).xd>= xcord && S(i).xd<(xcord+(xm/6)) && S(i).yd>= ycord && S(i).yd<(ycord+(ym/4)))
                        if S(i).E == maxenergy
                            for x=1:1:cluster1
                                if C(x).id==0
                                    C(x).id=x;
                                    C(x).xd=S(i).xd;
                                    C(x).yd=S(i).yd;
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        end
   end
 
    %determining number of nodes in each cluster
    member=zeros(cluster1,1);
    for xcord=0:xm/6:((5*xm)/6)
        for ycord=0:ym/4:((3*ym)/4)
            % INSIDE A CLUSTER
            for j=1:1:cluster1
                if (C(j).xd>= xcord & C(j).xd<=(xcord+(xm/6)) & C(j).yd>= ycord & C(j).yd<=(ycord+(ym/4)))
                    for i=1:1:n
                        if (S(i).xd>= xcord && S(i).xd<=(xcord+(xm/6)) && S(i).yd>= ycord && S(i).yd<=(ycord+(ym/4)))
                            if S(i).R~=1
                                if S(i).E>0
                                   member(j)=member(j)+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    end
    
               
gamma = 0.5;                  % weighted variable
w=zeros(n,1);
threshold = sqrt(((xm/18)^2)+((ym/12)^2));
% construction of sub-clusters in each cluster
for xcord=0:xm/6:((5*xm)/6)
    for ycord=0:ym/4:((3*ym)/4)
        % INSIDE A CLUSTER
        % Determine the CH in this cluster
        if q==1
        E_residual=0;
        for j=1:1:cluster1
            if (C(j).xd>= xcord & C(j).xd<=(xcord+(xm/6)) & C(j).yd>= ycord & C(j).yd<=(ycord+(ym/4)))
                break;
            end
        end
        E_residual=C(j).E;                  % residual energy of cluster head
        end
        s=1;
        t=1;
        % GOING INTO SUB-CLUSTERS
        for xcoor=xcord:xm/18:(xcord+(xm/6)-(xm/18))
            for ycoor=ycord:ym/12:(ycord+(ym/4)-(ym/12))
                % INSIDE A SUB-CLUSTER
                % count the number of nodes in sub-cluster
                
                nodes=zeros(n,1);               % nodes transferring data to CH through relay nodes
                for i=1:1:n
                    if (S(i).CH~=1 && S(i).R~=1)
                        if (S(i).xd>= xcoor && S(i).xd<=(xcoor+(xm/18)) && S(i).yd>= ycoor && S(i).yd<=(ycoor+(ym/12)))
                                nodes((10*s)+t)=nodes((10*s)+t)+1; 
                        end
                    end
                end
                % Calculation of w for each node in this sub-cluster
                if q==1
                min_w = INFINITY;
                for i=1:1:n
                    if (S(i).R~=1)
                    if (S(i).xd>= xcoor && S(i).xd<=(xcoor+(xm/18)) && S(i).yd>= ycoor && S(i).yd<=(ycoor+(ym/12)))
                        % Calculation of summation term
                        term=0;
                        for k=1:1:n
                            if (S(k).xd>= xcoor && S(k).xd<=(xcoor+(xm/18)) && S(k).yd>= ycoor && S(k).yd<=(ycoor+(ym/12)))
                                term=term+(((sqrt(((S(i).xd - S(k).xd)^2)+((S(i).yd - S(k).yd)^2)))-(sqrt(((S(k).xd - C(j).xd)^2)+((S(k).yd - C(j).yd)^2))))^2);
                            end
                        end
                        w(i) = (gamma/C(j).E) + ((1-gamma)/nodes((10*s)+t))*term;
                        if w(i) < min_w
                            min_w = w(i);
                        end
                    end
                    end
                end
                end
              
                % determining the relay node of the sub-cluster
                if q==2
                dis_to_relaynode=0;
                for i=1:1:n
                    if (S(i).xd>= xcoor && S(i).xd<=(xcoor+(xm/18)) && S(i).yd>= ycoor && S(i).yd<=(ycoor+(ym/12)))
                        if w(i)==min_w
                            S(i).relaynode=1;
                            % energy dissipated by relay nodes in inter-sub-cluster transmissions
                            distoCH= sqrt(((S(i).xd - C(j).xd)^2)+((S(i).yd - C(j).yd)^2));
                            if S(i).E>0
                            S(i).E = S(i).E - (nodes((10*s)+t))*L*(E_elec+E_DA);          % relay node receives and aggregates data
                            end
                            if S(i).E>0
                            S(i).E = S(i).E - ((nodes((10*s)+t)+1))*( (E_elec*(L)) + (E_fs*L*((distoCH)^2)));  % relay node transmits data to CH
                            end
                            % calculation of energy dissipated by nodes in intra-sub-cluster txns.
                            for k=1:1:n
                                if S(k).R~=1
                                if (S(k).xd>= xcoor && S(k).xd<=(xcoor+(xm/18)) && S(k).yd>= ycoor && S(k).yd<=(ycoor+(ym/12)))
                                    dis_to_relaynode = sqrt( ((S(i).xd - S(k).xd)^2) + ((S(i).yd - S(k).yd)^2) );
                                    if S(k).E>0
                                    S(k).E = S(k).E - ( E_elec*(L) + E_fs*L*((dis_to_relaynode)^2)); 
                                    end
                                end
                                end
                            end
                        end
                    end
                end
                end
              
                s=s+1;
                t=t+1;
            end
        end
    end
end
    
                
if q~=1
for o=1:1:cluster1                                % energy dissipated by CHs during reception of data
    if C(o).E >0
        C(o).E = C(o).E - (E_elec*L*member(o));
        for k=1:1:n
            if (S(k).xd==C(o).xd && S(k).yd==C(o).yd)
                S(k).E=C(o).E;
            end
        end
    end
end
end
 
%Energy dissipated in inter-cluster transmissions
if q~=1
dis_to_rn=0;
for l=1:1:cluster1
    min_dis_to_rn=INFINITY;
    % determination of RN at least distance
    for i=1:1:n
        if S(i).E>0
            if S(i).R==1
                dis_to_rn= sqrt(((S(i).xd-C(l).xd)^2)+((S(i).yd-C(l).yd)^2));
                if dis_to_rn < min_dis_to_rn
                    min_dis_to_rn= dis_to_rn;
                end
            end
        end
    end
    %Energy dissipated by RNs for reception
    for i=1:1:n
        if S(i).E>0
            if S(i).R==1
                if (sqrt(((S(i).xd-C(l).xd)^2)+((S(i).yd-C(l).yd)^2)) == min_dis_to_rn)
                    S(i).E = S(i).E - (E_elec+E_DA)*L;
                end
            end
        end
    end
    % energy dissipated by CHs for transmission of data to RNs
    if min_dis_to_rn > do
        if C(l).E>0
            C(l).E = C(l).E - ( E_elec*(L) + E_mp*L*((min_dis_to_rn)^4));
            for k=1:1:n
                if (S(k).xd==C(l).xd && S(k).yd==C(l).yd)
                    S(k).E=C(l).E;
                end
            end
        end
    else
        if C(l).E>0
            C(l).E = C(l).E - ( E_elec*(L) + E_fs*L*((min_dis_to_rn)^2));
            for k=1:1:n
                if (S(k).xd==C(l).xd && S(k).yd==C(l).yd)
                    S(k).E=C(l).E;
                end
            end
        end
    end
end
end
end
 
x=1:1:r;
y=1:1:r;
w=1:1:r;
 
for i=1:r;
    x(i)=i;
    y(i) = n- DEAD(i);   % remaining number of live nodes
    w(i)=RES(i);
end
figure(1)
plot(x,y,'r');
hold on;
figure(2)
plot(x,w,'r');
hold on;
