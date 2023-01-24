clear;

M1 = 200; % Number of realizations
N = 32; %Total Number of Chromosomes
M = 23600; %Total Number of Time Steps

dist_count = zeros(N/2,13,M1);
dist_count1 = zeros(N/2,13,M1);
dist_count_nh = zeros(N/2,13,M1);
dist_count1_nh = zeros(N/2,13,M1);

chrom_pair = zeros(13);
chrom_total = zeros(13);

figure(20)

for pp = 1:M1
    pp

% Parameters
Ca = .003; %Attractive Strength
Cr = .05; %Repuslive Strength (non-matching chromosome)
Cr1 = .05; %Repuslive Strength for matching chromosome enforcing body size
% la = 0.7; % Length scale for attractive force
% lr = 0.4; %Length scale for repuslive force (non-matching)
% lr1 = 0.05; %length scale for repuslive force (matching)

    la(1) = 0.7*(1.53/5.03);
    lr(1) = 0.4*(1.53/5.03);
    lr1(1) = 0.05*(1.53/5.03);
    la(17) = 0.7*(1.53/5.03);
    lr(17) = 0.4*(1.53/5.03);
    lr1(17) = 0.05*(1.53/5.03);
     la(2) = 0.7*(5.4/5.03);
    lr(2) = 0.4*(5.4/5.03);
    lr1(2) = 0.05*(5.4/5.03);
    la(18) = 0.7*(5.4/5.03);
    lr(18) = 0.4*(5.4/5.03);
    lr1(18) = 0.05*(5.4/5.03);
    la(3) = 0.7*(2.13/5.03);
    lr(3) = 0.4*(2.13/5.03);
    lr1(3) = 0.05*(2.13/5.03);
    la(19) = 0.7*(2.13/5.03);
    lr(19) = 0.4*(2.13/5.03);
    lr1(19) = 0.05*(2.13/5.03);
        la(4) = 0.7*(10.2/5.03);
    lr(4) = 0.4*(10.2/5.03);
    lr1(4) = 0.05*(10.2/5.03);
    la(20) = 0.7*(10.2/5.03);
    lr(20) = 0.4*(10.2/5.03);
    lr1(20) = 0.05*(10.2/5.03);
        la(5) = 0.7*(3.86/5.03);
    lr(5) = 0.4*(3.86/5.03);
    lr1(5) = 0.05*(3.86/5.03);
    la(21) = 0.7*(3.86/5.03);
    lr(21) = 0.4*(3.86/5.03);
    lr1(21) = 0.05*(3.86/5.03);
        la(6) = 0.7*(1.8/5.03);
    lr(6) = 0.4*(1.8/5.03);
    lr1(6) = 0.05*(1.8/5.03);
    la(22) = 0.7*(1.8/5.03);
    lr(22) = 0.4*(1.8/5.03);
    lr1(22) = 0.05*(1.8/5.03);
        la(7) = 0.7*(7.26/5.03);
    lr(7) = 0.4*(7.26/5.03);
    lr1(7) = 0.05*(7.26/5.03);
    la(23) = 0.7*(7.26/5.03);
    lr(23) = 0.4*(7.26/5.03);
    lr1(23) = 0.05*(7.26/5.03);
        la(8) = 0.7*(3.73/5.03);
    lr(8) = 0.4*(3.73/5.03);
    lr1(8) = 0.05*(3.73/5.03);
    la(24) = 0.7*(3.73/5.03);
    lr(24) = 0.4*(3.73/5.03);
    lr1(24) = 0.05*(3.73/5.03);
        la(9) = 0.7*(2.93/5.03);
    lr(9) = 0.4*(2.93/5.03);
    lr1(9) = 0.05*(2.93/5.03);
    la(25) = 0.7*(2.93/5.03);
    lr(25) = 0.4*(2.93/5.03);
    lr1(25) = 0.05*(2.93/5.03);
        la(10) = 0.7*(5.0/5.03);
    lr(10) = 0.4*(5.0/5.03);
    lr1(10) = 0.05*(5.0/5.03);
    la(26) = 0.7*(5.0/5.03);
    lr(26) = 0.4*(5.0/5.03);
    lr1(26) = 0.05*(5.0/5.03);
        la(11) = 0.7*(4.46/5.03);
    lr(11) = 0.4*(4.46/5.03);
    lr1(11) = 0.05*(4.46/5.03);
    la(27) = 0.7*(4.46/5.03);
    lr(27) = 0.4*(4.46/5.03);
    lr1(27) = 0.05*(4.46/5.03);
        la(12) = 0.7*(7.2/5.03);
    lr(12) = 0.4*(7.2/5.03);
    lr1(12) = 0.05*(7.2/5.03);
    la(28) = 0.7*(7.2/5.03);
    lr(28) = 0.4*(7.2/5.03);
    lr1(28) = 0.05*(7.2/5.03);
        la(13) = 0.7*(6.13/5.03);
    lr(13) = 0.4*(6.13/5.03);
    lr1(13) = 0.05*(6.13/5.03);
    la(29) = 0.7*(6.13/5.03);
    lr(29) = 0.4*(6.13/5.03);
    lr1(29) = 0.05*(6.13/5.03);
        la(14) = 0.7*(5.2/5.03);
    lr(14) = 0.4*(5.2/5.03);
    lr1(14) = 0.05*(5.2/5.03);
    la(30) = 0.7*(5.2/5.03);
    lr(30) = 0.4*(5.2/5.03);
    lr1(30) = 0.05*(5.2/5.03);
        la(15) = 0.7*(7.26/5.03);
    lr(15) = 0.4*(7.26/5.03);
    lr1(15) = 0.05*(7.26/5.03);
    la(31) = 0.7*(7.26/5.03);
    lr(31) = 0.4*(7.26/5.03);
    lr1(31) = 0.05*(7.26/5.03);
        la(16) = 0.7*(6.33/5.03);
    lr(16) = 0.4*(6.33/5.03);
    lr1(16) = 0.05*(6.33/5.03);
    la(32) = 0.7*(6.33/5.03);
    lr(32) = 0.4*(6.33/5.03);
    lr1(32) = 0.05*(6.33/5.03);


Crb = 0.05; %Boundary Repulsion Strength
lrb = 0.05; %Boundary Repulsive Distance
L = 3.25;  %Domain Radius
Diff = 0.0001; %Diffusion Constant
alpha = 0.3; %Movement speed
Time = 0.0;
count1 = -1;  %Count for determining when curves of distances are generated.
count1_nh = -1;  %Count for determining when curves of nh distances are generated.
countp = 0; %Count for determining intervals for plots
count_tj = 0;
count_tj1 = 0;
ctj1 = 0; %Zero before pairing, changes to 1 after pairing

count_sort = 0; %To do initial sorting, then keep that order
count_sort_nh = 0;

dt = .2; %time step
x1 = zeros(N,1);
y1 = zeros(N,1);


%Initial Conditions
%Zero out distance vector
for i =1:N/2
    for j = 1:6
        dist_count(i,j) = 0.0;
    end
end
%Random Initial Placement
for i=1:N
x(i) = rand()*(2*L)-L;
y(i) = rand()*(2*L)-L;
match(i) = 0;
test = 0;
while (test == 0) 
    pos = sqrt(x(i)^2+y(i)^2);
    if (pos < L)
        test = 1;
        vx(i) = rand()*2*1-1;
        vy(i) = rand()*2*1-1;
    elseif (pos >= L)
        x(i) = rand()*(2*L)-L;
        y(i) = rand()*(2*L)-L;
        test = 0;
    end
end
end

% Find Non-homologous Chromsomes with closest starting distance
for i = 1:N
     dist_chrom = 2*L;
       if (i <=N/2)
       dist_homolog = sqrt((x(i)-x(i+N/2))^2+(y(i)-y(i+N/2))^2);
       elseif (i > N/2)
       dist_homolog = sqrt((x(i)-x(i-N/2))^2+(y(i)-y(i-N/2))^2);
       end
   for j = 1:N
      if (i ~= j && abs(i-j) ~= N/2)
         dist_c1 = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
         dist_c2 = abs(dist_homolog-dist_c1);
         if (dist_c2 < dist_chrom)
             dist_chrom = dist_c2;
            chrom_nh_match(i) = j; 
         end
      end
   end
end


% Dynamics of Chromosomes
for i = 1:M
    Time = Time + 5.0*dt;
    %Loop over each chromosome
    for j = 1:N
        ch_repx(j) = 0.0;
        ch_repy(j) = 0.0;
        ch_intx(j) = 0.0;
        ch_inty(j) = 0.0;
        
        x1(j) = x(j);
        y1(j) = y(j);
        v1x(j) = alpha*(vx(j)/sqrt(vx(j)^2+vy(j)^2));
        v1y(j) = alpha*(vy(j)/sqrt(vx(j)^2+vy(j)^2));
        temp1 = v1x(j);
        temp2 = v1y(j);
        % Atttractive Interaction Between Pairs
        for k = 1:N
            if (j-k == N/2 || k-j == N/2)
            dist = sqrt((x1(j)-x1(k))^2+(y1(j)-y1(k))^2);
            ch_intx(j) = -(Ca/la(k))*(x1(j)-x1(k))*exp(-dist/la(k))/dist+(Cr1/lr1(k))*(x1(j)-x1(k))*exp(-dist/lr1(k))/dist;
            ch_inty(j) = -(Ca/la(k))*(y1(j)-y1(k))*exp(-dist/la(k))/dist+(Cr1/lr1(k))*(y1(j)-y1(k))*exp(-dist/lr1(k))/dist;
                if (dist < lr1(k) && j > N/2)
                    match(j) = 1;
                end
            elseif (j ~= k)
            dist = sqrt((x1(j)-x1(k))^2+(y1(j)-y1(k))^2);
            ch_repx(j) = ch_repx(j)+(Cr/lr(k))*(x1(j)-x1(k))*exp(-dist/lr(k))/dist;
            ch_repy(j) = ch_repy(j)+(Cr/lr(k))*(y1(j)-y1(k))*exp(-dist/lr(k))/dist;
            end
        end
%           % Boundary Repulsion       
             theta = atan2(y(j),x(j));
             distb = sqrt((x(j)-L*cos(theta))^2+(y(j)-L*sin(theta))^2);
             ch_repx(j) = ch_repx(j)+(Crb/lrb)*(x1(j)-L*cos(theta))*exp(-distb/lrb)/distb;
             ch_repy(j) = ch_repy(j)+(Crb/lrb)*(y1(j)-L*sin(theta))*exp(-distb/lrb)/distb;
             
        if (match(j) ~= 1)
            diffx(j) = Diff*normrnd(0,sqrt(dt));
            diffy(j) = Diff*normrnd(0,sqrt(dt));
        x(j) = x1(j) + dt*(temp1+ch_intx(j)+ch_repx(j))+diffx(j);
        y(j) = y1(j) + dt*(temp2+ch_inty(j)+ch_repy(j))+diffy(j);
        vx(j) = temp1+ch_intx(j)+ch_repx(j)+diffx(j);
        vy(j) = temp2+ch_inty(j)+ch_repy(j)+diffy(j);
        
        elseif (match(j) == 1 && j > N/2)
            x(j) = x1(j) + dt*(v1x(j-N/2)+ch_intx(j-N/2)+ch_repx(j-N/2))+diffx(j-N/2);
            y(j) = y1(j) + dt*(v1y(j-N/2)+ch_inty(j-N/2)+ch_repy(j-N/2))+diffy(j-N/2);
            vx(j) =v1x(j-N/2)+ch_intx(j-N/2)+ch_repx(j-N/2);
            vy(j) = v1y(j-N/2)+ch_inty(j-N/2)+ch_repy(j-N/2);
        end
    end
    
    %Boundary Conditions
    for j = 1:N
       dist1 = sqrt(x(j)*x(j)+y(j)*y(j));

       if (dist1 > L) %Reorienting at the boundary
%           vx(j) = -vx(j);
%           vy(j) = -vy(j);
          x(j) = .99*L*x(j)/sqrt(x(j)*x(j)+y(j)*y(j));
           y(j) = .99*L*y(j)/sqrt(x(j)*x(j)+y(j)*y(j));
       end 
    end
    
    if (pp == 1)
        if (match(N/2+1) == 0)
    count_tj = count_tj+1;
    
    xtj(count_tj) = x(1);
    ytj(count_tj) = y(1);
    x2tj(count_tj) = x(N/2+1);
    y2tj(count_tj) = y(N/2+1);
    x3tj(count_tj) = x(chrom_nh_match(1));
    y3tj(count_tj) = y(chrom_nh_match(1));
    
    elseif (match(N/2+1) == 1)
        ctj1 = 1; %parameter flips to 1 if pairing occurs
        count_tj1 = count_tj1+1;
    xtj1(count_tj1) = x(1);
    ytj1(count_tj1) = y(1);
    x2tj1(count_tj1) = x(N/2+1);
    y2tj1(count_tj1) = y(N/2+1);
    x3tj1(count_tj1) = x(chrom_nh_match(1));
    y3tj1(count_tj1) = y(chrom_nh_match(1));
        
        end
    end

       countp = countp + 1;
       if (countp > 59 || count1 == -1) %5 origially
           countp = 0;
       
    

    x11 = 1:N/2;
    
    count1 = count1+1;
    for p = 1:13
        if (count1 == (p-1)*30) 
            for ii = 1:N/2
                for jj = (N/2+1):N
                    if (jj-ii == N/2)
                        dist_count(ii,p,pp) = sqrt((x1(ii)-x1(jj))^2+(y1(ii)-y1(jj))^2);
                    end
                end
            end
            %Sort 1st time and then no more
            %if (count_sort == 0)
            %    count_sort = 1;
            %    [dist_count1(:,p,pp),I] = sort(dist_count(:,p,pp));
            %else
            dist_count1(:,p,pp) = dist_count(:,p,pp);
            %end
                        % Paired Percent per time
            for q = (N/2+1):N
                chrom_total(p) = chrom_total(p) + 1;
                if (match(q) == 1)
                    chrom_pair(p) = chrom_pair(p) + 1;
                end
            end
        end
    end

    x11 = 1:N/2;
    
    count1_nh = count1_nh+1;
    for p = 1:13
        if (count1_nh == (p-1)*30) 
            
            for ii = 1:N/2
                dist_count_nh(ii,p,pp) = 0.0; %2*L;
                % Add in non-homolog distances closest to initial distance
                % of homolog
                temp1 = sqrt((x1(ii)-x1(chrom_nh_match(ii)))^2+(y1(ii)-y1(chrom_nh_match(ii)))^2);
                temp2 = sqrt((x1(ii+N/2)-x1(chrom_nh_match(ii+N/2)))^2+(y1(ii+N/2)-y1(chrom_nh_match(ii+N/2)))^2);
                
                dist_count_nh(ii,p,pp) = .5*(temp1+temp2);
                
            end
            
             %Sort 1st time and then no more
            %if (count_sort_nh == 0)
            %    count_sort_nh = 1;
            %    [dist_count1_nh(:,p,pp),I1] = sort(dist_count_nh(:,p,pp));
            %else
            dist_count1_nh(:,p,pp) = dist_count_nh(:,p,pp);
            %end           
            
            %dist_count1_nh(:,p,pp) = sort(dist_count_nh(:,p,pp));
        end
    end


    
 %   hold off
 %   drawnow()
       end
   % pause
end

%plot(x(1),y(1),'ro',x(2),y(2),'bo')

    if (pp == 1)
    traj2(:,1) = xtj;
    traj2(:,2) = x2tj;
    traj2(:,3) = x3tj;
    traj2(:,4) = ytj;
    traj2(:,5) = y2tj;
    traj2(:,6) = y3tj;
    dlmwrite(['trajectory_before_pair_data.dat'],traj2,'delimiter','\t','precision',5)
        if (ctj1 == 1)
    traj1(:,1) = xtj1;
    traj1(:,2) = x2tj1;
    traj1(:,3) = x3tj1;
    traj1(:,4) = ytj1;
    traj1(:,5) = y2tj1;
    traj1(:,6) = y3tj1;
    dlmwrite(['trajectory_after_pair_data.dat'],traj1,'delimiter','\t','precision',5)
        end
    end

end


%Compute Average an Standard Deviation over data
 
      for ii = 1:N/2
          for p = 1:13
 homolog_avg_mc(ii,p) = mean(dist_count1(ii,p,:));
 homolog_sd_mc(ii,p) = std(dist_count1(ii,p,:))/sqrt(M1);
 nonhomolog_avg_mc(ii,p) = mean(dist_count1_nh(ii,p,:));
 nonhomolog_sd_mc(ii,p) = std(dist_count1_nh(ii,p,:))/sqrt(M1);
          end
      end
         

C1 = 1:N/2;

clf()
subplot(1,2,1)

hold on

    C = {'k','b','r','g','c','m',[0 .5 0],[.8 .2 .6],[.62 .66 .12],[.25 .25 0.25],[.2 .8 .6],[.8 .6 .2],[.6 .2 .8]};
for p = 1:13

errorbar(C1,homolog_avg_mc(:,p),homolog_sd_mc(:,p),'o-','MarkerSize',6,'CapSize',10,'color',C{p})
end
xlabel('Order','FontSize',20)
ylabel('Focus distance (\mu m)','FontSize',20)
title('Homologous Chromosome','FontSize',20)
legend('t = 3.0h','t = 3.5h','t = 4.0h','t = 4.5h','t = 5.0h','t = 5.5h','t = 6.0h','t = 6.5h','t = 7.0h','t = 7.5h','t = 8.0h','t = 8.5h','t = 9.0h','Location','NorthWest')
axis([0, N/2, 0, 8])

subplot(1,2,2)

hold on
for p = 1:13
errorbar(C1,nonhomolog_avg_mc(:,p),nonhomolog_sd_mc(:,p),'o-','MarkerSize',6,'CapSize',10,'color',C{p})
end
xlabel('Order','FontSize',20)
ylabel('Focus distance (\mu m)','FontSize',20)
title('Nonhomologous Chromosome','FontSize',20)
legend('t = 3.0h','t = 3.5h','t = 4.0h','t = 4.5h','t = 5.0h','t = 5.5h','t = 6.0h','t = 6.5h','t = 7.0h','t = 7.5h','t = 8.0h','t = 8.5h','t = 9.0h','Location','NorthWest')
axis([0 N/2 0 8])

drawnow()

for p = 1:13 
     err1(:,1) = C1;
     err1(:,2) = homolog_avg_mc(:,p);
     err1(:,3) = homolog_sd_mc(:,p);
 
 dlmwrite(['homolog_data_',num2str(p),'.dat'],err1,'delimiter','\t','precision',5)
end
for p = 1:13 
     err1(:,1) = C1;
     err1(:,2) = nonhomolog_avg_mc(:,p);
     err1(:,3) = nonhomolog_sd_mc(:,p);
 
 dlmwrite(['nonhomolog_data_',num2str(p),'.dat'],err1,'delimiter','\t','precision',5)
end

tt1 = 3:.5:9;
cpair(:,1) = tt1;
for q = 1:13
cpair(q,2) = chrom_pair(q)/chrom_total(q);
end
dlmwrite(['chrom_pair_v_time.dat'],cpair,'delimiter','\t','precision',5)


function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r*cos(th) + x;
yunit = r*sin(th) + y;
h = plot(xunit,yunit,'k-','LineWidth',3);
hold off
end
