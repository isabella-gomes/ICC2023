clear all
%close all
clc

potNoise = 10^(-100/10)/1000; %noise power
Mt = 40; %number of transmit antennas at the BS
N = 4; %number of receive antennas at the users
K = 2; %number of users
SNRdb = 0:3:21;
gammath = 10.^(SNRdb/10);
sigma2 = 10^(-100/10)/1000;

%% Creating the Grid
dpos1 = 0:100:1000;
dpos2 = 0:100:1000;
[A,B] = meshgrid(dpos1,dpos2);
c=cat(2,A',B');
coord=sortrows(reshape(c,[],2)); %all coordinates from grid (100x100)
NX = length(dpos1);
NY = length(dpos2);
nel = (NX-1)*(NY-1) ;       % Number of elements in the surface
nnel = 5 ;                  % Number of nodes per element
nodes = zeros(nnel,nel) ;   % Initialize nodal connectivity matrix
count = 0 ;
x = A;
y = B;
for i = 1:NX-1
    for j = 1:NY-1
        l = (j-1)*(NX)+i ;
        count = count+1 ;
        nodes(:,count) =[l l+1 l+NX+1 l+NX l];  
    end
end
 
xm = mean(x(nodes(1:4,:))) ; %mid of each cell 
ym = mean(y(nodes(1:4,:))) ; 
[theta,rad] = cart2pol(xm,ym); 

potnoiseest = -10; %noise estimation power in dBm
potNoiseEst = 10.^(potnoiseest/10)/1000;
BS = [0,0]; %coordinates base station
BStheta = 0;
grid = [theta' rad'];
sterAdv = exp(-1i*2*pi*(1/2)*(0:1:(Mt-1)).*sin(grid(:,1))); %the estimated grid of positions for the target by the adversary considers the number of antennas that the base station has
temp = 1500; %number of observations
opc = zeros(length(SNRdb),4); %possible options (Detection, false alarm, missdetection, and undetection)
userspos = [800 100; 750 300]; %coordinates of users 
targetpos = [550 400]; %coordinates of the target
[targettheta, ~] = cart2pol(targetpos(1,1),targetpos(1,2));
[~, usersrad] = cart2pol(userspos(:,1),userspos(:,2)); 
phi = 3; %path loss exponent
H = zeros(N, Mt, K);
for i = 1:K
H(:,:,i) = (sqrt(usersrad(i).^(-phi)/2)).*(randn(N, Mt)+1i*randn(N, Mt)); %Rayleigh channel coefficient
end
%% Monte Carlo
for snr=1:length(SNRdb)
    [u, Wc,Wr,asteer,R, Rk, thetal,alfa,dthetal] = SDRMIMO(Mt,N, sigma2,K,H,gammath(snr),targettheta);
for MC=1:100
adv =1;
M = 500; %number of particles
T = temp;
q = zeros(M,1); %weight of the particles
Xinicart = rand(1,M).*1000; %initial  coordinate of the particles
Yinicart = 1000.*rand(1,M);
[Xtheta, Xrad] = cart2pol(Xinicart,Yinicart); 
Xini = [Xtheta',Xrad'];
qaux= repmat(1/M, M,1);
Xest = Xini;
aloc = ones(M,length(sterAdv));

for t=1:T
    %Wr = chol(Rr, 'lower');
    Wcaux = ones(Mt,K,N).*Wc;
    recWcaux = Wcaux + normrnd(0,sqrt(potNoiseEst),Mt,1,N);
    recWc = mean(recWcaux,3);
    Wraux =  Wr*ones(Mt,1);
    recWraux =  ones(Mt,1,N).*Wraux + normrnd(0,sqrt(potNoiseEst),Mt,1,N);
    recWr = mean(recWraux,3); 
    wEst = [recWc recWr];
    %wEst = [recWr];
    Ry = wEst*wEst';
    for i=1:length(thetal)
         py(i) = real(asteer(:,i)'*Ry*asteer(:,i)); %estimated beampattern
         %pys(i) = asteer(:,i)'*Rys*asteer(:,i);
         pR(i) = real(asteer(:,i)'*R*asteer(:,i)); %original beampattern
    end
     %
     pyamp = 5*py;
        %thed = linspace(0,90,901);
          %figure
          %plot(thed,paux,'r-')
%          hold on
%          plot(thed,py,'b-')
    dthetalaux = ones(1,length(thetal));
    [D,I] = pdist2(thetal',grid(:,1),'euclidean','Smallest',1);
    [Xx,Xy] = pol2cart(Xest(:,1),Xest(:,2));
    for j = 1:length(sterAdv)
        aloc(:,j) = inpolygon(Xx,Xy,x(nodes(:,j)),y(nodes(:,j))); %verify the particles within each cell
    end      
    for i=1:M   
        q(i) = (1-exp(-(pyamp(I(aloc(i,:)==1))))).*qaux(i); %recalculating the weights
    end
    qaux = q./sum(q); 
    for i = 1: length(sterAdv)
        qmedia(i) = sum(qaux(aloc(:,i)==1)); %taking the mean of the new weights calculated
    end
        for i = 1:M
         indices = multinomial(qmedia); %employing the multinomial resample method
         Xnewc(i) = min(x(nodes(:,indices)))+(max(x(nodes(:,indices)))-min(x(nodes(:,indices))))*rand(1);
         Ynewc(i) = min(y(nodes(:,indices)))+(max(y(nodes(:,indices)))-min(y(nodes(:,indices))))*rand(1);
         [Xnewt(i), Xnewr(i)] = cart2pol(Xnewc(i),Ynewc(i)); 
        end
         Xnew = [Xnewt' Xnewr']; %new coordinates of the particles
         qnew = repmat(1/M, M,1); %new weights to each particle
    Xest = Xnew;
    qaux = qnew;

    [Xcart,Ycart] = pol2cart(Xest(:,1),Xest(:,2));
    Xesttime(:,t) = Xcart;
    Yesttime(:,t) = Ycart;
end
 [Xf,Yf] = pol2cart(Xest(:,1),Xest(:,2)); %final position of the particles 
 for j = 1:length(sterAdv)
     alocf(:,j) = inpolygon(Xf,Yf,x(nodes(:,j)),y(nodes(:,j))); %verify the particles within each cell
     alocrad(:,j) = inpolygon(targetpos(1,1),targetpos(1,2),x(nodes(:,j)),y(nodes(:,j))); %check in which cell is the radar
 end 
 [v,ind]= max(sum(alocf)/M);
angAdv(MC) = rad2deg(grid(alocf(ind,:)'));
%angRad(MC) = rad2deg(grid(alocrad(1,:)'));
angRad(MC) = rad2deg(targettheta);
if abs(angAdv(MC)-angRad(MC))<=10  && max(sum(alocf)/M)>=0.9
    opc(snr,1) = opc(snr,1)+1;
elseif abs(angAdv(MC)-angRad(MC))>10 && max(sum(alocf)/M)>=0.9
    opc(snr,2) = opc(snr,2)+1;
elseif abs(angAdv(MC)-angRad(MC))<=10 && max(sum(alocf)/M)<0.9
    opc(snr,3) = opc(snr,3)+1;
elseif abs(angAdv(MC)-angRad(MC))>10 && max(sum(alocf)/M)<0.9
    opc(snr,4) = opc(snr,4)+1;
end
end
end

% animatedimage(Xesttime,Yesttime,targetpos, userspos) %gif to evaluate the
% convergence of particles
