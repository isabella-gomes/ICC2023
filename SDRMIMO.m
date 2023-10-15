function [u, wi,Wr,asteer,R, Rk, thetal,alfa,dthetal] = SDRMIMO(Mt,N, sigma2,K,H,gammath,targettheta)
Pt = 1;
thetal = 0:pi/1800:pi/2; %angle grid sample
Delta = pi/18;
a = (1i*2*pi*(1/2)*(1:1:(Mt-1)).*sin(thetal'))';
%atargets = (i*2*pi*(1/2)*(1:1:(M-1)).*sin([theta2 theta1 theta3]'))';
asteer = [ones(1,length(thetal));exp(a)];
%asteertargets = [ones(1,3);exp(atargets)];
theta1 = targettheta;
dthetal = zeros(1,length(thetal));
thetap = [(theta1-Delta/2) (theta1+Delta/2)];
dthetal((thetal>=thetap(1)&thetal<=thetap(2)))=1; %desired beampattern
u = rand(N,K);
epsilon = 0.01;
normu = 10;
%u = 0;
beam = zeros(1,length(thetal));
iterac = 0;
while normu> epsilon 
    uaux = beam;
          cvx_clear;
            cvx_begin quiet
            variables R(Mt,Mt)  Rk(Mt,Mt,K)  alfa(1,1)
            expression beam(1,length(thetal)); 
                for l=1:length(thetal)
                    beam(l) = pow_abs(alfa.*dthetal(:,l)- real(asteer(:,l)'*R*asteer(:,l)),2);  
                end
            minimize ((1/length(thetal))*sum(beam));
            subject to
                R == semidefinite(Mt);
                diag(R)== Pt/Mt;
                for i = 1:K
                    Rk(:,:,i)==semidefinite(Mt);
                    %Raux = Raux + Rk;
                     ((1 - (1/gammath)).*real(u(:,i)'*H(:,:,i)*Rk(:,:,i)*H(:,:,i)'*u(:,i)))>= (real(u(:,i)'*H(:,:,i)*R*H(:,:,i)'*u(:,i)-u(:,i)'*H(:,:,i)*Rk(:,:,i)*H(:,:,i)'*u(:,i))+sigma2.*norm(u(:,i),2));
                    %Raux = Raux + rk;
                end
               (R-sum(Rk,3))==semidefinite(Mt);
            cvx_end;

        wi = zeros(Mt,K);
        Rktil = zeros(Mt,Mt,K);
    for i = 1:K
        %v = u(:,i)'*H(:,:,i)
        wi(:,i) = real(((u(:,i)'*H(:,:,i)*Rk(:,:,i)*H(:,:,i)'*u(:,i)))^(-1/2).*(Rk(:,:,i)*H(:,:,i)'*u(:,i)));
        aux(i) = ((u(:,i)'*H(:,:,i)*Rk(:,:,i)*H(:,:,i)'*u(:,i)))^(-1/2);
        %Wk(1+(i-1)*Mt:i*Mt,:) = wi;
        %wk=permute(Wk(i,:,:),[2 3 1]); % received signal at each user (M x N)
        Rktil(:,:,i) =  wi(:,i)*wi(:,i)';
        %Raux1 = Raux1 + Rk;
    end
    Rr = R - sum(Rktil,3);
   % Wr = ((real(H')*WrWrh*real(H))^(-1/2)).*(WrWrh*real(H));
    Wr = chol(Rr,'lower');
    sumWc = zeros(N, N, K);
    sumWr = zeros(N, N, K);
    for i = 1:K 
        waux = wi;
        waux(:,i)=[];
        for j = 1:K-1
        sumWc(:,:,i) = sumWc(:,:,i) + real(H(:,:,i)*waux(:,j)*waux(:,j)'*H(:,:,i)');% Wc Transmit precoder
        end
    end
    for i = 1:K
        for j = 1:Mt
        sumWr(:,:,i) = sumWr(:,:,i) + real(H(:,:,i)*Wr(:,j)*Wr(:,j)'*H(:,:,i)'); % Wr transmit precoder
        end
    end
    for i = 1:K
         u(:,i) = real(inv(sumWc(:,:,i)+sumWr(:,:,i)+sigma2*eye(N))*H(:,:,i)*wi(:,i)); %receive beamformer
    end 
    %display(beam)
    normu = norm(beam-uaux)/norm(beam)
    iterac = iterac + 1;
end
 end