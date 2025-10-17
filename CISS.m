period=100;
I=1120;  
N=16;   

d=0.049;
norm_h=zeros(I,period);
squared_error_KF = zeros(I, period);

px = zeros(I, period);
py = zeros(I, period);
pz = zeros(I, period);

px_est_KF = zeros(I, period);
py_est_KF = zeros(I,period);
pz_est_KF = zeros(I, period);
px_est_NCNKF = zeros(I, period);
py_est_NCNKF = zeros(I,period);
pz_est_NCNKF = zeros(I, period);
r = zeros(I, period);
abs_val = zeros(I, period);
squared_error_NCNKF= zeros(I, period);
normalized_error_KF= zeros(I, period);
normalized_error_NCNKF= zeros(I, period);
MSE1=zeros(I,period);
MSE2=zeros(I,period);
MSE_h3=zeros(I,period);
MSE_h4=zeros(I,period);

SNR = zeros(I, period);

squared_error_KF_music = zeros(I, period);
N_side = sqrt(N);
Fish1=zeros(period,1);
Fish4=zeros(period,1);
for u=1:period  
  
%% system setup& data generation
%N=16;           %number of elements
M=20;            %number of antanas
Pt=10;        %transmit power
fade_var = 0.1; % variance of  Rayleigh channel


weight_P=0.1;
sigma_W=sqrt(0.00001);  


R=1*sigma_W^2;
q=0.01;
Q=q*eye(N);

speed_light=300000000;
n=zeros(N,1);
phi = sqrtm(Pt) * dftmtx(N)/sqrt(N);      %get the optimal phi based on the lemma 1
Y=zeros(M,I);
N_side = sqrt(N);
X_positions = zeros(N,2);
speed_max=300/3.6;        %speed
f=300000000;     %carrier frequency
lambda=speed_light/f;
time_slot=1e-4;
h1=zeros(N,1);
varphi=zeros(I,1); theta1=zeros(I,1);
speed=zeros(I,1); epsilon=zeros(I,1);
aalpha1=zeros(I,1); bbeta1=zeros(I,1);
aalpha3=zeros(I,1); bbeta3=zeros(I,1);
row1=zeros(I,1); col1=zeros(I,1);
row3=zeros(I,1); col3=zeros(I,1);
C   = complex(zeros(M,      N,     I));   % C(:,:,i) = G · diag(φ)
C1  = zeros(2*M,  2*N,   I);              % C1(:,:,i) = [Re  −Im; Im  Re]
Y   = complex(zeros(M,      I));         
Y3  = zeros(2*M,     I);                 
%W2  = complex(zeros(M,      I));         
                  


h        = complex(zeros(N,     I));     
xn        = zeros(2*N,    I);             % [Re(h); Im(h)]


h_est_label1 = complex(zeros(N, I));      
h_est_label4 = complex(zeros(N, I));      
xn_est_label = zeros(2*N,   I);          



for i=1:N
    for jj=1:M
        G(jj,i) = sqrt(fade_var/2) * (normrnd(0,1,1,1) + 1i * normrnd(0,1,1,1)); 
    end
end

r(1,u)=50+100*rand(1);
varphi(1)=pi-2*pi*rand(1);
%varphi(1)=deg2rad(135);
theta1(1)=rand(1)*pi/2;
%theta1(1)=deg2rad(45);
 % px(1,u)=r(1,u)*cos(theta1(1))*cos(varphi(1));
 % py(1,u)=r(1,u)*cos(theta1(1))*sin(varphi(1));
 % pz(1,u)=r(1,u)*sin(theta1(1));
px(1,u)=20;
py(1,u)=20;
pz(1,u)=20;
ii=1;
for n=1:N_side
    for m=1:N_side
    xx=n*d-d*(1+sqrt(N))/2;
    yy=m*d-d*(1+sqrt(N))/2;
    h1(ii)=(1/r(1,u))*exp(-1j*2*pi*r(1,u)/lambda+1j*2*pi*(xx*px(1,u)+yy*py(1,u))/(lambda*r(1,u)));
    X_positions(ii,:) = [xx, yy];
    ii=ii+1;
    end
end


norm_h(1,u)=norm(h1);
g1=log(h1);  %Generate random channel h and corrsponding g for the first sub-time slot
xn1=[real(h1);imag(h1)];
h(:,1)=h1;
h_true(:, 1, u) = h(:,1);
g(:,1)=g1;
xn(:,1)=xn1;
speed(1)=speed_max * rand();
%speed(1)=80;
for i=2:I     
    speed(i)=speed_max * rand();
    %speed(i)=90;
   varphi(i)=-pi+2*pi*rand();
   % varphi(i)=varphi(i-1);
    theta1(i)= -pi / 2+ rand(1) * pi;
    %theta1(i)=theta1(i-1);

    epsilon(i)=speed(i)*time_slot;
    px(i,u)=px(i-1,u)+epsilon(i)*cos(theta1(i))*cos(varphi(i));
    py(i,u)=py(i-1,u)+epsilon(i)*cos(theta1(i))*sin(varphi(i));
    pz(i,u)=pz(i-1,u)+epsilon(i)*sin(theta1(i));
    r(i,u)=sqrt(px(i,u)^2+py(i,u)^2+pz(i,u)^2);
      
    ii=1;
for n=1:N_side
    for m=1:N_side
    %h1(i)=(0.01+0.01*(rand(1)-0.5))*exp(1j*2*pi*(rand(1)-0.5));
    xx=n*d-d*(1+sqrt(N))/2;
    yy=m*d-d*(1+sqrt(N))/2;
    h1(ii)=(1/r(i,u))*exp(-1j*2*pi*r(i,u)/lambda+1j*2*pi*(xx*px(i,u)+yy*py(i,u))/(lambda*r(i,u)));
    
    ii=ii+1;
    end
end
    g1=log(h1);
    xn1=[real(h1);imag(h1)];
    h(:,i)=h1;
    h_true(:, i, u) = h(:,i);
    norm_h(i,u)=norm(h1);
    g(:,i)=g1;
    xn(:,i)=xn1;
     W=(1/sqrt(2))*normrnd(0,sigma_W,M,1)+ (1/sqrt(2))*1i * normrnd(0,sigma_W,M,1); % Observation noise
    C(:,:,i)=G*diag(phi(mod(i,N)+1,:));
   
     C1(:,:,i)=[real(C(:,:,i)),-imag(C(:,:,i));imag(C(:,:,i)),real(C(:,:,i))];

     Y1=C(:,:,i)*h(:,i)+W;
     SNR(i,u)=(norm(C(:,:,i)*h(:,i)))^2/(M*sigma_W^2);
     Y(:,i)=Y1;
     
    D(:,:,i)=diag(h(:,i))*G';
    W2=(1/sqrt(2))*normrnd(0,sigma_W,M,N)+ (1/sqrt(2))*1i * normrnd(0,sigma_W,M,N);
    Y2(:,:,i)=phi.'*D(:,:,i)+W2.';
     Y3(:,i)=[real(Y1);imag(Y1)];
     

end




%% modified KF    

h_est1=zeros(N,1);      %channel estimation
P_pri1=zeros(N,N);      %prior covariance
P_post1=zeros(N,N);      %posterior covariance
K1=zeros(N,M);           %Kalman filter
%kalman filter algorithm
P_post1=weight_P*eye(N);
average_value = mean(h(:,1), 'all')/2;
%h_est1 = ones(N,1)*average_value;
h_est1 =0.5*h(:,1);
%h_est1 = zeros(N,1);


MSE1(1,u)=(norm(h(:,1)-h_est1))^2;
h_est_label1(:,1)=h_est1;
for i=2:I

       
        Q11=q*eye(N)*diag((h_est1.*conj(h_est1)));
            
        P_pri1=P_post1+Q11;
        K1=P_pri1*(C(:,:,i))'/((C(:,:,i))*P_pri1*(C(:,:,i))'+R*eye(M));

        h_est1=h_est1+K1*(Y(:,i)-(C(:,:,i))*h_est1);
        h_est_label1(:,i)=h_est1;
        P_post1=(eye(N)-K1*(C(:,:,i)))*P_pri1;
       
        
        MSE1(i,u)=(norm(h(:,i)-h_est1,2))^2;




end


%% Other paper
rho=0.97;
P_min=1e-12;
f_dopshift=speed_max*f/300000000; %maximum doppler shift frequency
Th=20e-6; %coherence time of h
alpha_ref=besselj(0,2*pi*f_dopshift*time_slot); %tried to use zero-order

h_est2=zeros(N,1); %channel estimation
P_pri2=zeros(N,N); %prior covariance
P_post2=zeros(N,N); %posterior covariance
K2=zeros(N,N); %Kalman filter
P_post2=weight_P*eye(16); 
h_est2 = 2*h(:,1);
D_est=ones(N,M)*0.1;
MSE2(1,u)=(norm(h(:,1)-h_est2))^2;
for o=1:16
    MSE2(o,u)=(norm(h(:,1)-h_est2))^2; %Initial value of MSE_h 
end
II=I/16;
for i=2:II
    D_est=alpha_ref*D_est;
    P_pri2=alpha_ref^2*P_post2+M*(1-alpha_ref^2)*eye(N);
    K2=P_pri2*conj(phi) *inv((phi.'*P_pri2*conj(phi)+sigma_W^2*eye(N)));
    D_est=D_est+K2*(Y2(:,:,i)-phi.'*D_est);
    P_post2=(eye(N)-K2*phi')*P_pri2; 
    P_post2=P_post2+P_min*eye(N);
    for o=1:16 
        temp=(i-1)*16+o; 
        MSE2(temp,u)=(norm(D(:,:,i)-D_est,'fro'))^2;
        if temp>800
            MSE2(temp,u)=MSE2(temp-1,u);
        end
       
    end 
end


%% EKF (Revised)

g_est = log(abs(h(:,1)) + 1e-8);  % Avoid log(0)
h_est3 = exp(g_est);             % Initial h estimate
P_post3 = weight_P * eye(N);     % Posterior covariance

MSE_h3(1,u) = norm(h(:,1) - h_est3)^2;
h_est_label3(:,1) = h_est3;

for i = 2:I
    % Q for g-space, small process noise
    Q_ekf = q * eye(N);

    % Jacobian H = C * diag(exp(g_est))
    H_ekf = C(:,:,i) * diag(exp(g_est));

    % EKF predict step
    P_pri3 = P_post3 + Q_ekf;

    % EKF update
    S = H_ekf * P_pri3 * H_ekf' + R * eye(M);
    K3 = P_pri3 * H_ekf' / S;

    y_pred = C(:,:,i) * exp(g_est);
    innovation = Y(:,i) - y_pred;

    % Update g_est (state)
    g_est = g_est + K3 * real(innovation);  % Approx: ignore phase wrapping

    % Update covariance
    P_post3 = (eye(N) - K3 * H_ekf) * P_pri3;

    % Update h and record
    h_est3 = exp(g_est);
    h_est_label3(:,i) = h_est3;
    MSE_h3(i,u) = norm(h(:,i) - h_est3)^2;
end


%% Unconstrainted KF
      
xn_est=zeros(2*N,1);
h_est4=zeros(N,1);      %channel estimation
P_pri4=zeros(2*N,2*N);      %prior covariance
P_post4=zeros(2*N,2*N);      %posterior covariance
K4=zeros(2*N,2*M);           %Kalman filter


%kalman filter algorithm initialization
P_post4=weight_P*eye(2*N);

average_value = mean(h(:,1), 'all')/2;
%h_est4 = zeros(N,1);
h_est4 =0.5*h(:,1);
%h_est4 = ones(N,1)*average_value;
xn_est=[real(h_est4);imag(h_est4)];
h_est_label4(:,1,u)=h_est4;


MSE_h4(1,u)=(norm(h(:,1)-h_est4))^2;

%rho=1;
P_min=1e-12;
qq=0.1;
q_ncnkf_scale = (sigma_W / 1e-5)^2;   
%q_ncnkf_scale=1;
for i=2:I
        Q1(:,:,i)=qq*(1/6)*speed_max^2*time_slot^2*[abs( h_est4); (2*pi/lambda)*ones(N,1)]*[abs( h_est4)', (2*pi/lambda)*ones(1,N)];
        B(:,:,i)=[diag(real( h_est4)),-diag(imag( h_est4));diag(imag( h_est4)),diag(real( h_est4))];

        Q4=q_ncnkf_scale*B(:,:,i)*Q1(:,:,i)*(B(:,:,i))'; %Q matrix
      % Q4=B(:,:,i)*Q1(:,:,i)*(B(:,:,i))'; %Q matrix
     
        P_pri4=(P_post4/rho)+Q4;
        K4=P_pri4*(C1(:,:,i))'/((C1(:,:,i))*P_pri4*(C1(:,:,i))'+R*eye(2*M));
       
        xn_est=xn_est+K4*(Y3(:,i)-(C1(:,:,i))*xn_est);
       
        xn_est_label(:,i)=xn_est;  %record the g_est
        for q2=1:N
            h_est4(q2)=xn_est(q2)+1i*xn_est(q2+N);
        end
       
        h_est_label4(:,i,u)=h_est4;
        P_post4=(eye(2*N)-K4*(C1(:,:,i)))*P_pri4;
        P_post4=P_post4+P_min*eye(2*N);
        
        MSE_h4(i,u)=(norm(h(:,i)-h_est4,2))^2;
       
end      
%% Normal KF
h_est11=zeros(N,1);      %channel estimation
P_pri11=zeros(N,N);      %prior covariance
P_post11=zeros(N,N);      %posterior covariance
K11=zeros(N,M);           %Kalman filter
%kalman filter algorithm
P_post11=weight_P*eye(N);
average_value = mean(h(:,1), 'all')/2;
%h_est1 = ones(N,1)*average_value;
h_est11 =0.5*h(:,1);
%h_est1 = zeros(N,1);


MSE11(1,u)=(norm(h(:,1)-h_est11))^2;
h_est_label11(:,1)=h_est11;
for i=2:I

        Q111=q*eye(N);
            
        P_pri11=P_post11+Q111;
        K11=P_pri11*(C(:,:,i))'/((C(:,:,i))*P_pri11*(C(:,:,i))'+R*eye(M));

        h_est11=h_est11+K11*(Y(:,i)-(C(:,:,i))*h_est11);
        h_est_label11(:,i)=h_est11;
        P_post11=(eye(N)-K11*(C(:,:,i)))*P_pri11;
       

        MSE11(i,u)=(norm(h(:,i)-h_est11,2))^2;




end
    
end



MSE1_1=mean(MSE1,2);
MSE2_2=mean(MSE2,2);
MSE3_3=mean(MSE_h3,2);
MSE4_4=mean(MSE_h4,2);
MSE5_5=mean(MSE11,2);
hnorm0=mean(norm_h,2);
hnorm00=mean(norm_h.^2,2);
NMSE1_1=MSE1_1./hnorm00;
NMSE2_2=MSE2_2./hnorm00;
NMSE3_3=MSE3_3./hnorm00;
NMSE4_4=MSE4_4./hnorm00;
NMSE5_5=MSE5_5./hnorm00;
SNR_W=mean(SNR,2);

filter_length = N; 
fir_coeff = ones(1, filter_length) / filter_length;


SNR_W_filtered = filter(fir_coeff, 1, SNR_W);


SNR_W_filtered_dB = 10 * log10(SNR_W_filtered);

%%
q3=1:1:I;

figure;
plot(q3, SNR_W_filtered_dB, 'r-', 'LineWidth', 1.5);
xlabel('Time');
ylabel('SNR (dB)');
title('FIR-filtered SNR (dB)');


pause

figure;
semilogy(q3,MSE1_1);
hold on
semilogy(q3,MSE2_2);
semilogy(q3,MSE3_3);
semilogy(q3,MSE4_4);
semilogy(q3,MSE5_5);
ylabel('Semilog MSE'); 
xlabel('Time'); 
%ylim([0 10]);
xlim([0 1000])
title('MSE Comparison with Time');
legend('KF','[Ref]','EKF','NCNKF','Standard KF');

hold off

pause



% figure;
% semilogy(q3,NMSE1_1);
% hold on
% %semilogy(q3,NMSE3_3);
% semilogy(q3,NMSE4_4);
% semilogy(q3,CRLB1);
% semilogy(q3,CRLB4);
% ylabel('Semilog NMSE'); 
% xlabel('Time'); 
% %ylim([0 10]);
% %xlim([0 1000])
% title('NMSE and CRLB with Time');
% legend('KF','NCNKF','CRLB1','CRLB4');
% 
% hold off
% 
% pause


