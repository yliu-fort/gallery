function [pot_S22] = generate_potential_fun (A,Tw,h,rev_period,x,y)
shift_con=0;
% A=22;
% Tw=16.6;

phi_shift=shift_con*pi;
ORDER=2;
DEBUG=0;

f_p = 1/Tw;%peak frequency
t_focus=0;
% delta_f=1/fs*2;
Nfft=900;
no_seg=100;
g=9.81;  
t_shift = 5e-3/2.0;

syms x1

options = optimset('Display','off');

w_p=2*pi*f_p;

myfun=@(x1)g*x1*tanh(x1*h)-w_p^2;

k_p=fsolve(myfun,0.01,options);

%{
k_max=6.5*k_p;

k3=linspace(10*k_max/Nfft,k_max,Nfft);
k3=k3.';
kvecx=k3;
ko=k_p;
Akp=A*ko;
sigma_f = (kvecx<=ko)*0.07+(kvecx>ko)*0.09;
%}
% % create spectrum
%
%{
gamma=3;
f = (kvecx*g.*tanh(kvecx*h)).^0.5./(2*pi);
alpha=exp(-(f-f_p).^2./(2*sigma_f.^2*f_p^2));   
Hs = 1;
C=1-0.287*log(gamma);    %JS normalising factor
Spm=0.3125*Hs^2*f_p^4*f.^-5.*exp(-1.25*(f_p./f).^4);    %PM pdf
Sjs=C*Spm.*(gamma.^alpha); %JONSWAP pdf
df = gradient(f);
%}
k_w=0.005606*(12/Tw)^2;
k_max=2.5*k_p;
k3=linspace(k_max/Nfft,k_max,Nfft);
k3=k3.';
kvecx=k3;
ko=k_p;
% 
s_k=exp(-(k3-k_p).^2./(2*k_w^2));

% plot(k3,s_f/max(s_f))
% hold on
% plot(k3,s_k)

% figure
t=-rev_period*Tw*ones(1,length(x)) + t_shift;

% t=0*ones(1,length(x));
% s_f=s_f_1.';

if ORDER~=1
    
mu1=0;
mu2=0;

delete_po=zeros(sum(1:Nfft),1);

count1=1;

for ii=1:Nfft
    for jj=(ii-1)*Nfft+1:(ii-1)*Nfft+1+ii-1
        delete_po(count1)=jj;
        count1=count1+1;
    end
end

end
leng=length(k3);

%am = Sjs.*df./(sum(Sjs.*df))*Akp/k_p;
am=s_k/sum(s_k)*A;

an=am.';

wm=sqrt(k3*g.*tanh(k3*h));
wn=wm.';

if ORDER~=1

am_mat=repmat(am,leng,1);

an_mat=repelem(an,leng).';

wm_mat=repmat(wm,leng,1);

wn_mat=repelem(wn,leng).';

k1=repmat(k3,leng,1);
k2=repelem(k3.',leng).';

Dp=(wm_mat+wn_mat).^2-g*abs(k1+k2).*tanh(abs(k1+k2)*h);


Dm=(wm_mat-wn_mat).^2-g*abs(k1-k2).*tanh(abs(k1-k2)*h);

Dp(delete_po)=[];
Dm(delete_po)=[];
k1(delete_po)=[];
k2(delete_po)=[];
wn_mat(delete_po)=[];
wm_mat(delete_po)=[];
am_mat(delete_po)=[];
an_mat(delete_po)=[];

Bp=(wm_mat.^2+wn_mat.^2)./(2*g)-(wm_mat.*wn_mat)/(2.*g).*(1-(cos(mu1-mu2))./(tanh(abs(k1)*h).*tanh(abs(k2)*h)))...
    .*(((wm_mat+wn_mat).^2+g*abs(k1+k2).*tanh(abs(k1+k2)*h))./(Dp))...
    +((wm_mat+wn_mat)./(2*g*Dp)).*((wm_mat.^3)./((sinh(abs(k1)*h)).^2)+(wn_mat.^3)./((sinh(abs(k2)*h)).^2));
Bm=(wm_mat.^2+wn_mat.^2)./(2*g)+(wm_mat.*wn_mat)/(2.*g).*(1+(cos(mu1-mu2))./(tanh(abs(k1)*h).*tanh(abs(k2)*h)))...
    .*(((wm_mat-wn_mat).^2+g*abs(k1-k2).*tanh(abs(k1-k2)*h))./(Dm))...
    +((wm_mat-wn_mat)./(2*g*Dm)).*((wm_mat.^3)./((sinh(abs(k1)*h)).^2)-(wn_mat.^3)./((sinh(abs(k2)*h)).^2));

Ap=-((wm_mat.*wn_mat).*(wm_mat+wn_mat))./(Dp).*(1-(cos(mu1-mu2))./(tanh(abs(k1)*h).*tanh(abs(k2)*h)))...
    +(1./(2*Dp)).*((wm_mat.^3)./((sinh(abs(k1)*h)).^2)+(wn_mat.^3)./((sinh(abs(k2)*h)).^2));
Am=((wm_mat.*wn_mat).*(wm_mat-wn_mat))./(Dm).*(1+(cos(mu1-mu2))./(tanh(abs(k1)*h).*tanh(abs(k2)*h)))...
    +(1./(2*Dm)).*((wm_mat.^3)./((sinh(abs(k1)*h)).^2)-(wn_mat.^3)./((sinh(abs(k2)*h)).^2));

X5=(am_mat.*an_mat.*Bp);
X6=(am_mat.*an_mat.*Bm);

XP5=(am_mat.*an_mat.*Ap./(cosh(abs(k1+k2)*h)));
XP6=(am_mat.*an_mat.*Am./(cosh(abs(k1-k2)*h)));
% XP5=(am_mat.*an_mat.*Ap);
% XP6=(am_mat.*an_mat.*Am);
process_po=1;
% size_seg=0;
size_seg=ceil(length(t)/no_seg);
eta_diff=zeros(length(t),1);
eta_sum=zeros(length(t),1);
end

phi3=k3*x-wm*(t+t_focus);

    phi3=phi3+phi_shift;

    
eta_linear=am.'*cos(phi3);
%y_vec=y;
y_vec=min(eta_linear+eta_diff.'+eta_sum.',max(-5/k_p,y));
disp('computing elevation...') 
% rand 
if ORDER~=1
    if k_p*h>5
    h=5/k_p;
end
process_po=1;
pot_diff=zeros(length(t),1);
pot_sum=zeros(length(t),1);
for a=1:no_seg
    break_flag=0;
        if process_po+size_seg>length(x)
        process_po=length(x)-size_seg;
        break_flag=1;
        end
    phi1=k1*x(process_po:process_po+size_seg)-wm_mat*(t(process_po:process_po+size_seg)+t_focus)+phi_shift;
    phi2=k2*x(process_po:process_po+size_seg)-wn_mat*(t(process_po:process_po+size_seg)+t_focus)+phi_shift;
    
    eta_diff(process_po:process_po+size_seg)=X6.'*cos(phi1-phi2);
    
    eta_sum(process_po:process_po+size_seg)=X5.'*cos(phi1+phi2);
    pot_diff(process_po:process_po+size_seg)=(XP6.*cosh(abs(k1-k2))).'*(sin(phi1-phi2).*(y_vec(process_po:process_po+size_seg)+h));
    pot_sum(process_po:process_po+size_seg)=(XP5.*cosh(abs(k1+k2))).'*(sin(phi1+phi2).*(y_vec(process_po:process_po+size_seg)+h));
    process_po=process_po+size_seg+1;
        if break_flag==1
            break
        end
   disp(['eta ' num2str(a)])
end
eta_S22=eta_linear.'+eta_diff+eta_sum;

if ORDER~=2
    heta=imag(hilbert(eta_linear));
    D31=(eta_linear.^2+heta.^2).*eta_linear;
    D33=(eta_linear.^2-3*heta.^2).*eta_linear;
    C=sech(2*k_p*h);
    S31=-(3*(1+3*C+3*C^2+2*C^3))/(8*(1-C)^3);
    S33=-S31;
    k_p_33=3*k_p;
    eta_31=(S31*D31)*k_p^2;
    eta_33=(S33*D33)*k_p^2;

end
end
if ORDER==1
    eta_S22=eta_linear.';
end

if ORDER==3
    eta_S22=eta_linear.'+eta_diff+eta_sum+eta_31.'+eta_33.';
    if DEBUG==1
        eta_S22=eta_33.';
    end
end


disp('computing velocity potential...') 

%pot_linear=(g*am./((wm).*cosh(abs(k3)*h))).'*(sin(phi3).*((cosh(abs(k3)*(eta_S22+h).'))));
%pot_linear=(g*am./((wm).*cosh(abs(k3)*h))).'*(sin(phi3).*((cosh(abs(k3)*(max(0,y_vec+h)))))).*(eta_S22.'>y);
pot_linear=(g*am./((wm).*cosh(abs(k3)*h))).'*(sin(phi3).*((cosh(abs(k3)*(y_vec+h)))));

% pot_linear=(g*am./((wm))).'*(sin(phi3));

if ORDER~=1

%pot_diff=pot_diff.*(eta_S22>y.');
%pot_sum=pot_sum.*(eta_S22>y.');
pot_S22=pot_linear.'+pot_diff+pot_sum;

if ORDER~=2
    pot_31=g/(w_p)*imag(hilbert(eta_31));
    
    w_p33=sqrt(k_p_33*g*tanh(k_p_33*h));
    pot_33=g/(w_p33)*imag(hilbert(eta_33));
    
end


end
if ORDER==1
    pot_S22=pot_linear.';
end

if ORDER==3
    
    pot_S22=pot_linear.'+pot_diff+pot_sum+pot_31.'+pot_33.';
    if DEBUG==1
        pot_S22=pot_33.';
    end
end

% if ORDER==2
% plot(x,eta_S22)
% hold on
% plot(x,eta_linear)
% xlabel('X (m)')
% ylabel('Elevation (m)')
% legend('Second order','Linear')
% figure
% plot(x,pot_S22)
% hold on
% plot(x,pot_linear)
% xlabel('X (m)')
% ylabel('Velocity potential (m^2/s)')
% 
% end


% end