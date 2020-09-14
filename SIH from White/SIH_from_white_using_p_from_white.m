clear all
close all
imax=100;
m=.5122;%.5122; %mosquito per human
a=.21; %bitting frequiency (per day)
g=1/10;      %mosquito death rate (per day)
n=12;        %duration of sporogony in mopsquito
b=.5;        %transmission probability: mosquito to human
c=.5;       %transmission probability: human to mosquito
alpha=1/332;  %hypnozoites activation rate

r=1/60;        %rate of blood stage infection clearance
omega=1/425  ;   %hypnozoites death rate
Nval=8.5;   %Number of hypnozoites per infection
N=Nval;



s_h0(1)=.95;          %initial
i_h0(1)=.05;

S0=.95;
I0=.05;
H0=0;

S_m0=.98*m;
E_m0=0;
I_m0=.02*m;




options = odeset('RelTol', 1e-5);
% options1 = odeset(options,'NonNegative', 1:7);

N=min(Nval,imax);
for i=2:imax+1
    s_h0(i)=0;
    i_h0(i)=0;
end

tmax=5000;   %max simulation time

stp=100;   %discretization step
part=20;     %how many part for pice wise construction fro $p$
time=linspace(0,tmax,stp);

[t,y]=ode45(@white,time,[s_h0,i_h0,S_m0,E_m0,I_m0],options,[m a b c r alpha omega n g N],imax);

% s_h=y(:,1:imax+1);            %Assigning variable
i_h=y(:,imax+2:2*imax+2);

for i=1:length(time)
    pp(i)=i_h(i,1)/sum(i_h(i,:),2);  % getting $p$ from white
end

for i=1:part
    p(i)=mean(pp((stp/part)*i-((stp/part)-1):(stp/part)*i));   % taking the average between pice
end


for i=1:part
    Time{i}=linspace((tmax/part)*i-(tmax/part),(tmax/part)*i,stp/part);   % 
end

ini=[S0,I0,H0,S_m0,E_m0,I_m0];     %initials for SIH model

for i=1:part
    [t1,y1]=ode45(@SEI_clasic, Time{i},ini,options,[m a b c r alpha omega N n g p(i)],imax);

    ini=y1(end,:);
    pr_sih=y1(:,2);
    pl=plot(t1,pr_sih,'-.g','linewidth',2);

    hold on
end
set(get(get(pl(end),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
 legend('SIH')

pr_white=sum(i_h,2);

pl1=plot(t,pr_white,'r','linewidth',2);
xlabel('Time (Days)','FontWeight','bold','Fontsize',12)
ylabel('Fraction of infected individual','FontWeight','bold','Fontsize',12)


function ydot=white(t,y,parameter,imax)



m=parameter(1); a=parameter(2); b=parameter(3);c=parameter(4);r=parameter(5);
alpha=parameter(6);omega=parameter(7);n=parameter(8);g=parameter(9);N=parameter(10);


ds_h=zeros(imax+1,1);
di_h=zeros(imax+1,1);
ds_m=0;
de_m=0;
di_m=0;

%geometric distribution of number of hypnozoites per infection
mat=zeros(imax+1,imax+1);
prob=1/(N+1);
for i=1:imax
    for j=1:i
        mat(i,j)=prob*(1-prob).^(i-j);
        
    end
end

mat(imax+1,:)=1-sum(mat,1); % for balanceing the total population (1-column sum(all previous)


s_h=y(1:imax+1);            %Assigning variable
i_h=y(imax+2:2*imax+2);

s_m=y(2*imax+3);
e_m=y(2*imax+4);
i_m=y(2*imax+5);

lambda=m*a*b*i_m;
lambda_mat=mat*lambda;


ds_h(1)=-lambda*s_h(1)+r*i_h(1)+omega*s_h(2);
di_h(1)=-lambda*i_h(1)-r*i_h(1)+(omega+alpha)*i_h(2)+alpha*s_h(2)...
   +(lambda_mat(1,1)*(s_h(1)+i_h(1)));

for i=2:imax
    
    ds_h(i)=-lambda*s_h(i)+r*i_h(i)-(i-1)*(omega+alpha)*s_h(i)+i*omega*s_h(i+1);
    di_h(i)=-lambda*i_h(i)-r*i_h(i)-(i-1)*(omega+alpha)*i_h(i)+i*(omega+alpha)*i_h(i+1)+i*alpha*s_h(i+1)...
        +sum(lambda_mat(i,:)*(s_h+i_h));
    
end

ds_h(imax+1)=-lambda*s_h(imax+1)+r*i_h(imax+1)-imax*(omega+alpha)*s_h(imax+1);
di_h(imax+1)=-lambda*i_h(imax+1)-r*i_h(imax+1)-imax*(omega+alpha)*i_h(imax+1)...
    +sum(lambda_mat(imax+1,:)*(s_h+i_h));

ds_m=g*(s_m+e_m+i_m)-a*c*sum(i_h)*s_m-g*s_m;
de_m=a*c*sum(i_h)*s_m-e_m*(1/n)-g*e_m;
di_m=e_m*(1/n)-g*i_m;


ydot=[ds_h;di_h;ds_m;de_m;di_m];


end

function ydot=SEI_clasic(t,y,parameter,imax)
m=parameter(1);a=parameter(2);b=parameter(3);c=parameter(4);r=parameter(5);alpha=parameter(6);omega=parameter(7);
N=parameter(8);n=parameter(9);g=parameter(10);p=parameter(11);


total=0;


for i=2:imax
    k(i)=((N/(N+1))^(i)*(1/(N+1)));    %calculation the probability with different number of hypnozoites
    total=total+(i-1)*k(i);            %  k2+2*k3+3*k4...
end
 k(1)=1-sum(k(2:imax));

  
  
  
dS=zeros;
dI=zeros;
dH=zeros;
dS_m=zeros;
dE_m=zeros;
dI_m=zeros;

S=y(1);
I=y(2);
H=y(3);

S_m=y(4);
E_m=y(5);
I_m=y(6);

lambda=m*a*b*I_m;

dS=-lambda*S+omega*(k(1))*H+(p)*r*I;  
dI=lambda*(H+S)+alpha*(1+total)*H-r*I;
dH=-lambda*H-omega*(k(1))*H-alpha*(1+total)*H+(1-p)*r*I; 

dS_m=g*(S_m+I_m+E_m)-a*c*I*S_m-g*S_m;
dE_m=a*c*I*S_m-(1/n)*E_m-g*E_m;
dI_m=(1/n)*E_m-g*I_m;

ydot=[dS;dI;dH;dS_m;dE_m;dI_m];


end
