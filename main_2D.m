clear all
clc;

tic

%------- model 1------------------
% aR=[15.06,1.64,10.84,3.13,4.25];
% % Q=[30,12,20,15,15];   
% Q=[15,15,15,15,12];

%------- model 2------------------
% aR=[20,6,12,3,4];
% % Q=[14,8,10,8,12];
% Q=[40,50,60,40,40];

%------- model 3------------------
% aR=[14.4,6.5,9.0,2.25,3.0];
% Q=[7,4,5,4,6];
% Q=[8,8,8,8,6];

aR=[14.40,4.50,9.00,2.25,3.00];
% Q=[5,5,5,5,7];
Q=[30,15,20,15,12];

a(1:5)=aR(1:5).*complex(1,-1./Q(1:5));

%------- Janlu's models ------------------
tic
[Vray,Aray,Qray]=TTI_CRRT(a);
toc
figure(1)
plot(Vray(:,1,1),Vray(:,1,3),'-g','linewidth',1.5);
hold on
plot(Vray(:,2,1),Vray(:,2,3),'-b','linewidth',1.5);
hold on
plot(Vray(:,3,1),Vray(:,3,3),'-r','linewidth',1.5);
hold on
title('V-ray');
legend('qP','qSV','qSH');
xlabel('x-component (s/km)');ylabel('z-component (s/km)');

figure(2)
plot(Aray(:,1,1),Aray(:,1,3),'-g','linewidth',1.5);
hold on
plot(Aray(:,2,1),Aray(:,2,3),'-b','linewidth',1.5);
hold on
plot(Aray(:,3,1),Aray(:,3,3),'-r','linewidth',1.5);
hold on
title('A-ray');
legend('qP','qSV','qSH');
xlabel('x-component (s/km)');ylabel('z-component (s/km)');

figure(3)
plot(Qray(:,1,1),Qray(:,1,3),'-g','linewidth',1.5);
hold on
plot(Qray(:,2,1),Qray(:,2,3),'-b','linewidth',1.5);
hold on
plot(Qray(:,3,1),Qray(:,3,3),'-r','linewidth',1.5);
hold on
title('Q-ray');
legend('qP','qSV','qSH');
xlabel('x-component (s/km)');ylabel('z-component (s/km)');


function [Vray,Aray,Qray]=TTI_CRRT(a)
% a=(a11,a13,a33,a44,a66)............five complex-valued moduli;
% (zt0,ph0)........angles of the symmetric axis of a TTI-medium;

rad=pi/180;
zt1=0:1:90;
n1=length(zt1);
Vray(1:n1,1:3,1:3)=0;
Aray(1:n1,1:3,1:3)=0;
Qray(1:n1,1:3,1:3)=0;
aR=real(a);
Q=-aR./imag(a);
ZT(1:n1,1:3)=0;

for i=1:n1
    zheta=rad*zt1(i);
    
    for k=2
        switch k
            case 1
                [~,p]=phase_velocity(a,k,zheta,0);
                [Tau1,Tau2]=equation_Tau(a,p);

                if Tau1>Tau2
                    L=1;
                    if Q(1)==Q(2)&&Q(2)==Q(3)&&Q(3)==Q(4)
                        Q_4=1;
                        Fun_zt=@(ztl) equation_zt(a,k,L,Q_4,zheta,ztl);
                        [ztl,~]=fzero(Fun_zt,0);
                        
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(k)); cI=imag(c(k));
                        U_unify=[coeff1*sin(zheta)+coeff2*cos(zheta),0,coeff3*cos(zheta)+coeff4*sin(zheta)];
                        
                        
                        UR=1/(cR^2+cI^2)*(cR-cI/Q(1))*U_unify;
                        UI=-1/(cR^2+cI^2)*(cR/Q(1)+cI)*U_unify;
                        dir=UR/norm(UR);
                        Vray(i,1,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,1,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,1,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                        ZT(i,1)=UR(1)*UI(3)-UR(3)*UI(1);

                    else
                        Q_4=2;
                        F1=@(zt) Pequation1(a,zheta,zt);
                        [ztl,~]=fzero(F1,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(1)); cI=imag(c(1));
                        U_unify1=[coeff1*cR+coeff2*cI,0,coeff3*cR+coeff4*cI];
                        U_unify2=[coeff2*cR-coeff1*cI,0,coeff4*cR-coeff3*cI];
                        
                        UR=1/(cR^2+cI^2)*U_unify1;
                        UI=1/(cR^2+cI^2)*U_unify2;

                        dir=UR/norm(UR);
                        Vray(i,1,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,1,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,1,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    end
                    
                else
                    L=2;
                    if Q(1)==Q(2)&&Q(2)==Q(3)&&Q(3)==Q(4)
                        Q_4=1;
                        Fun_zt=@(ztl) equation_zt(a,k,L,Q_4,zheta,ztl);
                        [ztl,~]=fzero(Fun_zt,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(1)); cI=imag(c(1));
                        U_unify=[coeff1*sin(zheta)+coeff2*cos(zheta),0,coeff3*cos(zheta)+coeff4*sin(zheta)];
                        
                        
                        UR=1/(cR^2+cI^2)*(cR-cI/Q(1))*U_unify;
                        UI=-1/(cR^2+cI^2)*(cR/Q(1)+cI)*U_unify;
                        dir=UR/norm(UR);
                        Vray(i,1,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,1,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,1,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    else
                        Q_4=2;
                        F1=@(zt) Pequation2(a,zheta,zt);
                        [ztl,~]=fzero(F1,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(1)); cI=imag(c(1));
                        U_unify1=[coeff1*cR+coeff2*cI,0,coeff3*cR+coeff4*cI];
                        U_unify2=[coeff2*cR-coeff1*cI,0,coeff4*cR-coeff3*cI];
                        
                        UR=1/(cR^2+cI^2)*U_unify1;
                        UI=1/(cR^2+cI^2)*U_unify2;

                        dir=UR/norm(UR);
                        Vray(i,1,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,1,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,1,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    end
                end
                
            case 2
                [~,p]=phase_velocity(a,k,zheta,0);
                [Tau1,Tau2]=equation_Tau(a,p);
                
                if Tau1>Tau2
                    L=1;
                    if Q(1)==Q(2)&&Q(2)==Q(3)&&Q(3)==Q(4)
                        Q_4=1;
                        Fun_zt=@(ztl) equation_zt(a,k,L,Q_4,zheta,ztl);
                        [ztl,~]=fzero(Fun_zt,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(2)); cI=imag(c(2));
                        U_unify=[coeff1*sin(zheta)+coeff2*cos(zheta),0,coeff3*cos(zheta)+coeff4*sin(zheta)];
                        
                        
                        UR=1/(cR^2+cI^2)*(cR-cI/Q(1))*U_unify;
                        UI=-1/(cR^2+cI^2)*(cR/Q(1)+cI)*U_unify;
                        dir=UR/norm(UR);
                        Vray(i,2,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,2,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,2,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    else
                        Q_4=2;
                        F1=@(zt) SVequation1(a,zheta,zt);
                        [ztl,~]=fzero(F1,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(2)); cI=imag(c(2));
                        U_unify1=[coeff1*cR+coeff2*cI,0,coeff3*cR+coeff4*cI];
                        U_unify2=[coeff2*cR-coeff1*cI,0,coeff4*cR-coeff3*cI];
                        
                        UR=1/(cR^2+cI^2)*U_unify1;
                        UI=1/(cR^2+cI^2)*U_unify2;
                        dir=UR/norm(UR);
                        Vray(i,2,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,2,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,2,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    end
                    
                else
                    L=2;
                    if Q(1)==Q(2)&&Q(2)==Q(3)&&Q(3)==Q(4)
                        Q_4=1;
                        Fun_zt=@(ztl) equation_zt(a,k,L,Q_4,zheta,ztl);
                        [ztl,~]=fzero(Fun_zt,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(2)); cI=imag(c(2));
                        U_unify=[coeff1*sin(zheta)+coeff2*cos(zheta),0,coeff3*cos(zheta)+coeff4*sin(zheta)];
                        
                        UR=1/(cR^2+cI^2)*(cR-cI/Q(1))*U_unify;
                        UI=-1/(cR^2+cI^2)*(cR/Q(1)+cI)*U_unify;
                        dir=UR/norm(UR);
                        Vray(i,2,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,2,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,2,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    else
                        Q_4=2;
                        F1=@(zt) SVequation2(a,zheta,zt);
                        [ztl,~]=fzero(F1,0);
                        [coeff1,coeff2,coeff3,coeff4]=matrix_coeff(a,k,L,Q_4,zheta,ztl);
                        [c,~]=phase_velocity(a,k,zheta,ztl);
                        
                        cR=real(c(2)); cI=imag(c(2));
                        U_unify1=[coeff1*cR+coeff2*cI,0,coeff3*cR+coeff4*cI];
                        U_unify2=[coeff2*cR-coeff1*cI,0,coeff4*cR-coeff3*cI];
                        
                        UR=1/(cR^2+cI^2)*U_unify1;
                        UI=1/(cR^2+cI^2)*U_unify2;
                        dir=UR/norm(UR);
                        Vray(i,2,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                        Aray(i,2,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                        Qray(i,2,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    end
                end
            case 3
                L_zt=2*(Q(4)-Q(5))*sin(zheta)*cos(zheta)/(1+Q(4)*Q(5));
                zt=asinh(L_zt)/2;
                [c,~]=phase_velocity(a,k,zheta,zt);
                cR=real(c(3)); cI=imag(c(3));
                
                if Q(4)==Q(5)
                    
                    coeff1=(cR-cI/Q(5))/(cR^2+cI^2);
                    coeff2=-(cR/Q(5)+cI)/(cR^2+cI^2);
                    U_unify=[aR(5)*sin(zheta),0,aR(4)*cos(zheta)];
                    
                    UR=coeff1*U_unify;
                    UI=coeff2*U_unify; 
                    dir=UR/norm(UR);
                    
                    Vray(i,3,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                    Aray(i,3,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                    Qray(i,3,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    
                else
                    
                    coeff1=aR(5)*cosh(zt)*sin(zheta)+aR(5)/Q(5)*sinh(zt)*cos(zheta);
                    coeff2=aR(5)*sinh(zt)*cos(zheta)-aR(5)/Q(5)*cosh(zt)*sin(zheta);
                    coeff3=aR(4)*cosh(zt)*cos(zheta)-aR(4)/Q(4)*sinh(zt)*sin(zheta);
                    coeff4=-aR(4)*sinh(zt)*sin(zheta)-aR(4)/Q(4)*cosh(zt)*cos(zheta);
                    
                    U_unify1=[coeff1*cR+coeff2*cI,0,coeff3*cR+coeff4*cI];
                    U_unify2=[coeff2*cR-coeff1*cI,0,coeff4*cR-coeff3*cI];
                    
                    UR=1/(cR^2+cI^2).*U_unify1;
                    UI=1/(cR^2+cI^2).*U_unify2;
                    dir=UR/norm(UR);
                    
                    Vray(i,3,:)=(dir*(norm(UR)^2+norm(UI)^2))/(norm(UR));
                    Aray(i,3,:)=(dir*(norm(UI)))/(norm(UR)^2+norm(UI)^2);
                    Qray(i,3,:)=(dir*((norm(UR))^2-(norm(UI)^2)))/(2*norm(UR)*norm(UI));
                    
                end
                
        end
    end
    
    
end
end

