function [A1,A2,A3,A4]=matrix_coeff(a,k,L,q4,zheta,zt)

zero=complex(0,0);one=complex(1,0);
deltz=eye(3,'like',complex(1,0));
n=[sin(zheta+1i*zt),0,cos(zheta+1i*zt)];

G11=a(1)*n(1)^2+a(4)*n(3)^2;
G22=a(5)*n(1)^2+a(4)*n(3)^2;
G33=a(4)*n(1)^2+a(3)*n(3)^2;
G13=(a(2)+a(4))*n(1)*n(3);
c(1)=sqrt(0.5*(G11+G33+sqrt((G11+G33)^2-4*(G11*G33-G13^2))));
c(2)=sqrt(0.5*(G11+G33-sqrt((G11+G33)^2-4*(G11*G33-G13^2))));
c(3)=sqrt(G22);

switch k
    case 1
        p1=n/c(1);
        G1(1,1)=a(1)*p1(1)*p1(1)+a(4)*p1(3)*p1(3);
        G1(1,3)=(a(2)+a(4))*p1(1)*p1(3);
        G1(2,2)=a(5)*p1(1)*p1(1)+a(4)*p1(3)*p1(3);
        G1(3,3)=a(4)*p1(1)*p1(1)+a(3)*p1(3)*p1(3);
        GM1=G1-deltz;
        
        B1(1)=GM1(1,1)*GM1(2,2);
        B1(2)=-GM1(1,3)*GM1(2,2);
        B1(3)=GM1(1,1)*GM1(3,3)-GM1(1,3)*GM1(1,3);
        B1(4)=GM1(2,2)*GM1(3,3);
        switch L
            case 1
                Pv1(1)=B1(2)/B1(1);  Pv1(2)=zero;  Pv1(3)=one;
                P_sqv1=sqrt(Pv1*Pv1'); P_g1=Pv1/P_sqv1;
                
            case 2
                Pv1(1)=one;  Pv1(2)=zero;  Pv1(3)=B1(2)/B1(4);
                P_sqv1=sqrt(Pv1*Pv1'); P_g1=Pv1/P_sqv1;
                
        end
    case 2
        p1=n/c(2);
        G1(1,1)=a(1)*p1(1)*p1(1)+a(4)*p1(3)*p1(3);
        G1(1,3)=(a(2)+a(4))*p1(1)*p1(3);
        G1(2,2)=a(5)*p1(1)*p1(1)+a(4)*p1(3)*p1(3);
        G1(3,3)=a(4)*p1(1)*p1(1)+a(3)*p1(3)*p1(3);
        GM1=G1-deltz;
        
        B1(1)=GM1(1,1)*GM1(2,2);
        B1(2)=-GM1(1,3)*GM1(2,2);
        B1(3)=GM1(1,1)*GM1(3,3)-GM1(1,3)*GM1(1,3);
        B1(4)=GM1(2,2)*GM1(3,3);
        switch L
            case 1
                Pv1(1)=B1(2)/B1(1);  Pv1(2)=zero;  Pv1(3)=one;
                P_sqv1=sqrt(Pv1*Pv1'); P_g1=Pv1/P_sqv1;
                
            case 2
                Pv1(1)=one;  Pv1(2)=zero;  Pv1(3)=B1(2)/B1(4);
                P_sqv1=sqrt(Pv1*Pv1'); P_g1=Pv1/P_sqv1;
                
        end
end

l_zt(1)=P_g1(1)*conj(P_g1(1));
l_zt(2)=P_g1(3)*conj(P_g1(3));
l_zt(3)=P_g1(1)*conj(P_g1(3))+P_g1(3)*conj(P_g1(1));

T1=l_zt(1)*a(1)+l_zt(2)*a(4); T1_R=real(T1); T1_I=imag(T1);
T2=l_zt(1)*a(5)+l_zt(2)*a(4); T2_R=real(T2); T2_I=imag(T2);
T3=l_zt(1)*a(4)+l_zt(2)*a(3); T3_R=real(T3); T3_I=imag(T3);
T4=l_zt(3)*(a(2)+a(4))/2; T4_R=real(T4); T4_I=imag(T4);

switch q4
    case 1
        A1=T1_R;
        A2=T4_R;
        A3=T3_R;
        A4=T4_R;
    case 2
        A1=cosh(zt)*(T1_R*sin(zheta)+T4_R*cos(zheta))+sinh(zt)*(-T1_I*cos(zheta)+T4_I*sin(zheta));
        A2=cosh(zt)*(T1_I*sin(zheta)+T4_I*cos(zheta))+sinh(zt)*(T1_R*cos(zheta)-T4_R*sin(zheta));
        A3=cosh(zt)*(T4_R*sin(zheta)+T3_R*cos(zheta))+sinh(zt)*(-T4_I*cos(zheta)+T3_I*sin(zheta));
        A4=cosh(zt)*(T4_I*sin(zheta)+T3_I*cos(zheta))+sinh(zt)*(T4_R*cos(zheta)-T3_R*sin(zheta));
        
end

end