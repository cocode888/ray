function F=equation_zt(a,k,L,q4,zheta,zt)

zero=complex(0,0);one=complex(1,0);
deltz=eye(3,'like',complex(1,0));
aR=real(a);
Q=-aR./imag(a);
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

switch q4
    case 1
        kapa1=0;
        kapa2=0;
        kapa3=0;
        kapa=(1+1/Q(1)^2)*(l_zt(3)^2*(aR(2)+aR(4))^2-(l_zt(1)*aR(1)+...
            l_zt(2)*aR(4))*(l_zt(1)*aR(4)+l_zt(2)*aR(3)));
        
        F=(kapa1+kapa2)*cosh(2*zt)+kapa*sinh(2*zt)+(kapa2-kapa1)*cos(2*zheta)+kapa3*sin(2*zheta);
        
    case 2
        kapa1=l_zt(1)*l_zt(3)*aR(1)*aR(2)*(1/Q(1)-1/Q(2))+l_zt(1)*l_zt(3)*aR(1)*aR(4)*(1/Q(1)-1/Q(4))+l_zt(2)*l_zt(3)*aR(2)*aR(4)*(1/Q(4)-1/Q(2));
        kapa2=l_zt(1)*l_zt(1)*aR(1)*aR(4)*(1/Q(1)-1/Q(4))+l_zt(1)*l_zt(2)*aR(1)*aR(3)*(1/Q(1)-1/Q(3))+l_zt(2)*l_zt(2)*aR(3)*aR(4)*(1/Q(4)-1/Q(3));
        kapa3=l_zt(1)*l_zt(3)*aR(2)*aR(4)*(1/Q(2)-1/Q(4))+l_zt(2)*l_zt(3)*aR(2)*aR(3)*(1/Q(2)-1/Q(3))+l_zt(2)*l_zt(3)*aR(3)*aR(4)*(1/Q(4)-1/Q(3));
        kapa=l_zt(3)*l_zt(3)*aR(2)*aR(2)*(1+1/(Q(2)*Q(2)))+l_zt(3)*l_zt(3)*2*aR(2)*aR(4)*(1+1/(Q(2)*Q(4)))+...
            l_zt(3)*l_zt(3)*aR(4)*aR(4)*(1+1/(Q(4)*Q(4)))-l_zt(1)*l_zt(1)*aR(1)*aR(4)*(1+1/(Q(1)*Q(4)))-l_zt(1)*l_zt(2)*aR(1)*aR(3)*(1+1/(Q(1)*Q(3)))-...
            l_zt(1)*l_zt(2)*aR(4)*aR(4)*(1+1/(Q(4)*Q(4)))-l_zt(2)*l_zt(2)*aR(3)*aR(4)*(1+1/(Q(3)*Q(4)));
        
        F=(kapa1+kapa2)*cosh(2*zt)+kapa*sinh(2*zt)+(kapa2-kapa1)*cos(2*zheta)+kapa3*sin(2*zheta);
        
end

end