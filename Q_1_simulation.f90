program main
    implicit none
    real:: e, pi, g, p, h, m0, mt, mb, S
    real:: l, k, Lx, mp, vf, vs, Tc1, Tc2, Tci, Tcj, Tt1, Tt2, Tt3, Tt4, Fs, Ff, F0
    real:: d, alpha1, alpha2, alphai, alphaj, theta0, theta1, theta2, theta3, theta4
    real:: y=0, yi, yj, x=0, xi, xj, dh, dv, C1, C2, C3, C4, t1, t2
    integer:: i=1, n
    e= 2.781
    pi= 3.14159
    g= 9.8
    p= 1.025*10**3
    h= 18.0
    m0= 1000
    mt= 10
    mb= 100
    S= 403
    write(*,*) '��������Լ���ˮ�����١�'
    read(*,*) vf, vs
    write(*,101)'vf=', vf, 'm/s, vs=', vs, 'm/s'
101 format (a4,f5.2,a8,f5.2,a3)  
    write(*,*) '���뺣ˮ��ȡ�'
    read(*,*) h
    write(*,102) 'h=', h, 'm'
102 format (a3,f4.1,a1)
    write(*,*) '����ê�������ܶȡ������Լ������������ȡ�'
    read(*,*) k, Lx, l
    n= ceiling((Lx/l)*1000)
    write(*,103) 'ê����', Lx, 'm��ÿ��ê������', k, 'kg����', n, '��������'
103 format (a6,f5.2,a15,f5.2,a6,i4,a8)    
    write(*,*) '���������������'
    read(*,*) mp
    write(*,104) '�����������Ϊ', mp, 'kg'
104 format (a14,f7.2,a2)
!   ���㸡����ˮƽ����
    d = ((m0+4*mt+mb+mp+Lx*k)*g)/(pi*p*g)
    Fs= 374*2*d*vs**2
    Ff= 0.625*2*(2-d)*vf**2
    F0= (m0+4*mt+mb+mp+Lx*k)*g
!   ������θֹܵ�ƫб��������
    Tt1= sqrt((F0-mt*g-m0*g)**2+(Ff+Fs)**2)
    theta1= atan((Ff+Fs)/(F0-mt-m0*g))
    Tt2= sqrt((Tt1*cos(theta1)-mt*g)**2+(Tt1*sin(theta1))**2)
    theta2= atan((Tt1*sin(theta1))/(Tt1*cos(theta1)-mt*g))
    Tt3= sqrt((Tt2*cos(theta2)-mt*g)**2+(Tt2*sin(theta2))**2)
    theta3= atan((Tt2*sin(theta2))/(Tt2*cos(theta2)-mt*g))
    Tt4= sqrt((Tt3*cos(theta3)-mt*g)**2+(Tt3*sin(theta3))**2)
    theta4= atan((Tt3*sin(theta3))/(Tt3*cos(theta3)-mt*g))
    dh= sin(theta1)+sin(theta2)+sin(theta3)+sin(theta4)
    dv= cos(theta1)+cos(theta2)+cos(theta3)+cos(theta4)
    Tc2= sqrt((Tt4*cos(theta4)-(mb+mp)*g)**2+(Tt4*sin(theta4))**2)
    theta0= atan((Tt4*sin(theta4))/(Tt4*cos(theta4)-(mb+mp)*g))
    dh= dh+sin(theta0)
    dv= dv+cos(theta0)
    alpha2= 1.570795-theta0
    open(UNIT=1001,ACTION='write',STATUS='old',FILE='D:\cordination_x.txt')
    open(UNIT=1002,ACTION='write',STATUS='old',FILE='D:\cordination_y.txt')
    Tci= sqrt((Tc2*sin(alpha2)-k*l*g/1000)**2+(Tc2*cos(alpha2))**2)
    alphai= atan((Tc2*sin(alpha2)-k*l*g/1000)/(Tc2*cos(alpha2)))
    xi= l*(cos(alphai)+cos(alpha2))/1000
    yi= l*(sin(alphai)+sin(alpha2))/1000
    i= i+1
    do while (alphai>0.AND.i<n)
        x= x+xi
        y= y+yi
        Tcj= sqrt((Tci*sin(alphai)-k*l*g/1000)**2+(Tci*cos(alphai))**2)
        alphaj= atan((Tci*sin(alphai)-k*l*g/1000)/(Tci*cos(alphai)))
        xj= l*cos(alphaj)/1000
        yj= l*sin(alphaj)/1000
        xi= xj
        yi= yj
        Tci= Tcj
        alphai= alphaj
        i= i+1
        write(1001,'(f10.6)') x
        write(1002,'(f10.6)') y
    end do
    alpha1= alphai
    C1= (Tci*cos(alphai))/(k*g)
    C2= log(tan(alpha1)+1/cos(alpha1))
    C3= C1/cos(alpha1)
    C4= C1*tan(alpha1)
    y= h-d-dv
    x= C1*(acosh((y+C3)/C1)-C2)
    alpha1= alpha1/pi*180
    theta0= theta0/pi*180
    dh= dh+x
    write(*,110) '��Ͱ����б�Ƕ�Ϊ��', theta0, '�ȡ�'
110 format (a18,f9.4,a4)
    write(*,111) '���ڸֹܵ���б�Ƕ�������������Ϊ��', theta1/pi*180, '�ȣ�', theta2/pi*180, '�ȣ�', theta3/pi*180, '�ȣ�', theta4/pi*180, '�ȡ�'
111 format (a34/,4(f9.4,a4))
    write(*,112) 'ê������Ϊ y=',C1,'*cosh(x/', C1,'+', C2, ')-', C3
112 format (a13,f9.6,a8,f9.6,a1,f9.6,a2,f9.6)
!    write(*,*) x
    write(*,113) '�����ˮ���', d, 'm��'
113 format(a12,f8.6,a3)
    write(*,114) '������ζ�����Ϊ��ê��ΪԲ�ģ�', dh, 'mΪ�뾶��Բ������'
114 format(a30,f5.2,a17)
end program