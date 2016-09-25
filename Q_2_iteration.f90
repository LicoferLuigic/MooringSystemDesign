program main
    implicit none
    real:: vf, vs, h, dh
    real:: mp, Lx, k, l
    real:: theta0, alpha1
    real:: C1, C2, C3, C4
    integer:: i, j, n
    real:: e, pi, g, p, m0, mt, mb, S
    real:: Tc1, Tc2, Tci, Tcj, Tt1, Tt2, Tt3, Tt4, Fs, Ff, F0
    real:: d, alpha2, alphai, alphaj, theta1, theta2, theta3, theta4
    write(*,*) '输入初始配重质量。'
    read(*,*) mp
    write(*,*) '输入风速以及海水的流速。'
    read(*,*) vf, vs
    write(*,101)'vf=', vf, 'm/s, vs=', vs, 'm/s'
101 format (a4,f5.2,a8,f5.2,a3)  
    write(*,*) '输入海水深度。'
    read(*,*) h
    write(*,102) 'h=', h, 'm'
102 format (a3,f4.1,a1)
    write(*,*) '输入锚链的线密度、长度以及单节链环长度。'
    read(*,*) k, Lx, l
    n= ceiling((Lx/l)*1000)
    write(*,103) '锚链长', Lx, 'm，每米锚链重量', k, 'kg，共', n, '节链环。'
103 format (a6,f5.2,a15,f5.2,a6,i4,a8)    
    do
        call sub1(mp, vf, vs, h, Lx, k, l, n, d, alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4)
        if (theta0>=5) then
            mp= mp+0.1
            j= 0
        else
            mp= mp-0.1
            j= 1
        end if
        if ((alpha1<=16).AND.(theta0<=5)) exit
    end do
    if (j==0) then
        write(*,*) '配重球质量至少为', mp, '公斤。'
    else
        write(*,*) '配重球质量至多为', mp, '公斤。'
    end if
    write(*,110) '钢桶的倾斜角度为：', theta0, '度。'
110 format (a18,f9.4,a4)
    write(*,111) '各节钢管的倾斜角度由上至下依次为：', theta1/pi*180, '度，', theta2/pi*180, '度，', theta3/pi*180, '度，', theta4/pi*180, '度。'
111 format (a34/,4(f9.4,a4))
    write(*,112) '锚链方程为 y=',C1,'*cosh(x/', C1,'+', C2, ')-', C3
112 format (a13,f9.6,a8,f9.6,a1,f9.6,a2,f9.6)
    write(*,113) '浮标吃水深度', d, 'm。'
113 format(a12,f18.6,a3)
    write(*,114) '浮标的游动区域为以锚点为圆心，', dh, 'm为半径的圆形区域。'
114 format(a30,f5.2,a17)
end program main
    
subroutine sub1(mp, vf, vs, h, Lx, k, l, n, d, alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4)
    implicit none
    real,intent(in):: mp, vf, vs, h, Lx, k, l
    integer,intent(in):: n
    real,intent(out):: d, alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4
    real:: e, pi, g, p, m0, mt, mb, S
    real:: Tc1, Tc2, Tci, Tcj, Tt1, Tt2, Tt3, Tt4, Fs, Ff, F0
    real:: alpha2, alphai, alphaj
    real:: y=0, yi, yj, x=0, xi, xj, dv
    integer:: i=1, j=1
    e= 2.781
    pi= 3.14159
    g= 9.8
    p= 1.025*10**3
    m0= 1000
    mt= 10
    mb= 100
    S= 403
    d = ((m0+4*mt+mb+mp+Lx*k)*g)/(pi*p*g)
    Fs= 374*2*d*vs**2
    Ff= 0.625*2*(2-d)*vf**2
    F0= (m0+4*mt+mb+mp+Lx*k)*g
    Tt1= sqrt((F0-mt*g-m0*g)**2+(Ff+Fs)**2)
    theta1= atan((Ff+Fs)/(F0-mt*g-m0*g))
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
    end do
    i= 1
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
    return
end subroutine sub1