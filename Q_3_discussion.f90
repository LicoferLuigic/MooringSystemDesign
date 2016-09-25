program main
    real:: mp, vf, vs, h, Lx, k, l, d
    real:: alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4
    integer:: num, n
    write(*,*) '输入锚链的线密度、长度以及单节链环长度。'
    read(*,*) k, Lx, l
    n= ceiling((Lx/l)*1000)
    write(1,103) '锚链长', Lx, 'm，每米锚链重量', k, 'kg，共', n, '节链环。'
103 format (a6,f5.2,a15,f5.2,a6,i4,a8)    
    write(*,*) '输入配重球的质量'
    read(*,*) mp
    write(1,104) '配重球的质量为', mp, 'kg'
104 format (a14,f7.2,a2)
    write(*,*) '选择需要讨论的参数：1-风速，2-流速，3-水深'
    read(*,*) num
    select case (num)
    case (1)
        write(*,*) '给定流速与水深：'
        read(*,*) vs, h
        open(unit=1,action='write',status='old',file='D:\windvelocity.txt')
        do vf=0,36,1
            call sub1(mp, vf, vs, h, Lx, k, l, n, d, alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4)
            write(1,110) vf, theta0, alpha1, d, dh
        end do
    case (2)
        write(*,*) '给定风速与水深：'
        read(*,*) vf, h
        open(unit=2,action='write',status='old',file='D:\streamvelocity.txt')
        do vs=0,1.5,0.1
            call sub1(mp, vf, vs, h, Lx, k, l, n, d, alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4)
            write(2,110) vs, theta0, alpha1, d, dh
        end do
    case (3)
        write(*,*) '给定风速与流速：'
        read(*,*) vf, vs
        open(unit=3,action='write',status='old',file='D:\depth.txt')
        do h=16,20,0.5
            call sub1(mp, vf, vs, h, Lx, k, l, n, d, alpha1, theta0, dh, C1, C2, C3, C4, theta1, theta2, theta3, theta4)
            write(3,110) h, theta0, alpha1, d, dh
        end do
110     format (f10.3,f10.3,f10.3,f10.3,f10.3)
    end select
end program
    
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