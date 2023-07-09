Program AtractorDeChen
  Implicit None
  real::h,x,y,z,k,L,M,t,limite
  real:: a,b,c

  !Constantes

  
 ! write(*,*)'Diga el valor de sigma'   
 ! read(*,*) sigma
  
 ! write(*,*)'Diga el valor de rho'   
  !read(*,*) rho
  
  !write(*,*)'Diga el valor de beta'   
  !read(*,*) beta

  a=5
  b=-10
  c=-0.38
  
    t=0
    h=0.003333
    limite=500
! Constantes iniciales
    x=1
    y=1                 
    z=2 

    Open(20, file='chen.dat')
    do while(t<=limite)
       call KLM(x,y,z,f,g,h,K,L,M)
       x=x+M
        y=y+K
        z=z+L
 
        t=t+h
        write(20,*)x,y,z
    enddo
    CLOSE(20)

    contains
    real function j(x,y,z) !derivative of x
        real::x,y,z
        j = a*x-y*z
    end function j

    real function f(x,y,z) !derivative of y
        real::x,y,z
        f = b*y+x*z
    end function f

    real function g(x,y,z) !derivative of z
        real::x,y,z
        g = c*z+(x*y)/3
    end function g

  

 subroutine KLM(x,y,z,f,g,h,K,L,M) !subroutine for Runge Kutta (order 4)
        real::x,f,g,y,z,k1,k2,k3,k4,h,K,L,l1,l2,l3,l4,M,m1,m2,m3,m4
        m1=h*j(x,y,z)
        k1=h*f(x,y,z)
        l1=h*g(x,y,z)

        m2=h*j(x+h/2,y+k1/2,z+l1/2)
        k2=h*f(x+h/2,y+k1/2,z+l1/2)
        l2=h*g(x+h/2,y+k1/2,z+l1/2)

        m3=h*j(x+h/2,y+k2/2,z+l2/2)
        k3=h*f(x+h/2,y+k2/2,z+l2/2)
        l3=h*g(x+h/2,y+k2/2,z+l2/2)

        m4 = h*j(x+h,y+k3,z+l3)
        k4 = h*f(x+h,y+k3,z+l3)
        l4 = h*g(x+h,y+k3,z+l3)

        K = (k1+k4)/6  + (k2+k3)/3 
        L = (l1+l4)/6 + (l2+l3)/3
        M = (m1+m4)/6 + (m2+m3)/3

    end subroutine KLM
  end program AtractorDeChen
  
