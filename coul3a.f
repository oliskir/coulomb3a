      PROGRAM TEST

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      logical l1
      character ftype
      character*25 histname 
      
      parameter(pi=3.141592654)
      parameter(Nmax=10000000)
      parameter(dt=1e-2) ! time step in units of 1e-23 sec
      
      double precision Malpha
      parameter(Malpha=3727.37897610664)
      double precision P(3,2),X(3,2),F(3,2)
     +     ,Pt(3,2),Xt(3,2)
     +     ,P0(3,2),X0(3,2),E0(3),Ecoulomb

      double precision Px_recoil,Py_recoil,Vx_recoil,Vy_recoil
     +     ,Vx_a2,Vy_a2,Vx_a3,Vy_a3,V_a1,V_a2,V_a3,cos12,cos13
     +     ,ang12,ang13         


      R2 = 4.5
      E3alpha = 7.275 
      ExC12 = 16.62

      R1 = 15.
      ExBe8 = 2.0 

      open(31,file='output/angles_e2000keV_r15fm.dat',
     +      status='unknown')
      
      open(32,file='output/energies_e2000keV_r15fm.dat',
     +      status='unknown')
      

      do ith=0,10

         theta = (dble(ith)*10 - 5.)*pi/180.
      
*     Positions at t=0
*     Alpha 1
         X0(1,1) = -R1
         X0(1,2) = 0.
*     Alpha 2
         X0(2,1) = R2/2.*cos(theta)
         X0(2,2) = R2/2.*sin(theta)
*     Alpha 3
         X0(3,1) = -R2/2.*cos(theta)
         X0(3,2) = -R2/2.*sin(theta)
         
         
*     Coulomb energy at t=0
         r12 = sqrt((X0(1,1)-X0(2,1))**2.+(X0(1,2)-X0(2,2))**2.)
         r13 = sqrt((X0(1,1)-X0(3,1))**2.+(X0(1,2)-X0(3,2))**2.)
         r23 = sqrt((X0(2,1)-X0(3,1))**2.+(X0(2,2)-X0(3,2))**2.)
         V1 = 1.4*2*2*(1./r12+1./r13)
         V23 = 1.4*2*2/r23
         
*     a1 and Be8 will share V1 while a2 and a3 will
*     share V23 after breakup.
         
*     Energy of a1 at time of 8Be breakup (R=10fm)
         Ea1 = 2./3.*( ExC12 - E3alpha - ExBe8 - V1 )
         if (Ea1.lt.0.) then
            write(6,*) 'Ea1 negative!'
            read(*,*)
         endif
         
*     Momentum
         Pxa1 = -sqrt(2.*Malpha*Ea1)
         Pya1 = 0.
         
*     Energy of Be8 at time of breakup
         EBe8 = 0.5*Ea1
*     Speed of Be8
         Vbe8 = sqrt(2.*EBe8/(2.*Malpha))
         
*     Energy of a2 and a3 in recoil system at distance R2
         Ea23 = 1./2.*( ExBe8+0.092 - V23 )
         if (Ea23.lt.0.) then
            write(6,*) 'Ea23 negative!'
            read(*,*)
         endif
         
*     Velocity of a2 and a3 in recoil system (8Be)
         w = sqrt(2.*Ea23/Malpha)
*     Transform to c.m.
         Vxa2 = Vbe8 + w*cos(theta)
         Vya2 = w*sin(theta)
         Vxa3 = Vbe8 - w*cos(theta)
         Vya3 = -w*sin(theta)
         
*     Momenta of a2 and a3 in c.m.
         Pxa2 = Malpha*Vxa2
         Pya2 = Malpha*Vya2
         Pxa3 = Malpha*Vxa3
         Pya3 = Malpha*Vya3

            
*     Momenta at t=0
*     Alpha 1
         P0(1,1) = Pxa1
         P0(1,2) = Pya1
*     Alpha 2
         P0(2,1) = Pxa2
         P0(2,2) = Pya2
*     Alpha 3
         P0(3,1) = Pxa3
         P0(3,2) = Pya3
*
*     Kinetic energies:
         do k=1,3
            E0(k) = (P0(k,1)**2.+P0(k,2)**2.)/(2.*Malpha)
         enddo
         
*     Total energy
         Etot0 = E0(1)+E0(2)+E0(3)+V1+V23
         
         write(6,*) 'Initial conditions:'
         write(6,*) 'x1,y1:',X0(1,1),X0(1,2)
         write(6,*) 'x2,y2:',X0(2,1),X0(2,2)
         write(6,*) 'x3,y3:',X0(3,1),X0(3,2)
         write(6,*) 'E1,E2,E3:',E0(1),E0(2),E0(3)
         write(6,*) 'Etot',Etot0
         

      
*
*     Numerical solution to equation of motion
*
         do k=1,3
            do i=1,2
               X(k,i) = X0(k,i)
               P(k,i) = P0(k,i)
            enddo
         enddo
         
      
         do istep=1,Nmax
         
*     distances:
            r12 = sqrt((X(1,1)-X(2,1))**2.+(X(1,2)-X(2,2))**2.)
            r13 = sqrt((X(1,1)-X(3,1))**2.+(X(1,2)-X(3,2))**2.)
            r23 = sqrt((X(2,1)-X(3,1))**2.+(X(2,2)-X(3,2))**2.)
            
*     Force on alpha 1
            F(1,1) = 1.4*2*2/r12**3.*(X(1,1)-X(2,1))
     +           + 1.4*2*2/r13**3.*(X(1,1)-X(3,1))
            F(1,2) = 1.4*2*2/r12**3.*(X(1,2)-X(2,2))
     +           + 1.4*2*2/r13**3.*(X(1,2)-X(3,2))
*     Force on alpha 2
            F(2,1) = 1.4*2*2/r12**3.*(X(2,1)-X(1,1))
     +           + 1.4*2*2/r23**3.*(X(2,1)-X(3,1))
            F(2,2) = 1.4*2*2/r12**3.*(X(2,2)-X(1,2))
     +           + 1.4*2*2/r23**3.*(X(2,2)-X(3,2))
*     Force on alpha 3
            F(3,1) = 1.4*2*2/r13**3.*(X(3,1)-X(1,1))
     +           + 1.4*2*2/r23**3.*(X(3,1)-X(2,1))
            F(3,2) = 1.4*2*2/r13**3.*(X(3,2)-X(1,2))
     +           + 1.4*2*2/r23**3.*(X(3,2)-X(2,2))
            
*     Calculate new momenta
            do k=1,3
               do i=1,2
                  Pt(k,i) = P(k,i) + F(k,i)*3.0*dt
               enddo
            enddo
            
*     Use new momenta to calculate new positions
            do k=1,3
               do i=1,2
                  Xt(k,i) = X(k,i) + Pt(k,i)/Malpha*3.0*dt
               enddo
            enddo
            
*     Update
            do k=1,3
               do i=1,2
                  X(k,i) = Xt(k,i)
                  P(k,i) = Pt(k,i)
               enddo
            enddo
            
         
            IF (MOD(istep,1000000).EQ.0) THEN
               WRITE(6,*) 'istep ... ',istep,ith
c     write(6,*) real(X(1,1)),real(X(1,2))
c     write(6,*) real(X(2,1)),real(X(2,2))
c     write(6,*) real(X(3,1)),real(X(3,2))
            ENDIF
                 
         
         enddo



*     Calculate end energies:
         Ea1 = (P(1,1)**2.+P(1,2)**2.)/(2.*Malpha)
         Ea2 = (P(2,1)**2.+P(2,2)**2.)/(2.*Malpha)
         Ea3 = (P(3,1)**2.+P(3,2)**2.)/(2.*Malpha)
         Etot = Ea1+Ea2+Ea3
         
         write(6,*) 'Final positions:'
         write(6,*) real(X(1,1)),real(X(1,2))
         write(6,*) real(X(2,1)),real(X(2,2))
         write(6,*) real(X(3,1)),real(X(3,2))
         write(6,*) 'Final energies (Ea1,Ea2,Ea3):'
         write(6,*) real(Ea1),real(Ea2),real(Ea3)
         write(6,*) 'Total energi before and after:'
         write(6,*) real(Etot0),real(Etot)
         
         
*     Energies
         write(32,*) real(theta*180./pi),',',real(Ea1),',',real(Ea2)
     +        ,',',real(Ea3)

*
*     Calculate angles of alphas at infinity:
*
*     Momentum of a2+a3
         Px_recoil = P(2,1) + P(3,1)
         Py_recoil = P(2,2) + P(3,2)
*     Velocity of recoil system
         Vx_recoil = Px_recoil/(2.*Malpha)
         Vy_recoil = Py_recoil/(2.*Malpha)
*     Transform velocities of a2 and a3 to recoil frame
         Vx_a2 = P(2,1)/Malpha - Vx_recoil
         Vy_a2 = P(2,2)/Malpha - Vy_recoil
         Vx_a3 = P(3,1)/Malpha - Vx_recoil
         Vy_a3 = P(3,2)/Malpha - Vy_recoil
*     Calculate angle relative to direction of a1
         V_a1 = sqrt((P(1,1)/Malpha)**2.+(P(1,2)/Malpha)**2.)
         V_a2 = sqrt(Vx_a2**2.+Vy_a2**2.)
         V_a3 = sqrt(Vx_a3**2.+Vy_a3**2.)
         cos12 = (P(1,1)/Malpha*Vx_a2 + P(1,2)/Malpha*Vy_a2)
     +        / (V_a1*V_a2)
         cos13 = (P(1,1)/Malpha*Vx_a3 + P(1,2)/Malpha*Vy_a3)
     +        / (V_a1*V_a3)
         ang12 = 180.-acos(cos12)*180./pi
         ang13 = 180.-acos(cos13)*180./pi
         
         write(31,*) real(theta*180./pi),',',real(ang12),',',real(ang13)
     +        ,',',real(ang12+ang13)
         

      enddo ! loop over angles
      close(31)
      close(32)

      END


  
