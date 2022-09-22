program Asteroids               ! Program Name
use, intrinsic :: omp_lib       ! Use library Open Multiprocessing
implicit none



interface Force                 ! Interfacing the Gravitation force calculating function with the program
    function Force(P, Pa, G, m1, m2, V, omega) result(F)
        
        Real(8) :: P(2,2), Pa(2), P13(2), P23(2), V(2)
        Real(8) :: G, m1, m2, r13, r23, omega
        Real(8), dimension(2) :: F

    end function 
end interface


! Allocating memory to all required variables
Real(8) :: G, m1, m2, d, h, t, msun, dscale, tscale                            
integer :: n, i, samples, j 
Real(8) :: x1, y1, x2, y2, x3, y3, rseed, thetaseed
Real(8) :: P(2,2), Pa(2), F(2), V(2), Vmid(2), Fmid(2), Pmid(2)
Real(8) :: omega


! initial values
G = 6.67e-11 !* msun / ((dscale**3) * (tscale**2))      ! gravitational constant
msun = 1.989e30                                         ! mass of sun    
dscale = 752.71e9                                       ! distance scaler for non dimensionalizing
tscale = 86400                                          ! time scaler for non dimensionalizing (seconds in a day)
m1 = 1.989e30!/msun                                     ! reduced mass of sun
m2 = 1.898e27!/msun                                     ! reduced mass of jupiter
d = 752.71e9!/dscale                                    ! sun-jupiter distance

samples = 100000                                        ! number of iterations

! Sun jupiter coordinates
x1 = -m2 * d / (m1 + m2) 
x2 =  m1 * d / (m1 + m2) 
y1 = 0; y2 = 0


P(1,1) = x1
P(1,2) = y1
P(2,1) = x2
P(2,2) = y2

! angular velocity for co rotating frame
omega = (G*(m1+m2)/d**3)**0.5 

! coordinates of asteroid at L4
!x3 = ( x2 + x1) / 2 * 1.01
!y3 = (-x1 + x2) * (3**0.5) / 2

Pa(1) = x3
Pa(2) = y3

! Asteroid velocity
V = (/0, 0/)            

! Simulation Parameters
h = tscale*14               ! timestep
t = 0                       ! initial time
n = 26*10000                ! 26 days





open(unit = 10, file =  "Asteroid.txt", status = "replace", action = "readwrite" )      ! opening output file
    write(10,*) P(1,1), P(1,2)
    write(10,*) P(2,1), P(2,2)
    write(10,*) Pa(1) , Pa(2)
    
    call omp_set_num_threads(4)
    !$OMP DO
    
        out: do j = 1,samples
            call random_number(rseed)
            rseed = rseed * (d) + 6.96e10

            call random_number(thetaseed)
            thetaseed = thetaseed * 360


            Pa = (/rseed* cos(thetaseed), rseed * sin(thetaseed)/)

         
            V(1) = ((G*Msun/rseed)**0.5 - omega*d) * cos(thetaseed)
            V(2) = ((G*Msun/rseed)**0.5 - omega*d) * sin(thetaseed)

            F = Force(P, Pa, G, m1, m2, V, omega)
            
            in: do i = 1 , n                                 ! Euler Richardson Integrator

                F = Force(P, Pa, G, m1, m2, V, omega)

                Vmid = V + F * h/2.0

                Pmid = Pa + V * h/2.0

                Fmid = Force(P, Pmid, G, m1, m2, Vmid, omega)

                V = V + Fmid * h

                Pa = Pa + Vmid * h 

                t = t + h

                if ((Pa(1)**2+Pa(2)**2)**0.5 >= 1e12) then   ! Kick out asteroids that go beyond orbit
                cycle out
                endif

            end do in
            !call cpu_time(finish)
            !print '("Time = ",f6.3," seconds.")',finish-start
            write(10,*) Pa(1), Pa(2)

            
        !print *, j
        end do out
    !$OMP END DO
close(10)    

end program Asteroids

! Gravitional force finding function
function Force(P, Pa, G, m1, m2, V, omega) result(F)
    implicit none
    Real(8) :: P(2,2), Pa(2), P13(2), P23(2), V(2)
    Real(8) :: G, m1, m2, r13, r23, omega
    Real(8), dimension(2) :: F 

    P13 = Pa - P(1,:) 
    P23 = Pa - P(2,:) 

    r13 = dot_product(P13, P13)**0.5
    r23 = dot_product(P23, P23)**0.5
    
    F(1) =  -( G * m1 * P13(1) / r13**3 ) - ( G * m2 * P23(1)/r23**3 ) + ( Pa(1)*omega**2 ) + ( 2*omega*V(2) )
    F(2) =  -( G * m1 * P13(2) / r13**3 ) - ( G * m2 * P23(2)/r23**3 ) + ( Pa(2)*omega**2 ) - ( 2*omega*V(1) )

end function 