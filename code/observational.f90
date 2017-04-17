! This module treats observational data
module observational   

    use files, only: getNumberOfLines, &
                     getNewUnit
    use commons, only: VW_PATH, &
                       UW_PATH, &
                       UV_PATH, &
                       OUTPUT_FORMAT 
    use astronomy, only: convertHoursToRad, &
                         convertDegreesToRad, &
                         convertEquatorToGalact, &
                         convertGalacticToXYZ, &
                         CTK, &
                         TRAN_MATR
    use math, only: multiplyMatrixByVector

    implicit none

    private :: OBSERV_PATH, &
               MBOL_CLOUD_PATH, &
               MBOL_AVG_PATH, &
               MBOL_MIN, &
               MBOL_MAX, &
               MBOL_INC, &
               NUM_OF_BINS, &
               VEL_RAD, &
               fillRotationMatrix, &
               CLOUD_FORMAT, &
               MBOL_AVG_FORMAT
        character(len = *), parameter :: OBSERV_PATH = './inputs&
                                                        &/observational.dat'
        character(len = *), parameter :: MBOL_CLOUD_PATH = './outputs&
                                                            &/mbol_cloud.dat'
        character(len = *), parameter :: MBOL_AVG_PATH = './outputs&
                                                          &/mbol_avg.dat'
        ! Binning parameters for bolometric magnitude (Mbol)
        real*8, parameter :: MBOL_MIN = 5.75d0                                      
        real*8, parameter :: MBOL_MAX = 20.75d0
        real*8, parameter :: MBOL_INC = 0.5d0
        integer, parameter :: NUM_OF_BINS = int((MBOL_MAX - MBOL_MIN) &
                                                / MBOL_INC)
        real*8, parameter :: VEL_RAD = 0.d0  ! Radial velocity 
        character(len = *), parameter :: CLOUD_FORMAT = '(4(f12.6,3x))' 
        character(len = *), parameter :: MBOL_AVG_FORMAT = '(7(f12.6,3x))'


    public :: treatObservationalData

contains

    subroutine treatObservationalData

        integer :: numberOfWDs
        ! In this array we record velocity components for each WD
        ! in order to calculate standart deviation 
        real*8, dimension(:, :), allocatable :: velocityArray
        ! In this array we record velocity components for each WD
        ! in every bin of Mbol in order to calculate standart 
        ! deviations per corresponding bin 
        ! NOTE: 2nd dimension should be dynamic - numberOfWDsInBin
        real*8, dimension(:, :, :), allocatable :: velocityArrayForMbol 
        real*8, dimension(:, :), allocatable :: magnitudeArray
        ! I/O units
        integer :: unitVW, &  ! Outputs with velocities
                   unitUW, &
                   unitUV, &
                   ! Observational data from Limoges et al. 2015
                   unitObs, &  
                   ! Output with individual (Mbol, U, V, W) for each WD
                   unitCloud, &  
                   ! Output with average velocities in bins of Mbol
                   unitMbolAvg
        integer :: i, &
                   binNumber
        ! Input data
        real*8 :: distance, &
                  motionInRA, &  ! Motion in Right Ascension (RA)
                  motionInDEC, &  ! Motion in Declination (DEC)
                  magnitude
        ! Input data in special formats to be converted
        character(len = 11) RA_inHours, &
                            DEC_inDegrees
        ! Data to be obtained from RA_inHours and DEC_inDegrees         
        real*8 :: rightAscension, &
                  declination
        ! Galactic angles
        real*8 :: l, &  ! Longitude
                  b  ! Lattitude
        ! WD coordinates in X, Y, Z
        real*8, dimension(3) :: coordOfWD
        ! WD velocity to be converted
        real*8, dimension(3) :: vel_motion = (/VEL_RAD, 0.d0, 0.d0/)
        real*8, dimension(3, 3) :: rotationMatrix 
        ! Heliocentric velocity of WD
        real*8, dimension(3) :: vel_hel
        integer, dimension(NUM_OF_BINS) :: numberOfWDsInBin = 0
        real*8, dimension(3, NUM_OF_BINS) :: sumOfVelocitiesInBin = 0.d0
        real*8 :: highestCoord
        integer, dimension(3) :: counterOfWD=(/0, 0, 0/)
        real*8, dimension(3) :: vel_sum = (/0.d0, 0.d0, 0.d0/)
        real*8, dimension(3) :: averageVelocityInBin
        real*8, dimension(3) :: sumOfRestsSquared = (/0.d0, 0.d0, 0.d0/)
        real*8, dimension(3, NUM_OF_BINS) :: magnitudeSigma
        real*8, dimension(3) :: vel_avg, &
                                sigma

        !-----------------------------------------------------------------------                        

        print*, "Observational data from Limoges et al. 2015"
        
        numberOfWDs = getNumberOfLines(OBSERV_PATH)
        print*, "Number of White Dwarfs:", numberOfWDs
        
        allocate(velocityArray(3, numberOfWDs))
        allocate(velocityArrayForMbol(3, NUM_OF_BINS, numberOfWDs))
        allocate(magnitudeArray(NUM_OF_BINS, numberOfWDs))
    
        open(getNewUnit(unitVW), file = VW_PATH, status='old')
        open(getNewUnit(unitUW), file = UW_PATH, status='old')
        open(getNewUnit(unitUV), file = UV_PATH, status='old')
        open(getNewUnit(unitObs), file = OBSERV_PATH, status='old')
        open(getNewUnit(unitCloud), file = MBOL_CLOUD_PATH, status='old')
        open(getNewUnit(unitMbolAvg), file = MBOL_AVG_PATH, status='old')
    
        do i = 1, numberOfWDs

            read(unitObs, *) distance, &
                             RA_inHours, &
                             DEC_inDegrees, &
                             motionInRA, &
                             motionInDEC, &
                             magnitude

            binNumber = ceiling((magnitude - MBOL_MIN) / MBOL_INC)
    
            ! Сonverting angles from char to radians
            rightAscension = convertHoursToRad(RA_inHours)      
            declination = convertDegreesToRad(DEC_inDegrees)
    
            ! Converting angles to lattitude and longitude
            call convertEquatorToGalact(rightAscension, declination, l, b) 
    
            ! Converting lattitude and longitude to X,Y,Z
            coordOfWD = convertGalacticToXYZ(distance, l, b)
    
            ! Motion in right ascension with asterisk
            motionInRA = motionInRA * dcos(declination)          
    
            !$\begin{pmatrix}U\\V\\W\end{pmatrix}=B\begin{pmatrix}v_{r}\\k\mu^{*}_{a}/\pi\\k\mu_{\delta}/\pi\end{pmatrix}$
            vel_motion(2) = CTK * motionInRA * distance              
            vel_motion(3) = CTK * motionInDEC * distance
    
            !$B=T\cdot{A}$
            call fillRotationMatrix(rightAscension, declination, rotationMatrix)         
            ! U,V,W         
            vel_hel = multiplyMatrixByVector(TRAN_MATR, &
                            multiplyMatrixByVector(rotationMatrix,&
                                                   vel_motion)) 
            write(unitCloud, CLOUD_FORMAT) magnitude, vel_hel
    
            if (binNumber .le. NUM_OF_BINS) then
                numberOfWDsInBin(binNumber) = numberOfWDsInBin(binNumber) + 1
                sumOfVelocitiesInBin(:, binNumber) &
                    = sumOfVelocitiesInBin(:, binNumber) + vel_hel
                velocityArrayForMbol(:, binNumber, numberOfWDsInBin(binNumber))&
                    = vel_hel
                magnitudeArray(binNumber, numberOfWDsInBin(binNumber)) &
                    = magnitude
            end if
    
            highestCoord = maxval(abs(coordOfWD))
    
            ! Highest is X
            if (abs(highestCoord - abs(coordOfWD(1))) .lt. 1.0d-5) then            
                ! Nº of WDs that will be taken into account 
                ! for plotting V and W
                counterOfWD(1) = counterOfWD(1) + 1
    
                ! Summing V and W components
                vel_sum(2) = vel_sum(2) + vel_hel(2)  
                vel_sum(3) = vel_sum(3) + vel_hel(3)
    
                ! Filling arrays for calculating SD
                velocityArray(2, counterOfWD(3)+counterOfWD(1)) = vel_hel(2)
                velocityArray(3, counterOfWD(2)+counterOfWD(1)) = vel_hel(3)
        
                ! Outputs/vw.dat
                write(unitVW, OUTPUT_FORMAT) vel_hel(2), vel_hel(3)  
                
            ! Highest is Y
            else if (abs( highestCoord - abs(coordOfWD(2) )) .lt. 1.0d-5) then
    
                ! Nº of WDs that will be taken into account 
                ! for plotting U and W
                counterOfWD(2) = counterOfWD(2) + 1  
        
                ! Summing U and W components 
                vel_sum(1) = vel_sum(1) + vel_hel(1) ! U 
                vel_sum(3) = vel_sum(3) + vel_hel(3) ! W
        
                ! Filling arrays for calculating SD
                velocityArray(1, counterOfWD(2)+counterOfWD(3)) = vel_hel(1)
                velocityArray(3,counterOfWD(2)+counterOfWD(1)) = vel_hel(3)
        
                ! outputs/uw.dat
                write(unitUW, OUTPUT_FORMAT) vel_hel(1), vel_hel(3)   
                
            ! Highest is Z
            else if (abs( highestCoord - abs(coordOfWD(3) )) .lt. 1.0d-5) then
        
                ! Nº of WDs that will be taken into account 
                ! for plotting U and V
                counterOfWD(3) = counterOfWD(3) + 1  
        
                ! Summing U and V components 
                vel_sum(1) = vel_sum(1) + vel_hel(1) 
                vel_sum(2) = vel_sum(2) + vel_hel(2)
    
                ! Filling arrays for calculating SD             
                velocityArray(1, counterOfWD(2)+counterOfWD(3)) = vel_hel(1)
                velocityArray(2, counterOfWD(3)+counterOfWD(1)) = vel_hel(2)
        
                !outputs/uv.dat
                write(unitUV, OUTPUT_FORMAT) vel_hel(1), vel_hel(2)  
                
            else
                print*, "Error: couldn't determine highest coordinate."
                
            end if    
    
        end do
    
        do binNumber = 1, NUM_OF_BINS
            
            !NOTE: If there are no WDs in bin then it will be zero
            averageVelocityInBin = sumOfVelocitiesInBin(:, binNumber) &
                                   / dfloat(numberOfWDsInBin(binNumber))
            do i = 1, numberOfWDsInBin(binNumber)
                sumOfRestsSquared = sumOfRestsSquared &
                                    + (velocityArrayForMbol(:, binNumber, i) &
                                    - averageVelocityInBin) ** 2
            end do

            if (numberOfWDsInBin(binNumber) .ne. 1) then
                magnitudeSigma(:, binNumber) &
                    = (sumOfRestsSquared / dfloat(numberOfWDsInBin(binNumber))&
                       - 1.d0) ** 0.5d0
            else
                magnitudeSigma(:, binNumber) = 100.d0
            end if

            write(unitMbolAvg, MBOL_AVG_FORMAT) &
                MBOL_MIN + MBOL_INC*(dfloat(binNumber) - 0.5d0), &
                averageVelocityInBin, &
                magnitudeSigma(:, binNumber)

            sumOfRestsSquared = 0

        end do
    
        !calculating average U,V,W
        vel_avg(1) = vel_sum(1) / dfloat(counterOfWD(3) + counterOfWD(2))  
        vel_avg(2) = vel_sum(2) / dfloat(counterOfWD(3) + counterOfWD(1)) 
        vel_avg(3) = vel_sum(3) / dfloat(counterOfWD(2) + counterOfWD(1))
        
        do i = 1, counterOfWD(3) + counterOfWD(2)
            sumOfRestsSquared(1) = sumOfRestsSquared(1) &
                                   + (velocityArray(1, i) - vel_avg(1)) ** 2
        end do
        do i = 1, counterOfWD(3) + counterOfWD(1)
            sumOfRestsSquared(2) = sumOfRestsSquared(2) &
                                   + (velocityArray(2, i) - vel_avg(2)) ** 2
        end do
        do i = 1, counterOfWD(2) + counterOfWD(1)
            sumOfRestsSquared(3) = sumOfRestsSquared(3) &
                                   + (velocityArray(3, i) - vel_avg(3)) ** 2
        end do
    
        sigma(1) = (sumOfRestsSquared(1) &
                    / dfloat(counterOfWD(3) + counterOfWD(2))) ** 0.5d0
        sigma(2) = (sumOfRestsSquared(2) &
                    / dfloat(counterOfWD(3) + counterOfWD(1))) ** 0.5d0
        sigma(3) = (sumOfRestsSquared(3) &
                    / dfloat(counterOfWD(2) + counterOfWD(1))) ** 0.5d0
        
        write(6,*) 'Average relative to Sun:', vel_avg
        write(6,*) 'Sigmas:                 ', sigma

        deallocate(velocityArray)

    end subroutine
    !---------------------------------------------------------------------------
    

    subroutine fillRotationMatrix(RA, DEC, A)
        real*8, intent(in) :: RA, &
                              DEC
        real*8, dimension(3, 3), intent(out) :: A
        A(1,1) = dcos(RA) * dcos(DEC)  
        A(2,1) = dsin(RA) * dcos(DEC)  
        A(3,1) = dsin(DEC)        

        A(1,2) = - dsin(RA)
        A(2,2) = dcos(RA)
        A(3,2) = 0.d0

        A(1,3) = - dcos(RA)*dsin(DEC)
        A(2,3) = - dsin(RA)*dsin(DEC)
        A(3,3) = dcos(DEC)
    end subroutine

end module
