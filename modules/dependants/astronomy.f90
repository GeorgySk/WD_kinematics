module astronomy

    use math, only: PI
    implicit none

    private :: RA_GPOLE, DEC_GPOLE, AUX_ANGLE
        ! Right ascension of Galactic pole
        real*8, parameter :: RA_GPOLE = 192.859508d0 * PI / 180.d0 
        ! Declination of Galactic pole
        real*8, parameter :: DEC_GPOLE = 27.128336d0 * PI /180.d0 
        ! Auxiliary angle
        real*8, parameter :: AUX_ANGLE = 122.932d0 * PI / 180.d0 
    public :: CTK, &
              TRAN_MATR, &
              convertGalacticToXYZ, &
              convertEquatorToGalact, &
              convertHoursToRad, &
              convertDegreesToRad
        ! Astronomical unit in km/s
        real*8, parameter :: CTK = 4.74047d0
        ! Transformation matrix from equatorial to galactic coordinates
        real*8, parameter, dimension(3, 3) :: TRAN_MATR &
            = reshape( (/-0.054875d0,-0.873437d0,-0.483835d0,&       
                          0.494109d0,-0.444829d0, 0.746982d0,&       
                         -0.867666d0,-0.198076d0, 0.455983d0/),&     
                         shape(TRAN_MATR), order = (/2,1/))

contains
    
    function convertGalacticToXYZ(r, l, b) result(coordinate)
        real*8, intent(in) :: r, l, b
        real*8, dimension(3) :: coordinate

        coordinate(1) = r * dcos(b) * dcos(l)
        coordinate(2) = r * dcos(b) * dsin(l)
        coordinate(3) = r * dsin(b)
    end function

    function convertHoursToRad(angleInHours) result(angleInRadians)
        character(len = 11), intent(in) :: angleInHours
        real*8 angleInRadians
        real*8 hours,minutes,seconds
        
        read(angleInHours(1:2), *) hours
        read(angleInHours(4:5), *) minutes
        read(angleInHours(7:11), *) seconds
        angleInRadians = (hours+(minutes+seconds/60.d0)/60.d0) * 15.d0 & 
                         * pi / 180.d0
    end function

    function convertDegreesToRad(angleInDegrees) result(angleInRadians)
        character(len = 11), intent(in) :: angleInDegrees
        real*8 angleInRadians
        real*8 degrees, arcmins, arcsecs
        
        read(angleInDegrees(1:2), *) degrees
        read(angleInDegrees(4:5), *) arcmins
        read(angleInDegrees(7:11), *) arcsecs
        angleInRadians = (degrees+(arcmins+arcsecs/60.d0)/60.d0) * pi / 180.d0
    end function

    subroutine convertEquatorToGalact(ra,dec,l,b)
        real*8,intent(in) :: ra,dec
        real*8,intent(out) :: l,b                 
        real*8 x,y

        !$b=arcsin\Big(cos(\delta)cos(\delta_{G})cos(\alpha-\alpha_{G})+\\+sin(\delta)sin(\delta_{G})\Big)$
        b = (dasin(dcos(dec)*dcos(DEC_GPOLE)*dcos(ra-RA_GPOLE)&
            + dsin(dec)*dsin(DEC_GPOLE))) !lattitude 
        x = dsin(dec) - dsin(b)*dsin(DEC_GPOLE)
        y = dcos(dec) * dsin(ra-RA_GPOLE) * dcos(DEC_GPOLE)
        l = datan(x/y) + AUX_ANGLE - PI/2.d0 !longitude    
        if (x.gt.0.d0 .AND. y.lt.0.d0) then
            l = l + PI
        else if (x.lt.0.d0 .AND. y.lt.0.d0) then
            l = l + PI
        else if (x.lt.0.d0 .AND. y.gt.0.d0) then
            l = l + 2.d0*PI
        else if (l .gt. 2.d0*PI) then 
            l = l - 2.d0*PI
        end if
    end subroutine

end module