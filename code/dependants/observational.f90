! This module treats observational data
module observational   
    use, intrinsic :: iso_fortran_env, dp=>real64
    use derived_types, only: Star
    use files, only: getNumberOfLines, &
                     getNewUnit
    use astronomy, only: convertHoursToRad, &
                         convertDegreesToRad
    use criterion, only: splitDataForUVWvsUVW, &
                         splitDataForUVWvsMbol
    use plot_uvw, only: plotUVWvsUVW
    use plot_mbol, only: plotUVWvsMbol
    use math, only: getSD

    implicit none

    private                                                    
        
    public :: treatObservData

    character(len=*), parameter :: INPUT_PATH = './inputs&
                                                 &/observational.dat'

contains


    subroutine treatObservData(limogesCriterionIsUsed)

        logical, intent(in) :: limogesCriterionIsUsed
        type (Star), dimension(:), allocatable :: whiteDwarfs
        type (Star), dimension(:), allocatable :: sampleUvsV, &
                                                  sampleUvsW, &
                                                  sampleVvsW, &
                                                  sampleUvsMbol, &
                                                  sampleVvsMbol, &
                                                  sampleWvsMbol
        character(len = 11), dimension(:), allocatable :: RA_inHours, &
                                                          DEC_inDegrees
        integer :: numberOfWDs, &
                   unitInput, &
                   i

        numberOfWDs = getNumberOfLines(INPUT_PATH)
        
        allocate(whiteDwarfs(numberOfWDs))

        allocate(RA_inHours(numberOfWDs))
        allocate(DEC_inDegrees(numberOfWDs))

        print *, "Observational data from Limoges et al. 2015"
        print *, "Number of White Dwarfs:", numberOfWDs

        open(getNewUnit(unitInput), file = INPUT_PATH, status='old')
        
        do i = 1, numberOfWDs
            read(unitInput, *) whiteDwarfs(i)%distance, &
                               RA_inHours(i), &
                               DEC_inDegrees(i), &
                               whiteDwarfs(i)%motionInRA, &
                               whiteDwarfs(i)%motionInDEC, &
                               whiteDwarfs(i)%magnitude
        end do

        whiteDwarfs(:)%rightAscension = convertHoursToRad(RA_inHours(:))
        whiteDwarfs(:)%declination = convertDegreesToRad(DEC_inDegrees(:))

        call whiteDwarfs(:)%equatToGalact
        call whiteDwarfs(:)%galactToXYZ
        call whiteDwarfs(:)%equatToUVW

        if (limogesCriterionIsUsed) then
            print *, "Limoges criterion is used"
            call splitDataForUVWvsUVW(whiteDwarfs, &
                                      sampleUvsV, &
                                      sampleUvsW, &
                                      sampleVvsW)
            call splitDataForUVWvsMbol(whiteDwarfs, &
                                       sampleUvsMbol, &
                                       sampleVvsMbol, &
                                       sampleWvsMbol)
            call plotUVWvsUVW(sampleUvsV, &
                              sampleUvsW, &
                              sampleVvsW)
            call plotUVWvsMbol(sampleUvsMbol, &
                               sampleVvsMbol, &
                               sampleWvsMbol)
            print *, "Average velocity components:", &
                sampleUvsMbol(:)%vel(1) / size(sampleUvsMbol), &
                sampleVvsMbol(:)%vel(2) / size(sampleVvsMbol), &
                sampleWvsMbol(:)%vel(3) / size(sampleWvsMbol)
            print *, "Standart deviations:        ", &
                getSD(sampleUvsMbol(:)%vel(1)), &
                getSD(sampleVvsMbol(:)%vel(2)), &
                getSD(sampleWvsMbol(:)%vel(3))
        else 
            call plotUVWvsUVW(whiteDwarfs)
            call plotUVWvsMbol(whiteDwarfs)
            print *, "Average velocity components:", &
                whiteDwarfs(:)%vel(1) / size(whiteDwarfs), &
                whiteDwarfs(:)%vel(2) / size(whiteDwarfs), &
                whiteDwarfs(:)%vel(3) / size(whiteDwarfs)
            print *, "Standart deviations:        ", &
                getSD(whiteDwarfs(:)%vel(1)), &
                getSD(whiteDwarfs(:)%vel(2)), &
                getSD(whiteDwarfs(:)%vel(3))
        end if
    end subroutine treatObservData
end module
