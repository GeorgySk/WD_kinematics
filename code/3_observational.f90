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
    character(len=*), parameter :: INPUT_PATH_WITH_WD_CLASSES &
                                        = './inputs/observ_DA_nonDA.dat'
    character(len=*), parameter :: OUTPUT_PATH_DA &
                                        = './outputs/DA_nonDA/DA_kinem.dat'
    character(len=*), parameter :: OUTPUT_PATH_NONDA &
                                        = './outputs/DA_nonDA/nonDA_kinem.dat'
contains


    subroutine treatObservData(limogesCriterionIsUsed, &
                               splittingNonDAFromDA)

        logical, intent(in) :: limogesCriterionIsUsed
        logical, intent(in) :: splittingNonDAFromDA
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
                   unitOutputDA, &
                   unitOutputNonDA, &
                   i

        numberOfWDs = getNumberOfLines(INPUT_PATH_WITH_WD_CLASSES)
        
        allocate(whiteDwarfs(numberOfWDs))

        allocate(RA_inHours(numberOfWDs))
        allocate(DEC_inDegrees(numberOfWDs))

        print *, "Observational data from Limoges et al. 2015"
        print *, "Number of White Dwarfs:", numberOfWDs

        open(getNewUnit(unitInput), file = INPUT_PATH_WITH_WD_CLASSES, &
             status='old')
        
        do i = 1, numberOfWDs
            read(unitInput, *) whiteDwarfs(i)%distance, &
                               RA_inHours(i), &
                               DEC_inDegrees(i), &
                               whiteDwarfs(i)%motionInRA, &
                               whiteDwarfs(i)%motionInDEC, &
                               whiteDwarfs(i)%magnitude, &
                               whiteDwarfs(i)%spectralType
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
                              sampleVvsW, &
                              "observational")
            call plotUVWvsMbol(sampleUvsMbol, &
                               sampleVvsMbol, &
                               sampleWvsMbol, &
                               "observational")
            print *, "Average velocity components:", &
                sum(sampleUvsMbol(:)%vel(1)) / size(sampleUvsMbol), &
                sum(sampleVvsMbol(:)%vel(2)) / size(sampleVvsMbol), &
                sum(sampleWvsMbol(:)%vel(3)) / size(sampleWvsMbol)
            print *, "Standart deviations:        ", &
                getSD(sampleUvsMbol(:)%vel(1)), &
                getSD(sampleVvsMbol(:)%vel(2)), &
                getSD(sampleWvsMbol(:)%vel(3))
        else 
            call plotUVWvsUVW(whiteDwarfs, "observational")
            call plotUVWvsMbol(whiteDwarfs, "observational")
            print *, "Average velocity components:", &
                sum(whiteDwarfs(:)%vel(1)) / size(whiteDwarfs), &
                sum(whiteDwarfs(:)%vel(2)) / size(whiteDwarfs), &
                sum(whiteDwarfs(:)%vel(3)) / size(whiteDwarfs)
            print *, "Standart deviations:        ", &
                getSD(whiteDwarfs(:)%vel(1)), &
                getSD(whiteDwarfs(:)%vel(2)), &
                getSD(whiteDwarfs(:)%vel(3))
        end if

        if (splittingNonDAFromDA) then
            open(getNewUnit(unitOutputDA), file = OUTPUT_PATH_DA)
            open(getNewUnit(unitOutputNonDA), file = OUTPUT_PATH_NONDA)
            do i = 1, numberOfWDs
                if (whiteDwarfs(i)%spectralType == "DA") then
                    write(unitOutputDA, *) whiteDwarfs(i)%vel
                else if (whiteDwarfs(i)%spectralType == "nonDA") then
                    write(unitOutputNonDA, *) whiteDwarfs(i)%vel
                else 
                    print *, "Error: while trying to write data about DA and &
                              &nonDA kinematics program encountered something &
                              & else:", whiteDwarfs(i)%spectralType
                end if
            end do
        end if
    end subroutine treatObservData
end module
