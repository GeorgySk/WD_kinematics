! This module treats observational data
module observational   
    use, intrinsic :: iso_fortran_env, dp=>real64
    use derived_types, only: Star, &
                             JaggedArray
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
                                        = './outputs/observ/&
                                            &DA_nonDA/DA_kinem.dat'
    character(len=*), parameter :: OUTPUT_PATH_NONDA &
                                        = './outputs/observ/&
                                            &DA_nonDA/nonDA_kinem.dat'
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
                                                  sampleWvsMbol, &
                                                  sampleDA, &
                                                  sampleNonDA
        character(len = 11), dimension(:), allocatable :: RA_inHours, &
                                                          DEC_inDegrees
        integer :: numberOfWDs, &
                   unitInput, &
                   unitOutputDA, &
                   unitOutputNonDA, &
                   i, &
                   counterDA, &
                   counterNonDA
        real(dp), dimension(3) :: sumOfDAVelocities, &
                                  sumOfNonDAVelocities


        type (JaggedArray), dimension(:), allocatable :: bins
        type (JaggedArray), dimension(:), allocatable :: binsNonDA
        real(dp), parameter :: MBOL_MIN = 5.75_dp
        real(dp), parameter :: MBOL_MAX = 20.75_dp
        real(dp), parameter :: MBOL_INC = 0.5_dp
        integer, parameter :: NUM_OF_BINS = int((MBOL_MAX - MBOL_MIN) &
                                                / MBOL_INC)
        integer :: binNumber
        integer, dimension(NUM_OF_BINS) :: wdInBinCounter = 0
        character(len = *), parameter :: MBOL_AVG_OBS_DA &
            = './outputs/observ/DA_nonDA/mbol_avg_da.dat'
        character(len = *), parameter :: MBOL_AVG_OBS_NONDA &
            = './outputs/observ/DA_nonDA/mbol_avg_nonda.dat'
        integer :: unitMbolAvgDA, unitMbolAvgNonDA
        character(len = *), parameter :: MBOL_AVG_FORMAT = '(7(f12.6,3x))'


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
            counterDA = 0
            counterNonDA = 0
            sumOfDAVelocities(:) = 0
            sumOfNonDAVelocities(:) = 0
            do i = 1, numberOfWDs
                if (whiteDwarfs(i)%spectralType == "DA") then
                    write(unitOutputDA, *) whiteDwarfs(i)%vel, &
                                           whiteDwarfs(i)%magnitude
                    counterDA = counterDA + 1
                else if (whiteDwarfs(i)%spectralType == "nonDA") then
                    write(unitOutputNonDA, *) whiteDwarfs(i)%vel, &
                                              whiteDwarfs(i)%magnitude
                    counterNonDA = counterNonDA + 1
                else 
                    print *, "Error: while trying to write data about DA and &
                              &nonDA kinematics program encountered something &
                              &else:", whiteDwarfs(i)%spectralType
                    stop
                end if
            end do

            allocate(sampleDA(counterDA))
            allocate(sampleNonDA(counterNonDA))

            counterDA = 0
            counterNonDA = 0
            do i = 1, numberOfWDs
                if (whiteDwarfs(i)%spectralType == "DA") then
                    counterDA = counterDA + 1
                    sampleDA(counterDA) = whiteDwarfs(i)
                else if (whiteDwarfs(i)%spectralType == "nonDA") then
                    counterNonDA = counterNonDA + 1
                    sampleNonDA(counterNonDA) = whiteDwarfs(i)
                else 
                    print *, "Error: while trying to write data about DA and &
                              &nonDA kinematics program encountered something &
                              &else:", whiteDwarfs(i)%spectralType
                    stop
                end if
            end do
            print *, "There are", counterDA, "of DA"
            print *, "There are", counterNonDA, "of nonDA"
            print *, "Average velocity components for DA:   ", &
                sum(sampleDA(:)%vel(1)) / size(sampleDA), &
                sum(sampleDA(:)%vel(2)) / size(sampleDA), &
                sum(sampleDA(:)%vel(3)) / size(sampleDA)
            print *, "Average velocity components for nonDA:", &
                sum(sampleNonDA(:)%vel(1)) / size(sampleNonDA), &
                sum(sampleNonDA(:)%vel(2)) / size(sampleNonDA), &
                sum(sampleNonDA(:)%vel(3)) / size(sampleNonDA)
            print *, "Standart deviations for DA:           ", &
                getSD(sampleDA(:)%vel(1)), &
                getSD(sampleDA(:)%vel(2)), &
                getSD(sampleDA(:)%vel(3))
            print *, "Standart deviations for nonDA:        ", &
                getSD(sampleNonDA(:)%vel(1)), &
                getSD(sampleNonDA(:)%vel(2)), &
                getSD(sampleNonDA(:)%vel(3))

            ! Filling bins
            allocate(bins(NUM_OF_BINS))
            do i = 1, size(sampleDA)
                if (sampleDA(i)%magnitude > MBOL_MIN .and. &
                        sampleDA(i)%magnitude < MBOL_MAX) then
                    binNumber = ceiling((sampleDA(i)%magnitude - MBOL_MIN) &
                                        / MBOL_INC)
                    wdInBinCounter(binNumber) = wdInBinCounter(binNumber) + 1
                end if
            end do
            do i = 1, NUM_OF_BINS
                if (wdInBinCounter(i) > 0) then
                    allocate(bins(i)%row(wdInBinCounter(i)))
                end if
            end do
            wdInBinCounter = 0
            do i = 1, size(sampleDA)
                if (sampleDA(i)%magnitude > MBOL_MIN .and. &
                        sampleDA(i)%magnitude < MBOL_MAX) then
                    binNumber = ceiling((sampleDA(i)%magnitude - MBOL_MIN) &
                                        / MBOL_INC)
                    wdInBinCounter(binNumber) = wdInBinCounter(binNumber) + 1
                    bins(binNumber)%row(wdInBinCounter(binNumber)) = sampleDA(i)
                end if
            end do
            open(getNewUnit(unitMbolAvgDA), file = MBOL_AVG_OBS_DA, status='old')
            do i = 1, NUM_OF_BINS
                if (allocated(bins(i)%row)) then
                    write(unitMbolAvgDA, MBOL_AVG_FORMAT) &
                        MBOL_MIN + MBOL_INC*(dfloat(i) - 0.5_dp), &
                        sum(bins(i)%row(:)%vel(1)) &
                            / size(bins(i)%row), &
                        sum(bins(i)%row(:)%vel(2)) &
                            / size(bins(i)%row), &
                        sum(bins(i)%row(:)%vel(3)) &
                            / size(bins(i)%row), &
                        getSD(bins(i)%row(:)%vel(1)), &
                        getSD(bins(i)%row(:)%vel(2)), &
                        getSD(bins(i)%row(:)%vel(3))
                end if
            end do

            wdInBinCounter(:) = 0

            ! Filling bins
            allocate(binsNonDA(NUM_OF_BINS))
            do i = 1, size(sampleNonDA)
                if (sampleNonDA(i)%magnitude > MBOL_MIN .and. &
                        sampleNonDA(i)%magnitude < MBOL_MAX) then
                    binNumber = ceiling((sampleNonDA(i)%magnitude - MBOL_MIN) &
                                        / MBOL_INC)
                    wdInBinCounter(binNumber) = wdInBinCounter(binNumber) + 1
                end if
            end do
            do i = 1, NUM_OF_BINS
                if (wdInBinCounter(i) > 0) then
                    allocate(binsNonDA(i)%row(wdInBinCounter(i)))
                end if
            end do
            wdInBinCounter = 0
            do i = 1, size(sampleNonDA)
                if (sampleNonDA(i)%magnitude > MBOL_MIN .and. &
                        sampleNonDA(i)%magnitude < MBOL_MAX) then
                    binNumber = ceiling((sampleNonDA(i)%magnitude - MBOL_MIN) &
                                        / MBOL_INC)
                    wdInBinCounter(binNumber) = wdInBinCounter(binNumber) + 1
                    binsNonDA(binNumber)%row(wdInBinCounter(binNumber)) = sampleNonDA(i)
                end if
            end do
            open(getNewUnit(unitMbolAvgNonDA), file = MBOL_AVG_OBS_NONDA, status='old')
            do i = 1, NUM_OF_BINS
                if (allocated(binsNonDA(i)%row)) then
                    write(unitMbolAvgNonDA, MBOL_AVG_FORMAT) &
                        MBOL_MIN + MBOL_INC*(dfloat(i) - 0.5_dp), &
                        sum(binsNonDA(i)%row(:)%vel(1)) &
                            / size(binsNonDA(i)%row), &
                        sum(binsNonDA(i)%row(:)%vel(2)) &
                            / size(binsNonDA(i)%row), &
                        sum(binsNonDA(i)%row(:)%vel(3)) &
                            / size(binsNonDA(i)%row), &
                        getSD(binsNonDA(i)%row(:)%vel(1)), &
                        getSD(binsNonDA(i)%row(:)%vel(2)), &
                        getSD(binsNonDA(i)%row(:)%vel(3))
                end if
            end do
        end if
    end subroutine treatObservData
end module
