module criterion
    use, intrinsic :: iso_fortran_env, dp=>real64
    use derived_types, only: Star
    implicit none

    private

    public :: splitDataForUVWvsUVW, &
              splitDataForUVWvsMbol

contains

    subroutine splitDataForUVWvsUVW(whiteDwarfs, &
                                    sampleUvsV, &
                                    sampleUvsW, &
                                    sampleVvsW)
        type (Star), dimension(:), intent(in) :: whiteDwarfs
        type (Star), dimension(:), allocatable, intent(out) :: sampleUvsV
        type (Star), dimension(:), allocatable, intent(out) :: sampleUvsW
        type (Star), dimension(:), allocatable, intent(out) :: sampleVvsW
        real(dp) :: highestCoordValue
        integer :: xHighestCounter = 0, &
                   yHighestCounter = 0, &
                   zHighestCounter = 0, &
                   counterUvsV = 1, &
                   counterUvsW = 1, &
                   counterVvsW = 1
        integer :: i

        ! Calculating number of WDs for each plot
        do i = 1, size(whiteDwarfs)
            highestCoordValue = maxval(abs(whiteDwarfs(i)%coords))
            if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(1))) &
                .lt. 1.0d-5) then
                xHighestCounter = xHighestCounter + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(2))) &
                .lt. 1.0d-5) then
                yHighestCounter = yHighestCounter + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(3))) &
                .lt. 1.0d-5) then
                zHighestCounter = zHighestCounter + 1
            end if
        end do
        
        allocate(sampleUvsV(zHighestCounter))
        allocate(sampleUvsW(yHighestCounter))
        allocate(sampleVvsW(xHighestCounter))

        ! Distributing WDs in samples for each plot
        do i = 1, size(whiteDwarfs)
            highestCoordValue = maxval(abs(whiteDwarfs(i)%coords))
            if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(1))) &
                .lt. 1.0d-5) then
                sampleVvsW(counterVvsW) = whiteDwarfs(i)
                counterVvsW = counterVvsW + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(2))) &
                .lt. 1.0d-5) then
                sampleUvsW(counterUvsW) = whiteDwarfs(i)
                counterUvsW = counterUvsW + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(3))) &
                .lt. 1.0d-5) then
                sampleUvsV(counterUvsV) = whiteDwarfs(i)
                counterUvsV = counterUvsV + 1
            end if
        end do
    end subroutine splitDataForUVWvsUVW


    subroutine splitDataForUVWvsMbol(whiteDwarfs, &
                                     sampleUvsMbol, &
                                     sampleVvsMbol, &
                                     sampleWvsMbol)
        type (Star), dimension(:), intent(in) :: whiteDwarfs
        type (Star), dimension(:), allocatable, intent(out) :: sampleUvsMbol
        type (Star), dimension(:), allocatable, intent(out) :: sampleVvsMbol
        type (Star), dimension(:), allocatable, intent(out) :: sampleWvsMbol
        real(dp) :: highestCoordValue
        integer :: xHighestCounter = 0, &
                   yHighestCounter = 0, &
                   zHighestCounter = 0, &
                   counterUvsMbol = 1, &
                   counterVvsMbol = 1, &
                   counterWvsMbol = 1
        integer :: i

        ! Calculating number of WDs for each plot
        do i = 1, size(whiteDwarfs)
            highestCoordValue = maxval(abs(whiteDwarfs(i)%coords))
            if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(1))) &
                .lt. 1.0d-5) then
                xHighestCounter = xHighestCounter + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(2))) &
                .lt. 1.0d-5) then
                yHighestCounter = yHighestCounter + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(3))) &
                .lt. 1.0d-5) then
                zHighestCounter = zHighestCounter + 1
            end if
        end do
        
        allocate(sampleUvsMbol(yHighestCounter + zHighestCounter))
        allocate(sampleVvsMbol(xHighestCounter + zHighestCounter))
        allocate(sampleWvsMbol(xHighestCounter + zHighestCounter))

        ! Distributing WDs in samples for each plot
        do i = 1, size(whiteDwarfs)
            highestCoordValue = maxval(abs(whiteDwarfs(i)%coords))
            if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(1))) &
                .lt. 1.0d-5) then
                sampleVvsMbol(counterVvsMbol) = whiteDwarfs(i)
                sampleWvsMbol(counterWvsMbol) = whiteDwarfs(i)
                counterVvsMbol = counterVvsMbol + 1
                counterWvsMbol = counterWvsMbol + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(2))) &
                .lt. 1.0d-5) then
                sampleUvsMbol(counterUvsMbol) = whiteDwarfs(i)
                sampleWvsMbol(counterWvsMbol) = whiteDwarfs(i)
                counterUvsMbol = counterUvsMbol + 1
                counterWvsMbol = counterWvsMbol + 1
            else if (abs(highestCoordValue - abs(whiteDwarfs(i)%coords(3))) &
                .lt. 1.0d-5) then
                sampleUvsMbol(counterUvsMbol) = whiteDwarfs(i)
                sampleVvsMbol(counterVvsMbol) = whiteDwarfs(i)
                counterUvsMbol = counterUvsMbol + 1
                counterVvsMbol = counterVvsMbol + 1
            end if
        end do
    end subroutine splitDataForUVWvsMbol

    
end module criterion