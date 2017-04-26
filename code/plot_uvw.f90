module plot_uvw
    use, intrinsic :: iso_fortran_env, dp=>real64
    use derived_types, only: Star
    use files, only: getNewUnit

    implicit none

    private

    public :: plotUVWvsUVW

    interface plotUVWvsUVW
            module procedure threeSamples_plotUVWvsUVW
            module procedure oneSample_plotUVWvsUVW
    end interface plotUVWvsUVW
    
    character(len = *), parameter :: VW_PATH = './outputs/observ&
                                                &/limoges/vw.dat'
    character(len = *), parameter :: UW_PATH = './outputs/observ&
                                                &/limoges/uw.dat'
    character(len = *), parameter :: UV_PATH = './outputs/observ&
                                                &/limoges/uv.dat'
    character(len = *), parameter :: UVW_PATH = './outputs/observ/no_crit/uvw.dat'
    character(len = *), parameter :: OUTPUT_FORMAT = '(2(f12.6,3x))'
    character(len = *), parameter :: UVW_FORMAT = '(3(f12.6,3x))'

contains


    subroutine oneSample_plotUVWvsUVW(sample)
        type (Star), dimension(:), intent(in) :: sample
        integer :: i
        integer :: unitUVW

        open(getNewUnit(unitUVW), file = UVW_PATH, status='old')
        do i = 1, size(sample)
            write(unitUVW, OUTPUT_FORMAT) sample(i)%vel(:)
        end do 
    end subroutine oneSample_plotUVWvsUVW


    subroutine threeSamples_plotUVWvsUVW(sampleUvsV, &
                                         sampleUvsW, &
                                         sampleVvsW)
        type (Star), dimension(:), allocatable, intent(in) :: sampleUvsV
        type (Star), dimension(:), allocatable, intent(in) :: sampleUvsW
        type (Star), dimension(:), allocatable, intent(in) :: sampleVvsW
        integer :: i
        integer :: unitUV, &
                   unitUW, &
                   unitVW

        if (allocated(sampleUvsV) .and. allocated(sampleUvsW) &
                                  .and. allocated(sampleVvsW)) then

            open(getNewUnit(unitUV), file = UV_PATH, status='old')
            do i = 1, size(sampleUvsV)
                write(unitUV, OUTPUT_FORMAT) sampleUvsV(i)%vel(1), &
                                             sampleUvsV(i)%vel(2)
            end do 

            open(getNewUnit(unitUW), file = UW_PATH, status='old')
            do i = 1, size(sampleUvsW)
                write(unitUW, OUTPUT_FORMAT) sampleUvsW(i)%vel(1), &
                                             sampleUvsW(i)%vel(2)
            end do 

            open(getNewUnit(unitVW), file = VW_PATH, status='old')
            do i = 1, size(sampleVvsW)
                write(unitVW, OUTPUT_FORMAT) sampleVvsW(i)%vel(1), &
                                             sampleVvsW(i)%vel(2)
            end do
        end if 
    end subroutine threeSamples_plotUVWvsUVW
end module plot_uvw
