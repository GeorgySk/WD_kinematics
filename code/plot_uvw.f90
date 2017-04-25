module plot_uvw
    use, intrinsic :: iso_fortran_env, dp=>real64
    use derived_types, only: Star
    use commons, only: OUTPUT_FORMAT, &
                       UVW_FORMAT
    use files, only: getNewUnit

    implicit none
    character(len = *), parameter :: VW_PATH = './outputs/observ&
                                                &/limoges/vw.dat'
    character(len = *), parameter :: UW_PATH = './outputs/observ&
                                                &/limoges/uw.dat'
    character(len = *), parameter :: UV_PATH = './outputs/observ&
                                                &/limoges/uv.dat'
    ! Following samples are used when plots use different data
    type (Star), dimension(:), allocatable, save :: sampleUvsV, &
                                                    sampleUvsW, &
                                                    sampleVvsW

    private :: UV_PATH, &
               UW_PATH, &
               VW_PATH

    public :: sampleUvsV, &
              sampleUvsW, &
              sampleVvsW, &
              plotUVWvsUVW

contains


    subroutine plotUVWvsUVW()
        integer :: i
        real(dp) :: vel_sum = 0.0_dp
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
    end subroutine plotUVWvsUVW


    subroutine fillFullDataForUVWPlot(vel_hel, &
                                      counter, &
                                      vel_sum, &
                                      vel_array, &
                                      unitUVW)
        
        real*8, dimension(3), intent(in) :: vel_hel
        integer, intent(in) :: counter
        real*8, dimension(3), intent(inout) :: vel_sum
        real*8, dimension(:, :), intent(inout) :: vel_array
        ! TODO: hide this somewhere
        integer, intent(in) :: unitUVW

        vel_sum = vel_sum + vel_hel
        vel_array(:, counter) = vel_hel

        write(unitUVW, UVW_FORMAT) vel_hel
        
    end subroutine fillFullDataForUVWPlot


end module plot_uvw
