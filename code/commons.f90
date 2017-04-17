module commons

    implicit none

    private

    public :: VW_PATH, &
              UW_PATH, &
              UV_PATH            
        character(len = *), parameter :: VW_PATH = './outputs/vw.dat'
        character(len = *), parameter :: UW_PATH = './outputs/uw.dat'
        character(len = *), parameter :: UV_PATH = './outputs/uv.dat'

contains

end module