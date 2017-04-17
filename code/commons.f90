module commons

    implicit none

    private

    public :: VW_PATH, &
              UW_PATH, &
              UV_PATH, &
              INPUT_PATH, &
              OUTPUT_FORMAT            
        character(len = *), parameter :: VW_PATH = './outputs/vw.dat'
        character(len = *), parameter :: UW_PATH = './outputs/uw.dat'
        character(len = *), parameter :: UV_PATH = './outputs/uv.dat'
        character(len = *), parameter :: INPUT_PATH = '/home/gemma&
        											  &/Documents/program&
        											  &/WD_population_40pc&
        											  &/output_data&
                                              		  &/boot_rowell_thin_1.out'
        character(len = *), parameter :: OUTPUT_FORMAT = '(2(f12.6,3x))'


contains

end module