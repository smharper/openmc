module geometry_global

  use, intrinsic :: ISO_C_BINDING

#ifdef MPIF08
  use mpi_f08
#endif

  !use bank_header,      only: Bank
  !use cmfd_header
  use constants
  use dict_header,      only: DictCharInt, DictIntInt
  use geometry_header,  only: Cell, Universe, Lattice, LatticeContainer
  !use material_header,  only: Material
  !use mesh_header,      only: RegularMesh
  !use mgxs_header,      only: Mgxs, MgxsContainer
  !use nuclide_header
  !use plot_header,      only: ObjectPlot
  !use sab_header,       only: SAlphaBeta
  !use set_header,       only: SetInt
  !use stl_vector,       only: VectorInt
  use surface_header,   only: SurfaceContainer
  !use source_header,    only: SourceDistribution
  !use tally_header,     only: TallyObject, TallyDerivative
  !use tally_filter_header, only: TallyFilterContainer, TallyFilterMatch
  !use trigger_header,   only: KTrigger
  !use timer_header,     only: Timer
  !use volume_header,    only: VolumeCalculation

  implicit none

  ! Main arrays
  type(Cell),             allocatable, target :: cells(:)
  type(Universe),         allocatable, target :: universes(:)
  type(LatticeContainer), allocatable, target :: lattices(:)
  type(SurfaceContainer), allocatable, target :: surfaces(:)

  ! Size of main arrays
  integer(C_INT32_T), bind(C) :: n_cells     ! # of cells
  integer :: n_universes ! # of universes
  integer :: n_lattices  ! # of lattices
  integer :: n_surfaces  ! # of surfaces

  ! These dictionaries provide a fast lookup mechanism -- the key is the
  ! user-specified identifier and the value is the index in the corresponding
  ! array
  type(DictIntInt) :: cell_dict
  type(DictIntInt) :: universe_dict
  type(DictIntInt) :: lattice_dict
  type(DictIntInt) :: surface_dict

  ! Number of lost particles
  integer :: n_lost_particles

contains

!===============================================================================
! FREE_MEMORY deallocates and clears  all global allocatable arrays in the
! program
!===============================================================================

  subroutine free_memory_geometry()

    integer :: i ! Loop Index

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)
    if (allocated(surfaces)) deallocate(surfaces)

    ! Deallocate dictionaries
    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()
    call surface_dict % clear()

  end subroutine free_memory_geometry

end module geometry_global
