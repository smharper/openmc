module geometry_xml

!  use algorithm,        only: find
!  use cmfd_input,       only: configure_cmfd
  use constants
!  use dict_header,      only: DictIntInt, DictCharInt, ElemKeyValueCI
!  use distribution_multivariate
!  use distribution_univariate
!  use endf,             only: reaction_name
  use error,            only: fatal_error, warning
  use geometry_global
  use geometry_header,  only: Cell, Lattice, RectLattice, HexLattice, &
                              get_temperatures, root_universe
  use global
!  use hdf5_interface
!  use list_header,      only: ListChar, ListInt, ListReal
!  use mesh_header,      only: RegularMesh
!  use message_passing
!  use mgxs_data,        only: create_macro_xs, read_mgxs
!  use multipole,        only: multipole_read
  use output,           only: write_message, title, header
!  use plot_header
!  use random_lcg,       only: prn, seed
  use surface_header
!  use set_header,       only: SetChar
!  use stl_vector,       only: VectorInt, VectorReal, VectorChar
!  use string,           only: to_lower, to_str, str_to_int, str_to_real, &
!                              starts_with, ends_with, tokenize, split_string, &
!                              zero_padded
!  use tally_header,     only: TallyObject
!  use tally_filter_header, only: TallyFilterContainer
!  use tally_filter
!  use tally_initialize, only: add_tallies
  use xml_interface

  implicit none
  save

contains

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    integer :: i, j, k, m, i_x, i_a, input_index
    integer :: n, n_mats, n_x, n_y, n_z, n_rings, n_rlats, n_hlats
    integer :: univ_id
    integer :: n_cells_in_univ
    integer :: coeffs_reqd
    integer :: i_xmin, i_xmax, i_ymin, i_ymax, i_zmin, i_zmax
    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer, allocatable :: temp_int_array(:)
    real(8) :: phi, theta, psi
    real(8), allocatable :: coeffs(:)
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    character(:), allocatable :: region_spec
    type(Cell),     pointer :: c
    class(Surface), pointer :: s
    class(Lattice), pointer :: lat
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_cell
    type(XMLNode) :: node_surf
    type(XMLNode) :: node_lat
    type(XMLNode), allocatable :: node_cell_list(:)
    type(XMLNode), allocatable :: node_surf_list(:)
    type(XMLNode), allocatable :: node_rlat_list(:)
    type(XMLNode), allocatable :: node_hlat_list(:)
    type(VectorInt) :: tokens
    type(VectorInt) :: rpn
    type(VectorInt) :: fill_univ_ids ! List of fill universe IDs
    type(VectorInt) :: univ_ids      ! List of all universe IDs
    type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                           ! universe contains


    ! Display output message
    call write_message("Reading geometry XML file...", 5)

    ! ==========================================================================
    ! READ CELLS FROM GEOMETRY.XML

    ! Check if geometry.xml exists
    filename = trim(path_input) // "geometry.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Geometry XML file '" // trim(filename) // "' does not &
           &exist!")
    end if

    ! Parse geometry.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Get pointer to list of XML <cell>
    call get_node_list(root, "cell", node_cell_list)

    ! Get number of <cell> tags
    n_cells = size(node_cell_list)

    ! Check for no cells
    if (n_cells == 0) then
      call fatal_error("No cells found in geometry.xml!")
    end if

    ! Allocate cells array
    allocate(cells(n_cells))

    if (check_overlaps) then
      allocate(overlap_check_cnt(n_cells))
      overlap_check_cnt = 0
    end if

    n_universes = 0
    do i = 1, n_cells
      c => cells(i)

      ! Initialize distribcell instances and distribcell index
      c % instances = 0
      c % distribcell_index = NONE

      ! Get pointer to i-th cell node
      node_cell = node_cell_list(i)

      ! Copy data into cells
      if (check_for_node(node_cell, "id")) then
        call get_node_value(node_cell, "id", c % id)
      else
        call fatal_error("Must specify id of cell in geometry XML file.")
      end if

      ! Copy cell name
      if (check_for_node(node_cell, "name")) then
        call get_node_value(node_cell, "name", c % name)
      end if

      if (check_for_node(node_cell, "universe")) then
        call get_node_value(node_cell, "universe", c % universe)
      else
        c % universe = 0
      end if
      if (check_for_node(node_cell, "fill")) then
        call get_node_value(node_cell, "fill", c % fill)
        if (find(fill_univ_ids, c % fill) == -1) &
             call fill_univ_ids % push_back(c % fill)
      else
        c % fill = NONE
      end if

      ! Check to make sure 'id' hasn't been used
      if (cell_dict % has_key(c % id)) then
        call fatal_error("Two or more cells use the same unique ID: " &
             // to_str(c % id))
      end if

      ! Read material
      if (check_for_node(node_cell, "material")) then
        n_mats = node_word_count(node_cell, "material")

        if (n_mats > 0) then
          allocate(sarray(n_mats))
          call get_node_array(node_cell, "material", sarray)

          allocate(c % material(n_mats))
          do j = 1, n_mats
            select case(trim(to_lower(sarray(j))))
            case ('void')
              c % material(j) = MATERIAL_VOID
            case default
              c % material(j) = int(str_to_int(sarray(j)), 4)

              ! Check for error
              if (c % material(j) == ERROR_INT) then
                call fatal_error("Invalid material specified on cell " &
                     // to_str(c % id))
              end if
            end select
          end do

          deallocate(sarray)

        else
          allocate(c % material(1))
          c % material(1) = NONE
        end if

      else
        allocate(c % material(1))
        c % material(1) = NONE
      end if

      ! Check to make sure that either material or fill was specified
      if (c % material(1) == NONE .and. c % fill == NONE) then
        call fatal_error("Neither material nor fill was specified for cell " &
             // trim(to_str(c % id)))
      end if

      ! Check to make sure that both material and fill haven't been
      ! specified simultaneously
      if (c % material(1) /= NONE .and. c % fill /= NONE) then
        call fatal_error("Cannot specify material and fill simultaneously")
      end if

      ! Check for region specification (also under deprecated name surfaces)
      if (check_for_node(node_cell, "surfaces")) then
        call warning("The use of 'surfaces' is deprecated and will be &
             &disallowed in a future release.  Use 'region' instead. The &
             &openmc-update-inputs utility can be used to automatically &
             &update geometry.xml files.")
        region_spec = node_value_string(node_cell, "surfaces")
        call get_node_value(node_cell, "surfaces", region_spec)
      elseif (check_for_node(node_cell, "region")) then
        region_spec = node_value_string(node_cell, "region")
      else
        region_spec = ''
      end if

      if (len_trim(region_spec) > 0) then
        ! Create surfaces array from string
        call tokenize(region_spec, tokens)

        ! Use shunting-yard algorithm to determine RPN for surface algorithm
        call generate_rpn(c%id, tokens, rpn)

        ! Copy region spec and RPN form to cell arrays
        allocate(c % region(tokens%size()))
        allocate(c % rpn(rpn%size()))
        c % region(:) = tokens%data(1:tokens%size())
        c % rpn(:) = rpn%data(1:rpn%size())

        call tokens%clear()
        call rpn%clear()
      end if
      if (.not. allocated(c%region)) allocate(c%region(0))
      if (.not. allocated(c%rpn)) allocate(c%rpn(0))

      ! Check if this is a simple cell
      if (any(c%rpn == OP_COMPLEMENT) .or. any(c%rpn == OP_UNION)) then
        c%simple = .false.
      else
        c%simple = .true.
      end if

      ! Rotation matrix
      if (check_for_node(node_cell, "rotation")) then
        ! Rotations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          call fatal_error("Cannot apply a rotation to cell " // trim(to_str(&
               &c % id)) // " because it is not filled with another universe")
        end if

        ! Read number of rotation parameters
        n = node_word_count(node_cell, "rotation")
        if (n /= 3) then
          call fatal_error("Incorrect number of rotation parameters on cell " &
               // to_str(c % id))
        end if

        ! Copy rotation angles in x,y,z directions
        allocate(c % rotation(3))
        call get_node_array(node_cell, "rotation", c % rotation)
        phi   = -c % rotation(1) * PI/180.0_8
        theta = -c % rotation(2) * PI/180.0_8
        psi   = -c % rotation(3) * PI/180.0_8

        ! Calculate rotation matrix based on angles given
        allocate(c % rotation_matrix(3,3))
        c % rotation_matrix = reshape((/ &
             cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta), &
             -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), &
             cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), &
             sin(phi)*cos(theta), &
             sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi), &
             -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi), &
             cos(phi)*cos(theta) /), (/ 3,3 /))
      end if

      ! Translation vector
      if (check_for_node(node_cell, "translation")) then
        ! Translations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          call fatal_error("Cannot apply a translation to cell " &
               // trim(to_str(c % id)) // " because it is not filled with &
               &another universe")
        end if

        ! Read number of translation parameters
        n = node_word_count(node_cell, "translation")
        if (n /= 3) then
          call fatal_error("Incorrect number of translation parameters on &
               &cell " // to_str(c % id))
        end if

        ! Copy translation vector
        allocate(c % translation(3))
        call get_node_array(node_cell, "translation", c % translation)
      end if

      ! Read cell temperatures.  If the temperature is not specified, set it to
      ! ERROR_REAL for now.  During initialization we'll replace ERROR_REAL with
      ! the temperature from the material data.
      if (check_for_node(node_cell, "temperature")) then
        n = node_word_count(node_cell, "temperature")
        if (n > 0) then
          ! Make sure this is a "normal" cell.
          if (c % material(1) == NONE) call fatal_error("Cell " &
               // trim(to_str(c % id)) // " was specified with a temperature &
               &but no material. Temperature specification is only valid for &
               &cells filled with a material.")

          ! Copy in temperatures
          allocate(c % sqrtkT(n))
          call get_node_array(node_cell, "temperature", c % sqrtkT)

          ! Make sure all temperatues are positive
          do j = 1, size(c % sqrtkT)
            if (c % sqrtkT(j) < ZERO) call fatal_error("Cell " &
                 // trim(to_str(c % id)) // " was specified with a negative &
                 &temperature. All cell temperatures must be non-negative.")
          end do

          ! Convert to sqrt(kT)
          c % sqrtkT(:) = sqrt(K_BOLTZMANN * c % sqrtkT(:))
        else
          allocate(c % sqrtkT(1))
          c % sqrtkT(1) = ERROR_REAL
        end if
      else
        allocate(c % sqrtkT(1))
        c % sqrtkT = ERROR_REAL
      end if

      ! Add cell to dictionary
      call cell_dict % add_key(c % id, i)

      ! For cells, we also need to check if there's a new universe --
      ! also for every cell add 1 to the count of cells for the
      ! specified universe
      univ_id = c % universe
      if (.not. cells_in_univ_dict % has_key(univ_id)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call universe_dict % add_key(univ_id, n_universes)
        call univ_ids % push_back(univ_id)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get_key(univ_id)
      end if
      call cells_in_univ_dict % add_key(univ_id, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ SURFACES FROM GEOMETRY.XML

    ! This variable is used to check whether at least one boundary condition was
    ! applied to a surface
    boundary_exists = .false.

    ! get pointer to list of xml <surface>
    call get_node_list(root, "surface", node_surf_list)

    ! Get number of <surface> tags
    n_surfaces = size(node_surf_list)

    ! Check for no surfaces
    if (n_surfaces == 0) then
      call fatal_error("No surfaces found in geometry.xml!")
    end if

    xmin = INFINITY
    xmax = -INFINITY
    ymin = INFINITY
    ymax = -INFINITY
    zmin = INFINITY
    zmax = -INFINITY

    ! Allocate cells array
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
      ! Get pointer to i-th surface node
      node_surf = node_surf_list(i)

      ! Copy and interpret surface type
      word = ''
      if (check_for_node(node_surf, "type")) &
           call get_node_value(node_surf, "type", word)
      select case(to_lower(word))
      case ('x-plane')
        coeffs_reqd  = 1
        allocate(SurfaceXPlane :: surfaces(i)%obj)
      case ('y-plane')
        coeffs_reqd  = 1
        allocate(SurfaceYPlane :: surfaces(i)%obj)
      case ('z-plane')
        coeffs_reqd  = 1
        allocate(SurfaceZPlane :: surfaces(i)%obj)
      case ('plane')
        coeffs_reqd  = 4
        allocate(SurfacePlane :: surfaces(i)%obj)
      case ('x-cylinder')
        coeffs_reqd  = 3
        allocate(SurfaceXCylinder :: surfaces(i)%obj)
      case ('y-cylinder')
        coeffs_reqd  = 3
        allocate(SurfaceYCylinder :: surfaces(i)%obj)
      case ('z-cylinder')
        coeffs_reqd  = 3
        allocate(SurfaceZCylinder :: surfaces(i)%obj)
      case ('sphere')
        coeffs_reqd  = 4
        allocate(SurfaceSphere :: surfaces(i)%obj)
      case ('x-cone')
        coeffs_reqd  = 4
        allocate(SurfaceXCone :: surfaces(i)%obj)
      case ('y-cone')
        coeffs_reqd  = 4
        allocate(SurfaceYCone :: surfaces(i)%obj)
      case ('z-cone')
        coeffs_reqd  = 4
        allocate(SurfaceZCone :: surfaces(i)%obj)
      case ('quadric')
        coeffs_reqd  = 10
        allocate(SurfaceQuadric :: surfaces(i)%obj)
      case default
        call fatal_error("Invalid surface type: " // trim(word))
      end select

      s => surfaces(i)%obj

      ! Copy data into cells
      if (check_for_node(node_surf, "id")) then
        call get_node_value(node_surf, "id", s%id)
      else
        call fatal_error("Must specify id of surface in geometry XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (surface_dict % has_key(s%id)) then
        call fatal_error("Two or more surfaces use the same unique ID: " &
             // to_str(s%id))
      end if

      ! Copy surface name
      if (check_for_node(node_surf, "name")) then
        call get_node_value(node_surf, "name", s%name)
      end if

      ! Check to make sure that the proper number of coefficients
      ! have been specified for the given type of surface. Then copy
      ! surface coordinates.

      n = node_word_count(node_surf, "coeffs")
      if (n < coeffs_reqd) then
        call fatal_error("Not enough coefficients specified for surface: " &
             // trim(to_str(s%id)))
      elseif (n > coeffs_reqd) then
        call fatal_error("Too many coefficients specified for surface: " &
             // trim(to_str(s%id)))
      end if

      allocate(coeffs(n))
      call get_node_array(node_surf, "coeffs", coeffs)

      select type(s)
      type is (SurfaceXPlane)
        s%x0 = coeffs(1)

        ! Determine outer surfaces
        xmin = min(xmin, s % x0)
        xmax = max(xmax, s % x0)
        if (xmin == s % x0) i_xmin = i
        if (xmax == s % x0) i_xmax = i
      type is (SurfaceYPlane)
        s%y0 = coeffs(1)

        ! Determine outer surfaces
        ymin = min(ymin, s % y0)
        ymax = max(ymax, s % y0)
        if (ymin == s % y0) i_ymin = i
        if (ymax == s % y0) i_ymax = i
      type is (SurfaceZPlane)
        s%z0 = coeffs(1)

        ! Determine outer surfaces
        zmin = min(zmin, s % z0)
        zmax = max(zmax, s % z0)
        if (zmin == s % z0) i_zmin = i
        if (zmax == s % z0) i_zmax = i
      type is (SurfacePlane)
        s%A = coeffs(1)
        s%B = coeffs(2)
        s%C = coeffs(3)
        s%D = coeffs(4)
      type is (SurfaceXCylinder)
        s%y0 = coeffs(1)
        s%z0 = coeffs(2)
        s%r = coeffs(3)
      type is (SurfaceYCylinder)
        s%x0 = coeffs(1)
        s%z0 = coeffs(2)
        s%r = coeffs(3)
      type is (SurfaceZCylinder)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%r = coeffs(3)
      type is (SurfaceSphere)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r = coeffs(4)
      type is (SurfaceXCone)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r2 = coeffs(4)
      type is (SurfaceYCone)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r2 = coeffs(4)
      type is (SurfaceZCone)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r2 = coeffs(4)
      type is (SurfaceQuadric)
        s%A = coeffs(1)
        s%B = coeffs(2)
        s%C = coeffs(3)
        s%D = coeffs(4)
        s%E = coeffs(5)
        s%F = coeffs(6)
        s%G = coeffs(7)
        s%H = coeffs(8)
        s%J = coeffs(9)
        s%K = coeffs(10)
      end select

      ! No longer need coefficients
      deallocate(coeffs)

      ! Boundary conditions
      word = ''
      if (check_for_node(node_surf, "boundary")) &
           call get_node_value(node_surf, "boundary", word)
      select case (to_lower(word))
      case ('transmission', 'transmit', '')
        s%bc = BC_TRANSMIT
      case ('vacuum')
        s%bc = BC_VACUUM
        boundary_exists = .true.
      case ('reflective', 'reflect', 'reflecting')
        s%bc = BC_REFLECT
        boundary_exists = .true.
      case ('periodic')
        s%bc = BC_PERIODIC
        boundary_exists = .true.

        ! Check for specification of periodic surface
        if (check_for_node(node_surf, "periodic_surface_id")) then
          call get_node_value(node_surf, "periodic_surface_id", &
               s % i_periodic)
        end if
      case default
        call fatal_error("Unknown boundary condition '" // trim(word) // &
             &"' specified on surface " // trim(to_str(s%id)))
      end select
      ! Add surface to dictionary
      call surface_dict % add_key(s%id, i)
    end do

    ! Check to make sure a boundary condition was applied to at least one
    ! surface
    if (run_mode /= MODE_PLOTTING) then
      if (.not. boundary_exists) then
        call fatal_error("No boundary conditions were applied to any surfaces!")
      end if
    end if

    ! Determine opposite side for periodic boundaries
    do i = 1, size(surfaces)
      if (surfaces(i) % obj % bc == BC_PERIODIC) then
        select type (surf => surfaces(i) % obj)
        type is (SurfaceXPlane)
          if (surf % i_periodic == NONE) then
            if (i == i_xmin) then
              surf % i_periodic = i_xmax
            elseif (i == i_xmax) then
              surf % i_periodic = i_xmin
            else
              call fatal_error("Periodic boundary condition applied to &
                   &interior surface.")
            end if
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        type is (SurfaceYPlane)
          if (surf % i_periodic == NONE) then
            if (i == i_ymin) then
              surf % i_periodic = i_ymax
            elseif (i == i_ymax) then
              surf % i_periodic = i_ymin
            else
              call fatal_error("Periodic boundary condition applied to &
                   &interior surface.")
            end if
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        type is (SurfaceZPlane)
          if (surf % i_periodic == NONE) then
            if (i == i_zmin) then
              surf % i_periodic = i_zmax
            elseif (i == i_zmax) then
              surf % i_periodic = i_zmin
            else
              call fatal_error("Periodic boundary condition applied to &
                   &interior surface.")
            end if
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        type is (SurfacePlane)
          if (surf % i_periodic == NONE) then
            call fatal_error("No matching periodic surface specified for &
                 &periodic boundary condition on surface " // &
                 trim(to_str(surf % id)) // ".")
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        class default
          call fatal_error("Periodic boundary condition applied to &
               &non-planar surface.")
        end select

        ! Make sure opposite surface is also periodic
        associate (surf => surfaces(i) % obj)
          if (surfaces(surf % i_periodic) % obj % bc /= BC_PERIODIC) then
            call fatal_error("Could not find matching surface for periodic &
                 &boundary on surface " // trim(to_str(surf % id)) // ".")
          end if
        end associate
      end if
    end do

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    ! Get pointer to list of XML <lattice>
    call get_node_list(root, "lattice", node_rlat_list)
    call get_node_list(root, "hex_lattice", node_hlat_list)

    ! Allocate lattices array
    n_rlats = size(node_rlat_list)
    n_hlats = size(node_hlat_list)
    n_lattices = n_rlats + n_hlats
    allocate(lattices(n_lattices))

    RECT_LATTICES: do i = 1, n_rlats
      allocate(RectLattice::lattices(i) % obj)
      lat => lattices(i) % obj
      select type(lat)
      type is (RectLattice)

      ! Get pointer to i-th lattice
      node_lat = node_rlat_list(i)

      ! ID of lattice
      if (check_for_node(node_lat, "id")) then
        call get_node_value(node_lat, "id", lat % id)
      else
        call fatal_error("Must specify id of lattice in geometry XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        call fatal_error("Two or more lattices use the same unique ID: " &
             // to_str(lat % id))
      end if

      ! Copy lattice name
      if (check_for_node(node_lat, "name")) then
        call get_node_value(node_lat, "name", lat % name)
      end if

      ! Read number of lattice cells in each dimension
      n = node_word_count(node_lat, "dimension")
      if (n == 2) then
        call get_node_array(node_lat, "dimension", lat % n_cells(1:2))
        lat % n_cells(3) = 1
        lat % is_3d = .false.
      else if (n == 3) then
        call get_node_array(node_lat, "dimension", lat % n_cells)
        lat % is_3d = .true.
      else
        call fatal_error("Rectangular lattice must be two or three dimensions.")
      end if

      ! Read lattice lower-left location
      if (node_word_count(node_lat, "lower_left") /= n) then
        call fatal_error("Number of entries on <lower_left> must be the same &
             &as the number of entries on <dimension>.")
      end if

      allocate(lat % lower_left(n))
      call get_node_array(node_lat, "lower_left", lat % lower_left)

      ! Read lattice pitches.
      ! TODO: Remove this deprecation warning in a future release.
      if (check_for_node(node_lat, "width")) then
        call warning("The use of 'width' is deprecated and will be disallowed &
             &in a future release.  Use 'pitch' instead.  The utility openmc/&
             &src/utils/update_inputs.py can be used to automatically update &
             &geometry.xml files.")
        if (node_word_count(node_lat, "width") /= n) then
          call fatal_error("Number of entries on <pitch> must be the same as &
               &the number of entries on <dimension>.")
        end if

      else if (node_word_count(node_lat, "pitch") /= n) then
        call fatal_error("Number of entries on <pitch> must be the same as &
             &the number of entries on <dimension>.")
      end if

      allocate(lat % pitch(n))
      ! TODO: Remove the 'width' code in a future release.
      if (check_for_node(node_lat, "width")) then
        call get_node_array(node_lat, "width", lat % pitch)
      else
        call get_node_array(node_lat, "pitch", lat % pitch)
      end if

      ! TODO: Remove deprecation warning in a future release.
      if (check_for_node(node_lat, "type")) then
        call warning("The use of 'type' is no longer needed.  The utility &
             &openmc/src/utils/update_inputs.py can be used to automatically &
             &update geometry.xml files.")
      end if

      ! Copy number of dimensions
      n_x = lat % n_cells(1)
      n_y = lat % n_cells(2)
      n_z = lat % n_cells(3)
      allocate(lat % universes(n_x, n_y, n_z))

      ! Check that number of universes matches size
      n = node_word_count(node_lat, "universes")
      if (n /= n_x*n_y*n_z) then
        call fatal_error("Number of universes on <universes> does not match &
             &size of lattice " // trim(to_str(lat % id)) // ".")
      end if

      allocate(temp_int_array(n))
      call get_node_array(node_lat, "universes", temp_int_array)

      ! Read universes
      do m = 1, n_z
        do k = 0, n_y - 1
          do j = 1, n_x
            lat % universes(j, n_y - k, m) = &
                 temp_int_array(j + n_x*k + n_x*n_y*(m-1))
            if (find(fill_univ_ids, lat % universes(j, n_y - k, m)) == -1) &
                 call fill_univ_ids % push_back(lat % universes(j, n_y - k, m))
          end do
        end do
      end do
      deallocate(temp_int_array)

      ! Read outer universe for area outside lattice.
      lat % outer = NO_OUTER_UNIVERSE
      if (check_for_node(node_lat, "outer")) then
        call get_node_value(node_lat, "outer", lat % outer)
        if (find(fill_univ_ids, lat % outer) == -1) &
             call fill_univ_ids % push_back(lat % outer)
      end if

      ! Check for 'outside' nodes which are no longer supported.
      if (check_for_node(node_lat, "outside")) then
        call fatal_error("The use of 'outside' in lattices is no longer &
             &supported.  Instead, use 'outer' which defines a universe rather &
             &than a material.  The utility openmc/src/utils/update_inputs.py &
             &can be used automatically replace 'outside' with 'outer'.")
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, i)

      end select
    end do RECT_LATTICES

    HEX_LATTICES: do i = 1, n_hlats
      allocate(HexLattice::lattices(n_rlats + i) % obj)
      lat => lattices(n_rlats + i) % obj
      select type (lat)
      type is (HexLattice)

      ! Get pointer to i-th lattice
      node_lat = node_hlat_list(i)

      ! ID of lattice
      if (check_for_node(node_lat, "id")) then
        call get_node_value(node_lat, "id", lat % id)
      else
        call fatal_error("Must specify id of lattice in geometry XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        call fatal_error("Two or more lattices use the same unique ID: " &
             // to_str(lat % id))
      end if

      ! Copy lattice name
      if (check_for_node(node_lat, "name")) then
        call get_node_value(node_lat, "name", lat % name)
      end if

      ! Read number of lattice cells in each dimension
      call get_node_value(node_lat, "n_rings", lat % n_rings)
      if (check_for_node(node_lat, "n_axial")) then
        call get_node_value(node_lat, "n_axial", lat % n_axial)
        lat % is_3d = .true.
      else
        lat % n_axial = 1
        lat % is_3d = .false.
      end if

      ! Read lattice lower-left location
      n = node_word_count(node_lat, "center")
      if (lat % is_3d .and. n /= 3) then
        call fatal_error("A hexagonal lattice with <n_axial> must have &
             &<center> specified by 3 numbers.")
      else if ((.not. lat % is_3d) .and. n /= 2) then
        call fatal_error("A hexagonal lattice without <n_axial> must have &
             &<center> specified by 2 numbers.")
      end if

      allocate(lat % center(n))
      call get_node_array(node_lat, "center", lat % center)

      ! Read lattice pitches
      n = node_word_count(node_lat, "pitch")
      if (lat % is_3d .and. n /= 2) then
        call fatal_error("A hexagonal lattice with <n_axial> must have <pitch> &
              &specified by 2 numbers.")
      else if ((.not. lat % is_3d) .and. n /= 1) then
        call fatal_error("A hexagonal lattice without <n_axial> must have &
             &<pitch> specified by 1 number.")
      end if

      allocate(lat % pitch(n))
      call get_node_array(node_lat, "pitch", lat % pitch)

      ! Copy number of dimensions
      n_rings = lat % n_rings
      n_z = lat % n_axial
      allocate(lat % universes(2*n_rings - 1, 2*n_rings - 1, n_z))

      ! Check that number of universes matches size
      n = node_word_count(node_lat, "universes")
      if (n /= (3*n_rings**2 - 3*n_rings + 1)*n_z) then
        call fatal_error("Number of universes on <universes> does not match &
             &size of lattice " // trim(to_str(lat % id)) // ".")
      end if

      allocate(temp_int_array(n))
      call get_node_array(node_lat, "universes", temp_int_array)

      ! Read universes
      ! Universes in hexagonal lattices are stored in a manner that represents
      ! a skewed coordinate system: (x, alpha) rather than (x, y).  There is
      ! no obvious, direct relationship between the order of universes in the
      ! input and the order that they will be stored in the skewed array so
      ! the following code walks a set of index values across the skewed array
      ! in a manner that matches the input order.  Note that i_x = 0, i_a = 0
      ! corresponds to the center of the hexagonal lattice.

      input_index = 1
      do m = 1, n_z
        ! Initialize lattice indecies.
        i_x = 1
        i_a = n_rings - 1

        ! Map upper triangular region of hexagonal lattice.
        do k = 1, n_rings-1
          ! Walk index to lower-left neighbor of last row start.
          i_x = i_x - 1
          do j = 1, k
            ! Place universe in array.
            lat % universes(i_x + n_rings, i_a + n_rings, m) = &
                 temp_int_array(input_index)
            if (find(fill_univ_ids, temp_int_array(input_index)) == -1) &
                 call fill_univ_ids % push_back(temp_int_array(input_index))
            ! Walk index to closest non-adjacent right neighbor.
            i_x = i_x + 2
            i_a = i_a - 1
            ! Increment XML array index.
            input_index = input_index + 1
          end do
          ! Return lattice index to start of current row.
          i_x = i_x - 2*k
          i_a = i_a + k
        end do

        ! Map middle square region of hexagonal lattice.
        do k = 1, 2*n_rings - 1
          if (mod(k, 2) == 1) then
            ! Walk index to lower-left neighbor of last row start.
            i_x = i_x - 1
          else
            ! Walk index to lower-right neighbor of last row start
            i_x = i_x + 1
            i_a = i_a - 1
          end if
          do j = 1, n_rings - mod(k-1, 2)
            ! Place universe in array.
            lat % universes(i_x + n_rings, i_a + n_rings, m) = &
                 temp_int_array(input_index)
            if (find(fill_univ_ids, temp_int_array(input_index)) == -1) &
                 call fill_univ_ids % push_back(temp_int_array(input_index))
            ! Walk index to closest non-adjacent right neighbor.
            i_x = i_x + 2
            i_a = i_a - 1
            ! Increment XML array index.
            input_index = input_index + 1
          end do
          ! Return lattice index to start of current row.
          i_x = i_x - 2*(n_rings - mod(k-1, 2))
          i_a = i_a + n_rings - mod(k-1, 2)
        end do

        ! Map lower triangular region of hexagonal lattice.
        do k = 1, n_rings-1
          ! Walk index to lower-right neighbor of last row start.
          i_x = i_x + 1
          i_a = i_a - 1
          do j = 1, n_rings - k
            ! Place universe in array.
            lat % universes(i_x + n_rings, i_a + n_rings, m) = &
                 temp_int_array(input_index)
            if (find(fill_univ_ids, temp_int_array(input_index)) == -1) &
                 call fill_univ_ids % push_back(temp_int_array(input_index))
            ! Walk index to closest non-adjacent right neighbor.
            i_x = i_x + 2
            i_a = i_a - 1
            ! Increment XML array index.
            input_index = input_index + 1
          end do
          ! Return lattice index to start of current row.
          i_x = i_x - 2*(n_rings - k)
          i_a = i_a + n_rings - k
        end do
      end do
      deallocate(temp_int_array)

      ! Read outer universe for area outside lattice.
      lat % outer = NO_OUTER_UNIVERSE
      if (check_for_node(node_lat, "outer")) then
        call get_node_value(node_lat, "outer", lat % outer)
        if (find(fill_univ_ids, lat % outer) == -1) &
             call fill_univ_ids % push_back(lat % outer)
      end if

      ! Check for 'outside' nodes which are no longer supported.
      if (check_for_node(node_lat, "outside")) then
        call fatal_error("The use of 'outside' in lattices is no longer &
             &supported.  Instead, use 'outer' which defines a universe rather &
             &than a material.  The utility openmc/src/utils/update_inputs.py &
             &can be used automatically replace 'outside' with 'outer'.")
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, n_rlats + i)

      end select
    end do HEX_LATTICES

    ! ==========================================================================
    ! SETUP UNIVERSES

    ! Allocate universes, universe cell arrays, and assign base universe
    allocate(universes(n_universes))
    do i = 1, n_universes
      associate (u => universes(i))
        u % id = univ_ids % data(i)

        ! Allocate cell list
        n_cells_in_univ = cells_in_univ_dict % get_key(u % id)
        allocate(u % cells(n_cells_in_univ))
        u % cells(:) = 0

        ! Check whether universe is a fill universe
        if (find(fill_univ_ids, u % id) == -1) then
          if (root_universe > 0) then
            call fatal_error("Two or more universes are not used as fill &
                 &universes, so it is not possible to distinguish which one &
                 &is the root universe.")
          else
            root_universe = i
          end if
        end if
      end associate
    end do

    do i = 1, n_cells
      ! Get index in universes array
      j = universe_dict % get_key(cells(i) % universe)

      ! Set the first zero entry in the universe cells array to the index in the
      ! global cells array
      associate (u => universes(j))
        u % cells(find(u % cells, 0)) = i
      end associate
    end do

    ! Clear dictionary
    call cells_in_univ_dict%clear()

    ! Close geometry XML file
    call doc % clear()

  end subroutine read_geometry_xml

!===============================================================================
! GENERATE_RPN implements the shunting-yard algorithm to generate a Reverse
! Polish notation (RPN) expression for the region specification of a cell given
! the infix notation.
!===============================================================================

  subroutine generate_rpn(cell_id, tokens, output)
    integer, intent(in) :: cell_id
    type(VectorInt), intent(in) :: tokens    ! infix notation
    type(VectorInt), intent(inout) :: output ! RPN notation

    integer :: i
    integer :: token
    integer :: op
    type(VectorInt) :: stack

    do i = 1, tokens%size()
      token = tokens%data(i)

      if (token < OP_UNION) then
        ! If token is not an operator, add it to output
        call output%push_back(token)

      elseif (token < OP_RIGHT_PAREN) then
        ! Regular operators union, intersection, complement
        do while (stack%size() > 0)
          op = stack%data(stack%size())

          if (op < OP_RIGHT_PAREN .and. &
               ((token == OP_COMPLEMENT .and. token < op) .or. &
               (token /= OP_COMPLEMENT .and. token <= op))) then
            ! While there is an operator, op, on top of the stack, if the token
            ! is left-associative and its precedence is less than or equal to
            ! that of op or if the token is right-associative and its precedence
            ! is less than that of op, move op to the output queue and push the
            ! token on to the stack. Note that only complement is
            ! right-associative.
            call output%push_back(op)
            call stack%pop_back()
          else
            exit
          end if
        end do

        call stack%push_back(token)

      elseif (token == OP_LEFT_PAREN) then
        ! If the token is a left parenthesis, push it onto the stack
        call stack%push_back(token)

      else
        ! If the token is a right parenthesis, move operators from the stack to
        ! the output queue until reaching the left parenthesis.
        do
          ! If we run out of operators without finding a left parenthesis, it
          ! means there are mismatched parentheses.
          if (stack%size() == 0) then
            call fatal_error('Mimatched parentheses in region specification &
                 &for cell ' // trim(to_str(cell_id)) // '.')
          end if

          op = stack%data(stack%size())
          if (op == OP_LEFT_PAREN) exit
          call output%push_back(op)
          call stack%pop_back()
        end do

        ! Pop the left parenthesis.
        call stack%pop_back()
      end if
    end do

    ! While there are operators on the stack, move them to the output queue
    do while (stack%size() > 0)
      op = stack%data(stack%size())

      ! If the operator is a parenthesis, it is mismatched
      if (op >= OP_RIGHT_PAREN) then
        call fatal_error('Mimatched parentheses in region specification &
             &for cell ' // trim(to_str(cell_id)) // '.')
      end if

      call output%push_back(op)
      call stack%pop_back()
    end do
  end subroutine generate_rpn

end module geometry_xml
