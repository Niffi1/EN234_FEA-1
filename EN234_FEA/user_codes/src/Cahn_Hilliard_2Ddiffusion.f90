!     Subroutines for Cahn-Hilliard 2D diffusion problem
!     element identifier == 600

!==========================SUBROUTINE Cahn_Hilliard_2Ddiffusion ==============================
subroutine Cahn_Hilliard_2Ddiffusion(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Globals, only: TIME, DTIME                  ! Total analysis time and time increment
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2, (x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points, kint, ii

    real (prec)  ::  D(6,6), B(6,length_dof_array), qvec(6)
    real (prec)  ::  dxidx(2,2), determinant                    ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)                  ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  Dcoef, kappa, theta                        ! Material properties
    real (prec)  ::  mu, dmu, dmudx1, dmudx2, c, dc, dcdx1, dcdx2
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 2D linear elastic elements
    !     El props are:
    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio


    fail = .false.
    
    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 4
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 4) n_points = 4
    if (n_nodes == 8) n_points = 9

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
    Dcoef = element_properties(1)
    kappa = element_properties(2)
    theta = element_properties(3)
!    dof_total_new = dof_total + dof_increment
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        B(1,1:2*n_nodes-1:2) = N(1:n_nodes)
        B(2,2:2*n_nodes-0:2) = N(1:n_nodes)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(4,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(5,2:2*n_nodes-0:2) = dNdx(1:n_nodes,1)
        B(6,2:2*n_nodes-0:2) = dNdx(1:n_nodes,2)
        ! map mua, ca to mu, c at integration point
        mu = 0.D0
        dmu = 0.D0
        dmudx1 = 0.D0
        dmudx2 = 0.D0
        c = 0.D0
        dc = 0.D0
        dcdx1 = 0.D0
        dcdx2 = 0.D0
        do ii = 1, n_nodes
            mu = mu + N(ii)*dof_total(2*ii-1)
            dmu = dmu + N(ii)*dof_increment(2*ii-1)
            dmudx1 = dmudx1 + dNdx(ii,1)*(dof_total(2*ii-1)+theta*dof_increment(2*ii-1))
            dmudx2 = dmudx2 + dNdx(ii,2)*(dof_total(2*ii-1)+theta*dof_increment(2*ii-1))
            c = c + N(ii)*dof_total(2*ii)
            dc = dc + N(ii)*dof_increment(2*ii)
            dcdx1 = dcdx1 + dNdx(ii,1)*(dof_total(2*ii)+dof_increment(2*ii))
            dcdx2 = dcdx2 + dNdx(ii,2)*(dof_total(2*ii)+dof_increment(2*ii))
        end do
        ! find q vector
        qvec = 0.D0
        qvec = (/ mu+dmu-(c+dc)*((c+dc)**2-1.D0), dc/DTIME, -kappa*dcdx1, -kappa*dcdx2, Dcoef*dmudx1, Dcoef*dmudx2 /)
        ! find D matrix
        D = 0.d0
        D(1,1) = 1.D0
        D(1,2) = 1.D0-3*(c+dc)**2
        D(2,2) = 1.D0/DTIME
        D(3,5) = -kappa
        D(4,6) = -kappa
        D(5,3) = theta*Dcoef
        D(6,4) = theta*Dcoef

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),qvec)*w(kint)*determinant
        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
            + matmul(transpose(B(1:6,1:2*n_nodes)),matmul(D,B(1:6,1:2*n_nodes)))*w(kint)*determinant

    end do

    return
end subroutine Cahn_Hilliard_2Ddiffusion

