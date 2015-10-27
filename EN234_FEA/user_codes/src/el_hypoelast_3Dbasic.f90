!     Subroutines for basic 3D hypo-elastic elements

!==========================SUBROUTINE el_hypoelast_3dbasic ==============================
subroutine el_hypoelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
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

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    ! Local Variables
    integer      :: n_points,kint,ii

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! the material tangent stiffness matrix
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  s0, e0, n0,K0                      ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    ! for B-bar element
    real (prec)  ::  el_vol          ! element volume
    real (prec)  ::  B_bar(6,length_dof_array), tempmatrix(6,length_dof_array)

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    s0 = element_properties(1)
    e0 = element_properties(2)
    n0 = element_properties(3)
    K0 = element_properties(4)

    D = 0.D0
    el_vol = 0.D0
    dNbardx = 0.D0
    ! hypoelasticity then use B-bar element with identifier 1004
    !     --  Loop over integration points
    ! compute average quantities: el_vol, dNbardx
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3)*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
    end do

    dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3)/el_vol
    ! modify the B matrix and strain
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        tempmatrix = 0.D0
        B_bar = 0.D0

        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        do ii = 1, n_nodes
            tempmatrix(1:3,3*ii-2) = dNbardx(ii,1) - dNdx(ii,1)
            tempmatrix(1:3,3*ii-1) = dNbardx(ii,2) - dNdx(ii,2)
            tempmatrix(1:3,3*ii-0) = dNbardx(ii,3) - dNdx(ii,3)
        end do

        B_bar = B + 1.D0/3.D0*tempmatrix

        strain = matmul(B_bar,dof_total)
        dstrain = matmul(B_bar,dof_increment)
        strain = strain + dstrain
        call mat_stiffness_3d(s0, e0, n0, K0, strain, stress, D)
!        stress = matmul(D,strain)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_bar),stress)*w(kint)*determinant
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B_bar(1:6,1:3*n_nodes)),matmul(D,B_bar(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do

    return
end subroutine el_hypoelast_3dbasic


!==========================SUBROUTINE fieldvars_hypoelast_3dbasic ==============================
subroutine fieldvars_hypoelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only: dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
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
    integer, intent( in )         :: n_field_variables                                      ! # field variables

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

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,ii

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  s0, e0, n0,K0                     ! Material properties
    real (prec)  ::  p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    ! for B-bar element
    real (prec)  ::  el_vol           ! element volume
    real (prec)  ::  B_bar(6,length_dof_array), tempmatrix(6,length_dof_array)

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    s0 = element_properties(1)
    e0 = element_properties(2)
    n0 = element_properties(3)
    K0 = element_properties(4)

    D = 0.D0
    el_vol = 0.D0
    dNbardx = 0.D0
    ! Hypoelasticity then use B-bar element with identifier 1004
    !     --  Loop over integration points
    ! compute average quantities: el_vol, dNbardx
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3)*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
    end do

    dNbardx(1:n_nodes,1:3) = dNbardx(1:n_nodes,1:3)/el_vol

    ! modify the B matrix and strain
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        tempmatrix = 0.D0
        B_bar = 0.D0

        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        do ii = 1, n_nodes
            tempmatrix(1:3,3*ii-2) = dNbardx(ii,1) - dNdx(ii,1)
            tempmatrix(1:3,3*ii-1) = dNbardx(ii,2) - dNdx(ii,2)
            tempmatrix(1:3,3*ii-0) = dNbardx(ii,3) - dNdx(ii,3)
        end do

        B_bar = B + 1.D0/3.D0*tempmatrix

        strain = matmul(B_bar,dof_total)
        dstrain = matmul(B_bar,dof_increment)
        strain = strain + dstrain
        call mat_stiffness_3d(s0, e0, n0, K0, strain, stress, D)
!        stress = matmul(D,strain)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1, n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'e11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'e22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'e33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'e12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'e13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'e23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(6)*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
    end do

    return
end subroutine fieldvars_hypoelast_3dbasic


!==========================SUBROUTINE mat_stiffness_3d ==============================
subroutine mat_stiffness_3d(s0, e0, n0, K0, strain, stress, D)
    real (kind=8) ::  s0, e0, n0,K0                      ! Material properties
    real (kind=8) ::  strain(6), dev_strain(6)           ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (kind=8) ::  stress(6)                          ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (kind=8) ::  D(6,6)                             ! the material tangent stiffness matrix
    real (kind=8) ::  e_dyadic_e(6,6), mat1(6,6), mat2(6,6)
    real (kind=8) ::  se, ee, dsde, e_kk
    integer ::  i, j

    ! variable initiation
    ee = 0.D0
    se = 0.D0
    dsde = 0.D0
    e_kk = 0.D0
    ! find deviation strain
    e_kk = strain(1) + strain(2) + strain(3)
    do i = 1, 3
        dev_strain(i) = strain(i) - e_kk/3.D0
        dev_strain(i+3) = strain(i+3)/2.D0
    end do
    ! find equiv. strain
    ee = dev_strain(1)**2 + dev_strain(2)**2 + dev_strain(3)**2
    ee = ee + 2.D0*( dev_strain(4)**2 + dev_strain(5)**2 + dev_strain(6)**2 )
    ee = dsqrt( 2.D0*ee/3.D0 )
    ! find equiv. stress
    if (ee .le. e0) then
        se = (1.D0+n0**2)/(n0-1.D0)**2 - ( n0/(n0-1.D0) - ee/e0 )**2
        se = s0*( dsqrt(se) - 1.D0/(n0-1) )
    else
        se = s0*( (ee/e0)**(1.D0/n0) )
    end if
    ! find uniaxial tangent
    if (ee .le. e0) then
        if (n0-1 .le. 1.0D-12) then
            dsde = s0/e0
        else
            dsde = s0*( n0/(n0-1.D0) -ee/e0 )/e0
            dsde = dsde/dsqrt( (1.D0+n0**2)/(n0-1.D0)**2 - ( -n0/(n0-1.D0) - ee/e0 )**2 )
        end if
    else
        dsde = s0*( (ee/e0)**(1.D0/n0) )/(n0*ee)
    end if
    ! find the tangent stiffness
    D = 0.D0
    e_dyadic_e = 0.D0
    mat1 = 0.D0
    mat2 = 0.D0
    do i = 1, 3
        mat1(i,i) = 2
        mat1(i+3,i+3) = 1
        do j = 1, 3
            mat2(i,j) = 1
        end do
    enddo
    e_dyadic_e = spread(dev_strain,dim=2,ncopies=6)*spread(dev_strain,dim=1,ncopies=6)
    if (ee .ge. 0) then
        D = se/ee/3.D0*mat1 + (K0-2.D0*se/ee/9.D0)*mat2 + 4*(dsde-se/ee)/(9*ee**2)*e_dyadic_e
    else
        D = dsde/3.D0*mat1 + (K0-2.D0*dsde/9.D0)*mat2
    end if
    ! find the stress
    do i = 1, 3
        stress(i) = 2.D0*se*dev_strain(i)/3/ee + K0*e_kk
        stress(i+3) = 2.D0*se*dev_strain(i+3)/3/ee
    end do

end subroutine mat_stiffness_3d


