!     Subroutines for basic 3D hyperelastic elements
!     element identifier == 1005
!==========================SUBROUTINE el_hyperelast_3dbasic ==============================
subroutine el_hyperelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_3D
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
    integer      ::  n_points,kint,i,ii,jj,kk

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6), G(6,9)                    ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)

    real (prec)  ::  F(3,3), lB(3,3), lBInv(3,3), rC(3,3), FInv(3,3), Bstar(9,length_dof_array)
    real (prec)  ::  total_dof(length_dof_array/3, 3), dNdy(1:n_nodes,1:3)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array), S(3,length_dof_array/3)
    real (prec)  ::  Sigma(length_dof_array,length_dof_array)
    real (prec)  ::  IVec(6), BVec(6), BVecInv(6), PVec(3*n_nodes), SVec(3*n_nodes)
    real (prec)  ::  tempmat1(6,6), tempmat2(6,6), tempmat3(6,6), tempmat4(6,6)
    real (prec)  ::  tau(3,3), EGreen(3,3), EEulerian(3,3)                           ! Kirchoff stress
    real (prec)  ::  delta(3,3)

    real (prec)  ::  dxidx(3,3), determinant, J, detB         ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu1, K1                          ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         shear modulus
    !     element_properties(2)         bulk modulus

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))
    total_dof = 0.D0
    total_dof = transpose( reshape( dof_total + dof_increment, (/ 3, length_coord_array/3 /) ) )
    delta = reshape((/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), shape(delta))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    mu1 = element_properties(1)
    K1 = element_properties(2)

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        D = 0.d0
        F = 0.D0
        tempmat1 = 0.D0
        tempmat2 = 0.D0
        tempmat3 = 0.D0
        tempmat4 = 0.D0
        ! find deformation gradient F
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, n_nodes
                    F(ii,jj) = F(ii,jj) + total_dof(kk,ii)*dNdx(kk,jj)
                end do
                F(ii,jj) = F(ii,jj) + delta(ii,jj)
            end do
            tempmat1(ii,ii) = 1.D0
            tempmat1(ii+3,ii+3) = 0.5D0
        enddo
        lB = matmul(F,transpose(F))
        BVec = (/ lB(1,1), lB(2,2), lB(3,3), lB(1,2), lB(1,3), lB(2,3) /)
        call invert_small(lB,lBInv,detB)
        BVecInv = (/ lBInv(1,1), lBInv(2,2), lBInv(3,3), lBInv(1,2), lBInv(1,3), lBInv(2,3) /)
        IVec = (/ 1.D0, 1.D0, 1.D0, 0.D0, 0.D0, 0.D0 /)
        tempmat2 = spread(IVec,dim=2,ncopies=6)*spread(BVecInv,dim=1,ncopies=6)
        tempmat3 = spread(IVec,dim=2,ncopies=6)*spread(IVec,dim=1,ncopies=6)
        tempmat4 = spread(BVec,dim=2,ncopies=6)*spread(BVecInv,dim=1,ncopies=6)
        call invert_small(F,FInv,J)
        J = dsqrt(detB) ! Guarantee positive J
        D = mu1*J**(-2.D0/3)*tempmat1 + K1*J*(J-0.5D0)*tempmat2
        D = D + mu1*J**(-2.D0/3)/3.D0*( (BVec(1)+BVec(2)+BVec(3))/3.D0*tempmat2 - tempmat3 - tempmat4 )
        ! find dNdy
        dNdy = 0.D0
        do ii = 1, 3
            do jj = 1, 3
                dNdy(1:n_nodes,ii) = dNdy(1:n_nodes,ii) + dNdx(1:n_nodes,jj)*FInv(jj,ii)
            end do
        end do
        ! find Kirchoff stress and true stress
        tau(1:3,1:3) = mu1*J**(-2.D0/3)*(lB(1:3,1:3)-(lB(1,1)+lB(2,2)+lB(3,3))*delta(1:3,1:3)/3.D0) + K1*J*(J-1)*delta(1:3,1:3)
        stress(1:6) = (/ tau(1,1), tau(2,2), tau(3,3), tau(1,2), tau(1,3), tau(2,3) /)/J

        ! find the Lagrange/Eulerian strain
!        rC = matmul(transpose(F), F)
!        EGreen = 0.5D0*(rC - delta)
!        strain(1:6) = (/ EGreen(1,1), EGreen(2,2), EGreen(3,3), EGreen(1,2), EGreen(1,3), EGreen(2,3) /)
        EEulerian = 0.5D0*(delta - matmul(transpose(FInv),FInv))
        strain(1:6) = (/ EEulerian(1,1), EEulerian(2,2), EEulerian(3,3), EEulerian(1,2), EEulerian(1,3), EEulerian(2,3) /)

        ! find G
        G = 0.D0
        G = reshape( (/ 2.D0*lB(1,1), 0.D0, 0.D0, 2.D0*lB(1,2), 2.D0*lB(1,3), 0.D0, &
                        0.D0, 2.D0*lB(2,2), 0.D0, 2.D0*lB(1,2), 0.D0, 2.D0*lB(2,3), &
                        0.D0, 0.D0, 2.D0*lB(3,3), 0.D0, 2.D0*lB(1,3), 2.D0*lB(2,3), &
                        2.D0*lB(1,2), 0.D0, 0.D0, 2.D0*lB(2,2), 2.D0*lB(2,3), 0.D0, &
                        0.D0, 2.D0*lB(1,2), 0.D0, 2.D0*lB(1,1), 0.D0, 2.D0*lB(1,3), &
                        2.D0*lB(1,3), 0.D0, 0.D0, 2.D0*lB(2,3), 2.D0*lB(3,3), 0.D0, &
                        0.D0, 0.D0, 2.D0*lB(1,3), 0.D0, 2.D0*lB(1,1), 2.D0*lB(1,2), &
                        0.D0, 2.D0*lB(2,3), 0.D0, 2.D0*lB(1,3), 0.D0, 2.D0*lB(3,3), &
                        0.D0, 0.D0, 2.D0*lB(2,3), 0.D0, 2.D0*lB(1,2), 2.D0*lB(2,2) /), shape(G))
        ! find B
        B = 0.D0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        ! find Bstar
        Bstar = 0.D0
        Bstar(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bstar(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bstar(3,3:3*n_nodes-0:3) = dNdy(1:n_nodes,3)
        Bstar(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        Bstar(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        Bstar(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        Bstar(7,3:3*n_nodes-0:3) = dNdy(1:n_nodes,1)
        Bstar(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        Bstar(9,3:3*n_nodes-0:3) = dNdy(1:n_nodes,2)

        ! find Sigma
        S = reshape(matmul(transpose(B),stress*J),(/3,length_dof_array/3/))
        do i = 1,n_nodes
            PVec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(PVec,dim=1,ncopies=3)
            SVec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*i-2:3*i,1:3*n_nodes) = spread(SVec,dim=1,ncopies=3)
        end do
        Sigma = Pmat*transpose(Smat)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress*J)*w(kint)*determinant
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul( transpose(B(1:6,1:3*n_nodes)), matmul( D(1:6,1:6), &
              matmul( G(1:6,1:9), Bstar(1:9,1:3*n_nodes) ) ) )*w(kint)*determinant &
            - Sigma(1:3*n_nodes,1:3*n_nodes)*w(kint)*determinant

    end do

    return
end subroutine el_hyperelast_3dbasic



!==========================SUBROUTINE fieldvars_hyperelast_3dbasic ==============================
subroutine fieldvars_hyperelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : principalvals33

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
    logical      ::  strcmp
  
    integer      ::  n_points,kint,i,ii,jj,kk,k

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6), G(6,9)                    ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)

    real (prec)  ::  F(3,3), lB(3,3), lBInv(3,3), rC(3,3), FInv(3,3), Bstar(9,length_dof_array)
    real (prec)  ::  total_dof(length_dof_array/3,3), dNdy(1:n_nodes,1:3)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array), S(3,length_dof_array/3)
    real (prec)  ::  Sigma(length_dof_array,length_dof_array)
    real (prec)  ::  IVec(6), BVec(6), BVecInv(6), PVec(3*n_nodes), SVec(3*n_nodes)
    real (prec)  ::  tempmat1(6,6), tempmat2(6,6), tempmat3(6,6), tempmat4(6,6)
    real (prec)  ::  tau(3,3), EGreen(3,3), EEulerian(3,3)              ! Kirchoff stress and Lagrange strain
    real (prec)  ::  pstrain(3), pstress(3)             ! principle strain and stress
    real (prec)  ::  delta(3,3)

    real (prec)  ::  dxidx(3,3), determinant, J, detB         ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu1, K1                          ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         shear modulus
    !     element_properties(2)         bulk modulus

    x = reshape(element_coords,(/3,length_coord_array/3/))
    total_dof = 0.D0
    total_dof = transpose( reshape( dof_total + dof_increment, (/ 3, length_coord_array/3 /) ) )
    delta = reshape((/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), shape(delta))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    mu1 = element_properties(1)
    K1 = element_properties(2)

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        D = 0.d0
        F = 0.D0
        ! find deformation gradient F
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, n_nodes
                    F(ii,jj) = F(ii,jj) + total_dof(kk,ii)*dNdx(kk,jj)
                end do
                F(ii,jj) = F(ii,jj) + delta(ii,jj)
            end do
            tempmat1(ii,ii) = 1.D0
            tempmat1(ii+3,ii+3) = 0.5D0
        enddo

        lB = matmul(F(1:3,1:3),transpose(F(1:3,1:3)))
        BVec = (/ lB(1,1), lB(2,2), lB(3,3), lB(1,2), lB(1,3), lB(2,3) /)
        call invert_small(lB,lBInv,detB)
        BVecInv = (/ lBInv(1,1), lBInv(2,2), lBInv(3,3), lBInv(1,2), lBInv(1,3), lBInv(2,3) /)
        IVec = (/ 1.D0, 1.D0, 1.D0, 0.D0, 0.D0, 0.D0 /)
        tempmat2 = spread(IVec,dim=2,ncopies=6)*spread(BVecInv,dim=1,ncopies=6)
        tempmat3 = spread(IVec,dim=2,ncopies=6)*spread(IVec,dim=1,ncopies=6)
        tempmat4 = spread(BVec,dim=2,ncopies=6)*spread(BVecInv,dim=1,ncopies=6)
        call invert_small(F,FInv,J)
        J = dsqrt(detB) ! Guarantee positive J
        D = mu1*J**(-2.D0/3)*tempmat1 + K1*J*(J-0.5D0)*tempmat2
        D = D + mu1*J**(-2.D0/3)/3.D0*( (BVec(1)+BVec(2)+BVec(3))/3.D0*tempmat2 - tempmat3 - tempmat4 )

        ! find dNdy
        dNdy = 0.D0
        do ii = 1, 3
            do jj = 1, 3
                dNdy(1:n_nodes,ii) = dNdy(1:n_nodes,ii) + dNdx(1:n_nodes,jj)*FInv(jj,ii)
            end do
        end do
        ! find Kirchoff stress and true stress
        tau(1:3,1:3) = mu1*J**(-2.D0/3)*(lB(1:3,1:3)-(lB(1,1)+lB(2,2)+lB(3,3))*delta(1:3,1:3)/3.D0) + K1*J*(J-1)*delta(1:3,1:3)
        stress(1:6) = (/ tau(1,1), tau(2,2), tau(3,3), tau(1,2), tau(1,3), tau(2,3) /)/J

        ! find the Lagrange/Eulerian strain
!        rC = matmul(transpose(F), F)
!        EGreen = 0.5D0*(rC - delta)
!        strain(1:6) = (/ EGreen(1,1), EGreen(2,2), EGreen(3,3), EGreen(1,2), EGreen(1,3), EGreen(2,3) /)
        EEulerian = 0.5D0*(delta - matmul(transpose(FInv),FInv))
        strain(1:6) = (/ EEulerian(1,1), EEulerian(2,2), EEulerian(3,3), EEulerian(1,2), EEulerian(1,3), EEulerian(2,3) /)

        ! find principle stress and strain
        pstress = principalvals33(stress)
        pstrain = principalvals33(strain)

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
            else if (strcmp(field_variable_names(k),'E11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'PS1',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + pstress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'PS2',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + pstress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'PS3',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + pstress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'PE1',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + pstrain(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'PE2',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + pstrain(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'PE3',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + pstrain(3)*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do

    return
end subroutine fieldvars_hyperelast_3dbasic

