module domain_class

  use contxt_class

  type geom_struc

     integer                              :: nb_dims   ! Dimension

     integer                              :: nb_elms
     integer                              :: nb_vers, nb_ver_elm

     real(rp)                             :: x_min , x_max
     real(rp)                             :: y_min , y_max
     real(rp)                             :: z_min , z_max

     integer    , pointer, dimension(:,:) :: ver_elm
     integer    , pointer, dimension(  :) :: ver_zon

     real(rp)   , pointer, dimension(:,:) :: ver_coo

     integer    , pointer, dimension(  :) :: elm_zon
     real(rp)  , pointer, dimension(  :) :: elm_mes
     integer    , pointer, dimension(:,:) :: elm_elm

     integer                              :: nb_facs, nb_fac_elm, nb_ver_fac

     integer    , pointer, dimension(:,:) :: fac_elm
     integer    , pointer, dimension(:,:) :: elm_fac
     integer    , pointer, dimension(:,:) :: ver_fac

     real(rp)   , pointer, dimension(:,:) :: fac_coo

     integer    , pointer, dimension(  :) :: fac_zon
     real(rp)  , pointer, dimension(  :) :: fac_mes

     integer                              :: nb_arts, nb_art_elm, nb_ver_art

     integer    , pointer, dimension(:,:) :: art_elm
     integer    , pointer, dimension(:,:) :: ver_art


     real(rp)   , pointer, dimension(:,:) :: art_coo

     integer    , pointer, dimension(  :) :: art_zon
     real(rp)  , pointer, dimension(  :) :: art_mes

  end type geom_struc

contains

  !\copyright

  ! Equipe Analyse Numerique
  ! Universite de Paris Sud
  ! 91405 ORSAY CEDEX
  ! FRANCE

  ! J. Laminie  F. Pascal

  ! Version 1.04 Juillet 1999

  !\purpose

  ! GM RECT 2D : Compute the mesh of a rectangle

  !\calling_sequence

  subroutine gm_rect_2d(geom, x1, y1, x2, y2, x3, y3, x4, y4, &
       nbseg1, nbseg2)
    
    !\parameters

    implicit none

    type(geom_struc), intent(inout) :: geom   ! The mesh
    integer         , intent(in   ) :: nbseg1 ! Number of segments between S1 S2
    integer         , intent(in   ) :: nbseg2 ! Number of segments between S1 S3
    real(rp)        , intent(in   ) :: x1, y1 ! Coordonnates of S1
    real(rp)        , intent(in   ) :: x2, y2 ! Coordonnates of S2
    real(rp)        , intent(in   ) :: x3, y3 ! Coordonnates of S3
    real(rp)        , intent(in   ) :: x4, y4 ! Coordonnates of S4
    
    !\guide

    ! The four vertices are defined as into the following figure.

    !  S4 +---------------------+ S3
    !     |                      \
    !     |                       \
    !     |                        \
    !     |                         \
    !  S1 +--------------------------+ S2

    ! nbseg1 is the number of segments between S1-S2 and S4-S3
    ! nbseg2 is the number of segments between S1-S4 and S2-S3

    !\end

    !--------------------------------
    ! Declarations of local variables
    !--------------------------------

    integer  :: i, ei, vi, is0
    integer  :: j
    
    real(rp) :: jx1, jx2, jy1, jy2, pi
    real(rp) :: pas_ix, pas_iy, pas_jx1, pas_jx2, pas_jy1, pas_jy2

    !---------------------
    ! Beginning of routine
    !---------------------

    !---------------------------------------------
    ! Compute the numbers of elements and vertices
    !---------------------------------------------

    geom%nb_dims    =  2
    geom%nb_vers    =    (nbseg1+1) * (nbseg2+1)
    geom%nb_elms    =     nbseg1    *  nbseg2 
    geom%nb_ver_elm = 4

    write(*,*) " Number of vertices             : ", geom%nb_vers
    write(*,*) " Number of elements             : ", geom%nb_elms
    write(*,*) " Number of vertices par element : ", geom%nb_ver_elm

    nullify(geom%ver_elm, geom%ver_coo, geom%ver_zon)

    allocate(geom%ver_elm(geom%nb_ver_elm,geom%nb_elms))
    allocate(geom%ver_coo(geom%nb_dims   ,geom%nb_vers))
    allocate(geom%ver_zon(                geom%nb_vers))

    !---------------------------------------------------
    ! Compute the number of the vertices of each element
    !---------------------------------------------------

    ei = 0

    do j = 1, nbseg2
       do i = 1, nbseg1

          vi = (j-1) * (nbseg1+1) +   i
          ei = ei + 1
          
          geom%ver_elm(1,ei) = vi
          geom%ver_elm(2,ei) = vi + 1
          geom%ver_elm(3,ei) = vi + nbseg1+2
          geom%ver_elm(4,ei) = vi + nbseg1+1

       end do
    end do

    !-----------------------------------------
    ! Compute the coordinnates of the vertices
    !-----------------------------------------

    ! ------------
    ! Uniform Mesh
    ! ------------

    pas_jx1 = (x4 - x1) / nbseg2 ; pas_jy1 = (y4 - y1) / nbseg2
    pas_jx2 = (x3 - x2) / nbseg2 ; pas_jy2 = (y3 - y2) / nbseg2

    vi = 0

    do j = 1, nbseg2+1

       jx1 = x1 + pas_jx1 * (j-1) ; jy1 = y1 + pas_jy1 * (j-1)
       jx2 = x2 + pas_jx2 * (j-1) ; jy2 = y2 + pas_jy2 * (j-1)

       pas_ix = (jx2 - jx1) / nbseg1 ; pas_iy = (jy2 - jy1) / nbseg1 

       do i = 1, nbseg1+1
          vi               = vi + 1
          geom%ver_coo(1,vi) = jx1 + pas_ix * (i-1)
          geom%ver_coo(2,vi) = jy1 + pas_iy * (i-1)
       end do

    end do

    !----------------------------
    ! Area number of the vertices
    !----------------------------

    !   8 3 3 3 3 3 7
    !   4 0 0 0 0 0 2 
    !   4 0 0 0 0 0 2
    !   4 0 0 0 0 0 2
    !   5 1 1 1 1 1 6

    geom%ver_zon(:) = 0
    
    is0 = geom%nb_vers - nbseg1-1
    
    geom%ver_zon(1    :nbseg1+1    ) = 1
    geom%ver_zon(1+is0:nbseg1+1+is0) = 3

    geom%ver_zon(nbseg1+1:geom%nb_vers-nbseg1:nbseg1+1) = 2
    geom%ver_zon(nbseg1+2:geom%nb_vers-nbseg1:nbseg1+1) = 4

    geom%ver_zon(1           ) = 5
    geom%ver_zon(nbseg1+1    ) = 6
    geom%ver_zon(geom%nb_vers) = 7
    geom%ver_zon(is0+1       ) = 8

    ! --------------
    ! End of routine
    ! --------------

  end subroutine gm_rect_2d

  subroutine gm_tria_2d(geom, x1, y1, x2, y2, x3, y3, x4, y4, &
       nbseg1, nbseg2)
    
    !\parameters

    implicit none

    type(geom_struc), intent(inout) :: geom   ! The mesh
    integer         , intent(in   ) :: nbseg1 ! Number of segments between v1 v2
    integer         , intent(in   ) :: nbseg2 ! Number of segments between v1 v3
    real(rp)        , intent(in   ) :: x1, y1 ! Coordonnates of v1
    real(rp)        , intent(in   ) :: x2, y2 ! Coordonnates of v2
    real(rp)        , intent(in   ) :: x3, y3 ! Coordonnates of v3
    real(rp)        , intent(in   ) :: x4, y4 ! Coordonnates of v4
    
    !\guide

    ! The four vertices are defined as into the following figure.

    !  v4 +---------------------+ v3
    !     |                      \
    !     |                       \
    !     |                        \
    !     |                         \
    !  v1 +--------------------------+ v2

    ! nbseg1 is the number of segments between v1-v2 and v4-v3
    ! nbseg2 is the number of segments between v1-v4 and v2-v3

    !\end

    !--------------------------------
    ! Declarations of local variables
    !--------------------------------

    integer  :: i, ei, vi, ei0, is0
    integer  :: j
    
    real(rp) :: jx1, jx2, jy1, jy2, pi
    real(rp) :: pas_ix, pas_iy, pas_jx1, pas_jx2, pas_jy1, pas_jy2

    !---------------------
    ! Beginning of routine
    !---------------------

    !---------------------------------------------
    ! Compute the numbers of elements and vertices
    !---------------------------------------------

    geom%nb_dims    =  2
    geom%nb_vers    =  (nbseg1+1) * (nbseg2+1)
    geom%nb_elms    =  2 *  nbseg1    *  nbseg2 
    geom%nb_ver_elm =  3

    write(*,*) " Number of vertices             : ", geom%nb_vers
    write(*,*) " Number of elements             : ", geom%nb_elms
    write(*,*) " Number of vertices par element : ", geom%nb_ver_elm

    nullify(geom%ver_elm, geom%ver_coo, geom%ver_zon)

    allocate(geom%ver_elm(geom%nb_ver_elm,geom%nb_elms))
    allocate(geom%ver_coo(geom%nb_dims   ,geom%nb_vers))
    allocate(geom%ver_zon(                geom%nb_vers))

    !---------------------------------------------------
    ! Compute the number of the vertices of each element
    !---------------------------------------------------

    ei = 0

    do j = 1, nbseg2
       do i = 1, nbseg1

          vi = (j-1) * (nbseg1+1) +   i
          ei = ei + 2
		  ei0 = ei -1
          
          geom%ver_elm(2,ei0) = vi
          geom%ver_elm(3,ei0) = vi + nbseg1+2
		  geom%ver_elm(1,ei0) = vi + nbseg1+1
          
		  geom%ver_elm(3,ei) = vi
		  geom%ver_elm(1,ei) = vi + 1
		  geom%ver_elm(2,ei) = vi + nbseg1+2

       end do
    end do

    !-----------------------------------------
    ! Compute the coordinnates of the vertices
    !-----------------------------------------

    ! ------------
    ! Uniform Mesh
    ! ------------

    pas_jx1 = (x4 - x1) / nbseg2 ; pas_jy1 = (y4 - y1) / nbseg2
    pas_jx2 = (x3 - x2) / nbseg2 ; pas_jy2 = (y3 - y2) / nbseg2

    vi = 0

    do j = 1, nbseg2+1

       jx1 = x1 + pas_jx1 * (j-1) ; jy1 = y1 + pas_jy1 * (j-1)
       jx2 = x2 + pas_jx2 * (j-1) ; jy2 = y2 + pas_jy2 * (j-1)

       pas_ix = (jx2 - jx1) / nbseg1 ; pas_iy = (jy2 - jy1) / nbseg1 

       do i = 1, nbseg1+1
          vi               = vi + 1
          geom%ver_coo(1,vi) = jx1 + pas_ix * (i-1)
          geom%ver_coo(2,vi) = jy1 + pas_iy * (i-1)
       end do

    end do

    !----------------------------
    ! Area number of the vertices
    !----------------------------

    !   8 3 3 3 3 3 7
    !   4 0 0 0 0 0 2 
    !   4 0 0 0 0 0 2
    !   4 0 0 0 0 0 2
    !   5 1 1 1 1 1 6

    geom%ver_zon(:) = 0
    
    is0 = geom%nb_vers - nbseg1-1
    
    geom%ver_zon(1    :nbseg1+1    ) = 1
    geom%ver_zon(1+is0:nbseg1+1+is0) = 3

    geom%ver_zon(nbseg1+1:geom%nb_vers-nbseg1:nbseg1+1) = 2
    geom%ver_zon(nbseg1+2:geom%nb_vers-nbseg1:nbseg1+1) = 4

    geom%ver_zon(1           ) = 5
    geom%ver_zon(nbseg1+1    ) = 6
    geom%ver_zon(geom%nb_vers) = 7
    geom%ver_zon(is0+1       ) = 8

    ! --------------
    ! End of routine
    ! --------------

  end subroutine gm_tria_2d


  !\resume

  !\purpose

  ! GM-PRINT : Print the Current Geom

  !\calling_sequence

  subroutine gm_print(txt, geom) !Generic name pr_var

    implicit none

    character(*)    , intent(in) :: txt
    type(geom_struc), intent(in) :: geom

    !\skip

    !-------------------
    ! Local Declarations
    !-------------------

    !------------------------
    ! Beginning of Subroutine
    !-----------------------

    write(*,*) txt//' Definition of the mesh'

    write(*,*) ' Space Dimension                : ', geom%nb_dims
  
    write(*,*) ' Number of elements             : ', geom%nb_elms
    write(*,*) ' Number of vertices             : ', geom%nb_vers
 !~   write(*,*) ' Number of faces                : ', geom%nb_facs
 !~   write(*,*) ' Number of edges                : ', geom%nb_arts

    write(*,*) ' Number of vectices per element : ', geom%nb_ver_elm
 !~   write(*,*) ' Number of faces    per element : ', geom%nb_fac_elm
 !~   write(*,*) ' Number of vectices per faces   : ', geom%nb_ver_fac
 !~   write(*,*) ' Number of edges    per element : ', geom%nb_art_elm
 !~   write(*,*) ' Number of vectices per edges   : ', geom%nb_ver_art
  
    if (associated(geom%ver_elm)) then ! we need first to see that the
       !                               ! pointer is associated with a 
       !                               ! target (see on pointers)
       write(*,*) ' Table of vectices per element  : ', geom%ver_elm
    end if
    if (associated(geom%ver_coo)) then  
       write(*,*) ' Tables of vectices coord       : ', geom%ver_coo
    end if
    if (associated(geom%ver_zon)) then  
       write(*,*) ' Tables of vectices zone        : ', geom%ver_zon
    end if

!~    if (associated(geom%elm_zon)) then  
!~       write(*,*) ' Tables of elements zone        : ', geom%elm_zon
!~    end if
!~    if (associated(geom%elm_mes)) then  
!~       write(*,*) ' Tables of elements mesure      : ', geom%elm_mes
!~    end if
!~    if (associated(geom%elm_elm)) then  
!~       write(*,*) ' Tables of elements around elem : ', geom%elm_elm
!~    end if
!~
!~    if (associated(geom%fac_elm)) then
!~       write(*,*) ' Table of faces per elements    : ', geom%fac_elm
!~    end if
!~    if (associated(geom%elm_fac)) then  
!~       write(*,*) ' Tables of elements per faces   : ', geom%elm_fac
!~    end if
!~
!~    if (associated(geom%ver_fac)) then  
!~       write(*,*) ' Tables of vectices par faces   : ', geom%ver_fac
!~    end if
!~    
!~    if (associated(geom%fac_zon)) then  
!~       write(*,*) ' Tables of faces    zone        : ', geom%fac_zon
!~    end if
!~    if (associated(geom%fac_mes)) then  
!~       write(*,*) ' Tables of faces    mesure      : ', geom%fac_mes
!~    end if
!~    if (associated(geom%fac_coo)) then  
!~       write(*,*) ' Tables of coordinates per face : ', geom%fac_coo
!~    end if
!~  
!~    if (associated(geom%art_elm)) then
!~       write(*,*) ' Table of vectices per edges    : ', geom%art_elm
!~    end if
!~    if (associated(geom%ver_art)) then  
!~       write(*,*) ' Tables of vectices par edges   : ', geom%ver_art
!~    end if
!~    
!~    if (associated(geom%art_zon)) then  
!~       write(*,*) ' Tables of edges    zone        : ', geom%art_zon
!~    end if
!~    if (associated(geom%art_mes)) then  
!~       write(*,*) ' Tables of edges    mesure      : ', geom%art_mes
!~    end if
!~    if (associated(geom%art_coo)) then  
!~       write(*,*) ' Tables of coordinates per edges: ', geom%art_coo
!~    end if

  !------------------
  ! End of Subroutine
  !------------------

  end subroutine gm_print
 
end module domain_class
