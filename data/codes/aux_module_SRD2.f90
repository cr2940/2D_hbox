! CODES FOR SRD:


module aux_module_SRD2
  use aux_module_hbox
  use hbox_layer
  implicit none

contains
  subroutine N_hboxes(area_hboxes,N_hbox,whose,mx,dx,dy)
    ! the overlap count for hboxes when applying SRD for transverse fluctuationn calucltions
    implicit none
    real(8) :: area_hboxes(:),dx,dy
    integer :: mx
    integer :: i, N, N_hbox(mx*4),whose(mx*4)
    real(8):: area_chk
    N = size(area_hboxes(:),1)
    N_hbox = 1
    do i=1,N-1
      if (area_hboxes(i) .lt. 0.5d0*dx*dy) then
        area_chk = area_hboxes(i) + area_hboxes(i+1)
        if (area_chk .ge. 0.5d0*dx*dy) then
          whose(i) =  i+1
          N_hbox(i+1) = N_hbox(i+1) + 1
        end if
      else
        whose(i) = i
      end if
    end do
    if (area_hboxes(N) .lt. 0.5d0*dx*dy) then
      area_chk = area_hboxes(i) + area_hboxes(i-1)
      if (area_chk .ge.0.5d0*dx*dy ) then
        whose(N) = N-1
        N_hbox(N-1) = N_hbox(N-1) + 1
      end if
    end if

  end subroutine

  subroutine ind_in_inds2(ind,inds,truea,loc)! result(truea)
    ! finding index in list of indices and the location if it exists
    implicit none
    integer :: ind(2),inds(:,:), i,N,loc
    logical :: truea
    N = size(inds,2)
    truea=.false.
    loc = 0
    do i=1,N
      if (ind(1)==inds(1,i) .and. ind(2)==inds(2,i)) then
        truea = .true.
        loc = i
        exit
      end if
    end do
  end subroutine


   subroutine frag_tallies(m,mx,num_frags,area_frags,index_frags,uniq_indices,indices_data,indices_data_area,&
         counts,i2)
     implicit none
     ! gets all the indices affected by hbox and their associated data (the number of hboxes covering it, which hboxes covering it, their areas)
     ! loop through each index, make a list of unique indices
     ! as you go through, record the index of hbox and the area frags associated
     ! if you find reoccuring index, go to its slot and record the index of hbox and area associated
    integer :: i,j,m,i2,mx ! i2  is the number of affected cells by the hboxes , the limit of uniq_indices
      integer :: uniq_indices(2,8*mx) ! all unique fragment indices
      integer :: indices_data(4,8*mx) ! indices of hboxes covering the unique index of uniq_indices
      real(8) :: indices_data_area(4,8*mx) ! area of frags associated with indices of hboxes of indices_data
      integer :: num_frags(:), index_frags(:,:,:) ! number of fragments in each hbox ; the index location of each fragment (<=4)
      real(8) :: area_frags(:,:)
      logical :: tof
      integer :: loc,counts(8*mx) ! for counting how many hboxes cover the fragment index cell
      counts = 1
      indices_data_area = 0.d0
      i2 = 1 ! unique frag indices counter i.e. counter for 8*mx
      do i = 1,m ! loop over each hbox and unbox its frags
        if (i.eq.1) then
          do j=1,num_frags(1)
            uniq_indices(:,i2) = index_frags(:,j,i)
            indices_data(1,i2) = i
            indices_data_area(1,i2) = area_frags(j,i)
            i2=i2+1
          end do
        else
          do j=1,num_frags(i)
            call ind_in_inds2(index_frags(:,j,i),uniq_indices(:,1:i2-1),tof,loc)
            if (.not. tof) then
              uniq_indices(:,i2) = index_frags(:,j,i)
              indices_data(1,i2) = i
              indices_data_area(1,i2) = area_frags(j,i)
              i2 = i2 + 1
            else
              ! print *, "LOC:",loc
              counts(loc) = counts(loc)+ 1
              indices_data(counts(loc),loc) = i
              indices_data_area(counts(loc),loc) = area_frags(j,i)
            end if
          end do
        end if
      enddo
      i2 = i2 -1
    end subroutine

    subroutine residual_increments_inds(dx,dy,mx,mbc,area_ij,indices_data_area,hbox_area,&
      counts,i2,uniq_indices,indices_data,comp_cov_inds,C2)
      ! gets the fully covered (by hbox) cell indices
      ! go through the unique indices of cells affected by hbox, add up the indices data area for each cell
      ! compare with their actual area (area_ij), and if equal to the actual area, then include in output list of inds.
      ! in step2ds, as you compute the qhbox valus, also get the residual increments (if not in this list)
      ! and use them as hbox outflux.
      implicit none
      integer :: mx,mbc,uniq_indices(:,:),counts(:),indices_data(:,:)
      integer :: i2,i,j
      real(8) :: indices_data_area(:,:),hbox_area(:)
      real(8) :: area_ij(1-mbc:mx+mbc,1-mbc:mx+mbc),dx,dy
      integer :: ind(2),N,C2,comp_cov_inds(2,mx*8) ! the indices of fully covered cells, if as you calcualate hbox values, if the index is not found here, must use as part of outer flux
      real(8) :: area_sum, area_actual
      C2=1
      do i=1,i2
        ind = uniq_indices(:,i)
        N = counts(i)
        area_sum = 0.d0
        ! print *, "count:", N
        ! print *, "indices data:",indices_data(:,i)
        do j=1,N
          ! if (indices_data(j,i).ne.0) then
            area_sum = area_sum + hbox_area(indices_data(j,i))*indices_data_area(j,i)
          ! end if
        end do
        area_actual = area_ij(ind(1),ind(2))*dx*dy
        ! print *, "diff:", area_sum,area_actual,abs(area_sum-area_actual)
        if (abs(area_sum - area_actual).le.1d-10) then
          comp_cov_inds(:,C2) = ind
          C2 = C2 + 1
        end if
      end do
      C2 = C2-1
    end subroutine

  subroutine cut_cells_find(x_0,y_0,x_e,y_e,dx,dy,xlower,ylower,mx,my,ii,jj,intersections,N_cells)   !  O
    ! returns ii, jj which are 4*mx (long) arrays that contain i-index and j-index respectively
    ! of cut cells, and N which is the actual number of cut cells (so go from ii(1:N),jj(1:N))
    implicit none

    real(8) :: x_0,y_0,x_e,y_e,dx,dy,xlower,ylower
    integer :: mx,my,ii(4*mx),jj(4*mx),N_cells,N
    real(8) :: intersections(2,4*mx)

    ! local
    integer :: i,j,k,i_0,j_0,size1,mbc
    integer, allocatable :: indexes(:)
    real(8) :: xe(-2:mx+2),ye(-2:my+2),slope_bar,theta,dist_x,dist_y
    real(8), allocatable :: intersect_0(:,:), intersect_top(:,:),midpoints(:,:)
    xe = (/(xlower+i*dx, i=-2,mx+2)/)
    ye = (/(ylower+i*dy, i=-2,my+2)/)
    mbc =2
    ! wall parameters:
    if (x_e-x_0 .ne. 0.d0) then
      slope_bar = (y_e-y_0)/(x_e-x_0)
      theta = atan(slope_bar)
    else
      i=-2
      j=-2
      do while (x_0 > xe(i))
        i = i +1
      end do
      do while  (y_0 > ye(j))
        j=j+1
      end do
      ii = i-1
      jj = (/(j-1+k,k=1,my)/)
      N_cells = my
      return
    end if

    ! get intersections between barrier and grid lines :
    do k=-2,mx+2
      if (xe(k)-tol .le. x_e .and. xe(k)+tol .ge. x_0) then
        dist_x = x_e - xe(k)
        dist_y = dist_x * tan(pi-theta)
        call AddToList_verts(intersect_0,(/xe(k),y_e+dist_y/))
      end if
    end do

    do k=-2,my+2
      if (ye(k)-tol.le.max(y_e,y_0) .and. ye(k)+tol.ge.min(y_0,y_e)) then
        dist_y = y_0 - ye(k)
        dist_x = dist_y/tan(pi-theta)
        call AddToList_verts(intersect_0,(/x_0+dist_x,ye(k)/))
      end if
    end do
    ! treat array (get rid of duplicate and sort)
       intersect_top = remove_dups(intersect_0)
       ! print*, "Intersections: x", intersect_top(1,:)
       ! print*, "intersections: y", intersect_top(2,:)
       size1 = size(intersect_top,2)
       N=size1
      allocate(indexes(size1))
      call KB07AD(intersect_top(1,:),size1,indexes)
      intersect_top(2,:) = intersect_top(2,indexes)
      ! print*,"Intersections; ", intersect_top(1,:)
      ! print*,"Intersections: ", intersect_top(2,:)
      intersections(:,1:N) = intersect_top(:,1:N)
      ! for each pair of neighboring intersections, find which cell it cuts
         !find the midpoint of the intersections and see which cell they fit in
         ! start with i_0,j_0 and walk up or down
        allocate(midpoints(N-1,2))
        do i=1, N-1
          midpoints(i,1:2) = (/0.5d0*(intersect_top(1,i)+intersect_top(1,i+1)),&
          0.5d0*(intersect_top(2,i)+intersect_top(2,i+1))/)
        end do
        ! write(*,*) "midpoint ", midpoints

      ! put them into ii and jj
      ! the first cell where the barrier starts is always cut:
      do i=1-mbc,mx+mbc
        if (midpoints(1,1) .gt. xe(i) .and. midpoints(1,1) .lt. xe(i+1)) then
          i_0 = i
          exit
        end if
      end do
      do i = 1-mbc,my+mbc
        if (midpoints(1,2) .gt. ye(i) .and. midpoints(1,2) .lt. ye(i+1)) then
          j_0 = i
          exit
        end if
      end do

      ii(1) = i_0
      jj(1) = j_0
      j = 1 ! counter for cut cells
      do i = 1, N-2
        if (midpoints(i+1,1) .lt. xe(ii(i)) .and. midpoints(i+1,1) .gt. xe(ii(i)-1)) then
          ii(i+1) = ii(i) - 1
        else if (midpoints(i+1,1) .gt. xe(ii(i)+1) .and. midpoints(i+1,1) .lt. xe(ii(i)+2)) then
          ii(i+1) = ii(i) + 1
        else if (midpoints(i+1,1) .gt. xe(ii(i)) .and. midpoints(i+1,1) .lt. xe(ii(i)+1)) then
          ii(i+1) = ii(i)
        end if
        if (midpoints(i+1,2) .lt. ye(jj(i)) .and. midpoints(i+1,2) .gt. ye(jj(i)-1)) then
          jj(i+1) = jj(i) - 1
        else if (midpoints(i+1,2) .gt. ye(jj(i)+1) .and. midpoints(i+1,2) .lt. ye(jj(i)+2)) then
          jj(i+1) = jj(i) + 1
        else if (midpoints(i+1,2) .gt. ye(jj(i)) .and. midpoints(i+1,2) .lt. ye(jj(i)+1)) then
          jj(i+1) = jj(i)
        end if
      end do
      ! N is number of intersections, N_cells number of cut cells.
      N_cells=N-1
      ! fortran indexing
      ii = ii + 1
      jj = jj + 1
   end subroutine


    subroutine small_cells_geom(ii,jj,intersections,mx,my,dx,dy,xlower,ylower,N_cells,&        ! O
      type_supper,type_sunder,lengths_supper,lengths_sunder,area_supper,area_sunder)

      implicit none
      integer :: ii(:),jj(:),mx,my,N_cells,i,i_0,j_0,j,k1,k2
      real(8) :: intersections(:,:),xlower,ylower,x1,y1,x2,y2,m,dx,dy
      real(8) :: x_0,x_e,y_0,y_e,y_int
      integer :: type_supper(4*mx), type_sunder(4*mx)
      real(8) :: area_supper(4*mx), area_sunder(4*mx)
      real(8) :: lengths_supper(5,4*mx), lengths_sunder(5,4*mx)
      real(8) :: coord_supper(2,6,4*mx), coord_sunder(2,6,4*mx)
      real(8) :: coords(2,4),  xe(-2:mx+2),ye(-2:my+2)
      integer :: uo_array(4)

      !initialization:
      lengths_supper = 0.d0
      lengths_sunder = 0.d0
      uo_array = (/3,4,1,2/) ! the order of cells' vertices checking for under-small cells' vertices

      ! edges:
      xe = (/(xlower+i*dx, i=-2,mx+2)/)
      ye = (/(ylower+i*dy, i=-2,my+2)/)

      do i=1,N_cells
        k1 = 1 ! counter for num sides for upper small cell
        k2 = 1 ! counter for num sides for under small cell
        ! the index
        i_0 = ii(i)
        j_0 = jj(i)
        ! the box
        x1 = xe(i_0-1)
        x2 = xe(i_0)
        y1 = ye(j_0-1)
        y2 = ye(j_0)
        coords(1,:) = (/x1,x1,x2,x2/)
        coords(2,:) = (/y1,y2,y2,y1/)
        ! the bar
        x_0 = intersections(1,i)
        y_0 = intersections(2,i)
        x_e = intersections(1,i+1)
        y_e = intersections(2,i+1)
        m = (y_e-y_0)/(x_e-x_0)
        y_int = y_0-m*x_0
        ! the coordinates of small cells (upper and under)
        coord_supper(:,1,i) = (/x_0,y_0/)
        coord_supper(:,2,i) = (/x_e,y_e/)
        k1 = k1+2
        coord_sunder(:,1,i) = (/x_0,y_0/)
        coord_sunder(:,2,i) = (/x_e,y_e/)
        k2 = k2+2
        do j=1,4
          if (m*coords(1,4-(j-1))+y_int < coords(2,4-(j-1))) then
            coord_supper(:,k1,i) = coords(:,4-(j-1))
            k1= k1 + 1
          end if
          if (m*coords(1,uo_array(j))+y_int > coords(2,uo_array(j))) then
            coord_sunder(:,k2,i) = coords(:,uo_array(j))
            k2 = k2 + 1
          end if
        end do
        coord_supper(:,k1,i) = (/x_0,y_0/)
        coord_sunder(:,k2,i) = (/x_0,y_0/)

        ! side lengths
        do j=1,k1-1
          lengths_supper(j,i) = dist(coord_supper(:,j,i),coord_supper(:,j+1,i))/dx
        end do
        do j=1,k2-1
          lengths_sunder(j,i) = dist(coord_sunder(:,j,i),coord_sunder(:,j+1,i))/dx
        end do
        k1 = k1-1 ! actual number of sides
        k2 = k2-1
        ! type 1-8
        if (k1.eq.5 .and. m>0) then
          type_supper(i) = 1
          type_sunder(i) = 1
        else if (k1.eq.4 .and. m>0) then
          if (abs(m)<1) then
            type_supper(i) = 3
            type_sunder(i) = 3
          else
            type_supper(i) = 2
            type_sunder(i) = 2
          endif
        else if (k1.eq.3 .and. m>0) then
          type_supper(i) = 4
          type_sunder(i) = 4
        else if (k1.eq.5 .and. m<0) then
          type_supper(i) = 8
          type_sunder(i) = 8
        else if (k1.eq.4 .and. m<0) then
          if (abs(m)<1) then
            type_supper(i) = 6
            type_sunder(i) = 6
          else
            type_supper(i) = 7
            type_sunder(i) = 7
          end if
        else if (k1.eq.3 .and. m<0) then
          type_supper(i) = 5
          type_sunder(i) = 5
        end if
        ! area of small cells
        area_supper(i) = area_polygon(coord_supper(1,1:k1+1,i), coord_supper(2,1:k1+1,i))/(dx*dy)
        area_sunder(i) = area_polygon(coord_sunder(1,1:k2+1,i), coord_sunder(2,1:k2+1,i))/(dx*dy)
      end do


    end subroutine


    subroutine SRD_undercells(N_cells,ii,jj,mx,my,area_sunder,x_0,y_0,x_e,y_e,dx,dy,&           !  [      ]
      unS_cells_i,unS_cells_j,N_ij,all_undercells_i,all_undercells_j,k_count)
      implicit none

      integer :: mx,my,ii(:),jj(:),N_cells,N_ij(-1:mx+2,-1:my+2)
      real(8) :: area_sunder(:),x_0,y_0,x_e,y_e,m,dx,dy
      integer :: unS_cells_i(mx*4), unS_cells_j(mx*4) ! THESE will tell you who the neighbor is for small cells
      integer :: all_undercells_i(mx*4), all_undercells_j(mx*4) ! THESE will just give you in order the indices of all affected cells by nhood inclusions

      integer :: i,j,k_count
      logical :: TF
      ! N_ij is the matrix of number of neighborhodds each cell belongs to
      N_ij = 1 ! i.e. everybody is its own neighbor

      ! slope of barrier :
      m = (y_e-y_0)/(x_e-x_0)
      ! initialization of indices that will have to be looped over for SRD updates
      unS_cells_i = huge(1)
      unS_cells_j = huge(1)

      if (abs(m).le.1.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.5d-13))then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i)
          else if (area_sunder(i) .lt. 0.5d0) then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i) - 1
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(unS_cells_i(i),unS_cells_j(i)) = N_ij(unS_cells_i(i),unS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.lt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.5d-13))then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i)
          elseif (area_sunder(i) .lt. 0.5d0) then
            unS_cells_i(i) = ii(i) - 1
            unS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(unS_cells_i(i),unS_cells_j(i)) = N_ij(unS_cells_i(i),unS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.gt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_sunder(i) .gt. 0.5d0 ).or.(abs(area_sunder(i)-0.5d0).lt.5d-13))then
            unS_cells_i(i) = ii(i)
            unS_cells_j(i) = jj(i)
          elseif (area_sunder(i) .lt. 0.5d0) then
            unS_cells_i(i) = ii(i) + 1
            unS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(unS_cells_i(i),unS_cells_j(i)) = N_ij(unS_cells_i(i),unS_cells_j(i)) + 1
        enddo
      end if

      k_count = 1
      do i = 1,N_cells
        all_undercells_i(k_count) = ii(i)
        all_undercells_j(k_count) = jj(i)
        k_count = k_count + 1
        if (area_sunder(i) .lt. 0.5d0 .and. abs(area_sunder(i)-0.5d0).gt.5d-13) then
          TF = check_intpair_in_set((/unS_cells_i(i),unS_cells_j(i)/),&
          ii(max(1,i-1):min(i+1,N_cells)),jj(max(1,i-1):min(i+1,N_cells)))
          if (.not. TF) then
            all_undercells_i(k_count) = unS_cells_i(i)
            all_undercells_j(k_count) = unS_cells_j(i)
            k_count = k_count + 1
          end if
        end if
      end do

      k_count = k_count - 1 ! number of actually affected cells for which SRD update applies

    end subroutine

    subroutine SRD_uppercells(N_cells,ii,jj,mx,my,area_supper,x_0,y_0,x_e,y_e,dx,dy,&        !  [     ]
      upS_cells_i,upS_cells_j,N_ij,all_uppercells_i,all_uppercells_j,k_count)
      implicit none

      integer :: mx,my,ii(:),jj(:),N_cells,N_ij(-1:mx+2,-1:my+2)
      real(8) :: area_supper(:),x_0,y_0,x_e,y_e,m,dx,dy
      integer :: upS_cells_i(mx*4), upS_cells_j(mx*4) ! THESE will tell you who the neighbor is for small cells
      integer :: all_uppercells_i(mx*4), all_uppercells_j(mx*4) ! THESE will just give you in order the indices of all affected cells by nhood inclusions  (AKA the YELLOW ARRAY OF INDICES)

      integer :: i,j,k_count
      logical :: TF
      ! N_ij is the matrix of number of neighborhodds each cell belongs to
      N_ij = 1

      ! slope of barrier :
      m = (y_e-y_0)/(x_e-x_0)
      ! initialization of indices that will have to be looped over for SRD updates
      upS_cells_i = huge(1)
      upS_cells_j = huge(1)

      if (abs(m).le.1.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.5d-13))then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i)
          else if (area_supper(i) .lt. 0.5d0) then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i) + 1
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(upS_cells_i(i),upS_cells_j(i)) = N_ij(upS_cells_i(i),upS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.lt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.5d-13))then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i)
          elseif (area_supper(i) .lt. 0.5d0) then
            upS_cells_i(i) = ii(i) + 1
            upS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(upS_cells_i(i),upS_cells_j(i)) = N_ij(upS_cells_i(i),upS_cells_j(i)) + 1
        enddo
      else if (abs(m).gt.1.d0 .and. m.gt.0.d0) then
        do i=1,N_cells
          N_ij(ii(i),jj(i)) = 0
          if ((area_supper(i) .gt. 0.5d0 ).or.(abs(area_supper(i)-0.5d0).lt.5d-13))then
            upS_cells_i(i) = ii(i)
            upS_cells_j(i) = jj(i)
          elseif (area_supper(i) .lt. 0.5d0) then
            upS_cells_i(i) = ii(i) - 1
            upS_cells_j(i) = jj(i)
            N_ij(ii(i),jj(i)) = N_ij(ii(i),jj(i)) + 1
          end if
          N_ij(upS_cells_i(i),upS_cells_j(i)) = N_ij(upS_cells_i(i),upS_cells_j(i)) + 1
        enddo
      end if

      k_count = 1
      do i = 1,N_cells
        all_uppercells_i(k_count) = ii(i)
        all_uppercells_j(k_count) = jj(i)
        k_count = k_count + 1
        if (area_supper(i) .lt. 0.5d0.and. abs(area_supper(i)-0.5d0).gt.5d-13) then
          TF = check_intpair_in_set((/upS_cells_i(i),upS_cells_j(i)/),&
          ii(max(1,i-1):min(i+1,N_cells)),jj(max(1,i-1):min(i+1,N_cells)))
          if (.not. TF) then
            all_uppercells_i(k_count) = upS_cells_i(i)
            all_uppercells_j(k_count) = upS_cells_j(i)
            k_count = k_count + 1
          end if
        end if
      end do

      k_count = k_count - 1 ! number of actually affected cells for which SRD update applies

    end subroutine


    function check_intpair_in_set(intpair,setx,sety) result (TF)
      implicit none
      integer :: intpair(2), setx(:),sety(:),i,n
      logical :: TF
      TF = .false.
      n = size(setx)
      do i =1,n
        if (intpair(1) .eq. setx(i) .and. intpair(2) .eq. sety(i)) then
          TF = .true.
          return
        end if
      end do
    end function


    subroutine area_cells(ii,jj,mx,my,mbc,area_sunder,area_supper,up_area_ij,un_area_ij)
      implicit none
      integer :: ii(:),jj(:),mx,my,mbc
      real(8) :: area_sunder(:),area_supper(:)
      real(8),intent(out) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
      ! local
      integer :: N,i

      ! initialization of area matrix
      up_area_ij = 1.d0
      un_area_ij = 1.d0

      ! change the small cell ones
      N = size(ii)
      do i=1,N
        up_area_ij(ii(i),jj(i)) = area_supper(i)
        un_area_ij(ii(i),jj(i)) = area_sunder(i)
      end do
    end subroutine


    subroutine SRD_update_correct(qnew,qnew2,all_undercells_i,all_undercells_j,&      !  [     ]
                    all_uppercells_i,all_uppercells_j,ii,jj,N_ij_up,N_ij_un,&
                    unS_cells_i,unS_cells_j,upS_cells_i,upS_cells_j,up_area_ij,&
                    un_area_ij,mx,my,mbc)
      implicit none
      real(8) :: qnew(3,1-mbc:mx+mbc,1-mbc:my+mbc),qnew2(3,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: all_undercells_i(:),all_uppercells_i(:),all_undercells_j(:)
      integer :: all_uppercells_j(:),ii(:),jj(:),N_ij_up(1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: N_ij_un(1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: unS_cells_i(:),unS_cells_j(:),upS_cells_i(:),upS_cells_j(:)
      real(8) :: up_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc),un_area_ij(1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc
      ! local
      integer :: i,N,i0,j0,num_ij,i_n,j_n, array(mx*4),k,array2(mx*4)
      real(8) :: beta,alpha,nhood_val(3)

      ! DO UNDER SIDE FIRST:
      N = size(all_undercells_i)
      k = 1
      do i=1,N
        i0 = all_undercells_i(i)
        j0 = all_undercells_j(i)
        if (i0.eq.ii(k) .and. j0.eq.jj(k)) then
          array(i) = k
          k = k +1
        else
          array(i)=k-1
        end if
        ! print *, " I0,J0:", i0,j0
        ! print *, "ARRAY:", array(i)
        num_ij = N_ij_un(i0,j0)
        ! print *, "NUM_IJ:", num_ij
        if (num_ij .eq. 2) then
          qnew(:,i0,j0) = qnew(:,i0,j0)/2.d0
        else if (num_ij .eq. 1) then
          ! get neighborhood value:
          ! print *, "area: UN", un_area_ij(i0,j0)
          if (un_area_ij(i0,j0).lt.0.5d0 .and. abs(un_area_ij(i0,j0)-0.5d0).gt.5d-13)then
            i_n = unS_cells_i(array(i))
            j_n = unS_cells_j(array(i))
            ! print*, "IN,JN:", i_n,j_n
            beta = un_area_ij(i_n,j_n)
            alpha = un_area_ij(i0,j0)
            nhood_val = (alpha+beta/2.d0)**(-1) * (alpha*qnew(:,i0,j0)+&
                 beta/2.d0*qnew(:,i_n,j_n))
            qnew(:,i0,j0) = nhood_val
            qnew(:,i_n,j_n) = qnew(:,i_n,j_n) + nhood_val/2.d0
          end if
        end if
      end do

      ! DO UPPER SIDE FIRST:
      N = size(all_uppercells_i)
      k=1
      do i=1,N
        i0 = all_uppercells_i(i)
        j0 = all_uppercells_j(i)
        if (i0.eq.ii(k) .and. j0.eq.jj(k)) then
          array2(i) = k
          k =k +1
        else
          array2(i)=k-1
        end if
        ! print*, "ARRAY2:", array2(i)
        num_ij = N_ij_up(i0,j0)
        ! print *, "NUM_IJ:", num_ij
        if (num_ij .eq. 2) then
          qnew2(:,i0,j0) = qnew2(:,i0,j0)/2.d0
        else if (num_ij .eq. 1) then
          ! print *, "area: UP", up_area_ij(i0,j0)
          ! get neighborhood value:
          if (up_area_ij(i0,j0).lt.0.5d0 .and. abs(up_area_ij(i0,j0)-0.5d0).gt.5d-13)then
            i_n = upS_cells_i(array2(i))
            j_n = upS_cells_j(array2(i))
            ! print *, "In,jn:",i_n,j_n
            beta = up_area_ij(i_n,j_n)
            alpha = up_area_ij(i0,j0)
            nhood_val = (alpha+beta/2.d0)**(-1) * (alpha*qnew2(:,i0,j0)+&
                 beta/2.d0*qnew2(:,i_n,j_n))
            qnew2(:,i0,j0) = nhood_val
            qnew2(:,i_n,j_n) = qnew2(:,i_n,j_n) + nhood_val/2.d0
          end if
        end if
      end do


    end subroutine

    subroutine rotate_state(q,q_rot,n_vec,t_vec)
      ! n_vec is the normal direction unit vector
      ! t_vec is the transverse direction unit vector, OG to n_vec
      ! q is the Cartesian coordinate aligned state vec
      ! q_rot is the rotated state vec
      implicit none
      real(8) :: q(3),q_rot(3),n_vec(2),t_vec(2)
      real(8) :: vel(2)

      q_rot(1) = q(1)
      vel = q(2:3)
      q_rot(2) = vel(1)*n_vec(1) + vel(2)*n_vec(2)
      q_rot(3) = vel(1)*t_vec(1) + vel(2)*t_vec(2)
    end subroutine










end module
