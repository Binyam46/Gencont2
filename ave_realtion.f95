program ave_relation

!------------------------------------------------------------
!	This subroutine calcualtes relationship matrix 
!
!------------------------------------------------------------

implicit none

! define input parameters for the subroutines
	INTEGER					:: n, ncand,i, j,kd, ks
	integer 				:: iid, sid, did 
	real(8), dimension (:,:), allocatable 	:: Acand, A
	integer, dimension (:), allocatable	:: cand
	character(len=100)			:: ped_file, id_file, n_cand, n_ped, A_name
	integer, dimension (:,:), allocatable	:: ped
	real(8)					:: ave


      	CALL getarg(1, ped_file)
      	CALL getarg(2, id_file)
	CALL getarg(3, n_ped)
 	CALL getarg(4, n_cand)
	CALL getarg(5, A_name)

	read(n_ped,*) n
	read(n_cand,*) ncand
	
	allocate (Acand(ncand,ncand), A(n,n), cand(ncand), ped(2,n))


! read a pedigree file
	open(12,file=ped_file,status='old')
   	
		do 
 			read(12,*,end=9)iid,sid,did
	 		ped(1:2,iid)=(/sid,did/)
		end do

		9 continue
	close(12)

! read a ID file
	open(13,file=id_file,status='old')
   	
		do i=1,ncand
 			read(13,*,end=9)ks
	 		cand(i)=ks
		end do
!print*, cand
!stop
		
	close(13)

! get an identity matrix
	do i=1,n
   		A(i,i)=1
	end do  

! the program begin here
	do i=1,n
    	do j= i,n
        	ks=ped(1,j)
            	kd=ped(2,j)
            
            if (ks > 0 .and. kd > 0) then
              if(i==j) then 
                A(i,i)=1+0.5*A(ks,kd)
                else 
                  A(i,j)=0.5*(A(i,ks)+ A(i,kd))
                  A(j,i)=A(i,j)
                end if
            end if

            if (ks == 0 .and. kd >0) then
              if (i==j) then
                A(i,i)=1
                else
                  A(j,i)=0.5*A(i,kd)
                  A(j,i)=A(i,j)
              end if
           end if

           if (kd==0 .and. ks >0) then
             if (i==j) then
               A(i,i)=1
               else
                 A(i,j)=0.5*A(i,ks)
                 A(j,i)=A(i,j)
               end if
            end if
         end do
       end do
             
	Acand=A(cand, cand)


if(A_name .NE. '') then
	open(10, file=A_name)
		do i=1,ncand
			write(10, *) (Acand(i,j), j=1,ncand)
		end do
	close (10)
  		
end if 	



CALL aveR (Acand, ave, ncand)
print*, ave
    
contains

!***********************************************************
! This subroutine calculates average relationship 
!	the calculation don't include relationship of an animal to itself
!
!***********************************************************

	 subroutine aveR (Amat, ave, n)

		implicit none

		integer, intent (in)					:: n
		real(8), intent(inout)					:: ave
		real(8), dimension (n, n),intent(in) 			:: Amat
        	integer							:: i, j
        
		ave=0
		do i=1,(n-1)
  			do j=i+1, n
    				ave=ave + Amat(i,j)
!print*, Amat(i,j)
  			end do
		end do
		ave=ave/((n*(1 + n)/2) - n)

	end subroutine aveR
    
!============================================================
! End of subroutine 'aveR' 
!============================================================

end program ave_relation
!=============================================================