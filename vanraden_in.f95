module vanrnden_mod

implicit none

contains

subroutine vanraden_s (nped, ncand, Acand, cand, ped_file)

	!implicit none

! define input parameters for the subroutines
    		integer, intent(in)				:: nped, ncand
        	integer, dimension (ncand), intent(in) 		:: cand
        	real(8), dimension (ncand, ncand), intent(out)	:: Acand
       		character(len=100), intent(in) 			:: ped_file

! define inside sburoutine parameters
		integer 				::iid, sid, did, year
		integer					:: i, j, n, jj, k, jk, ks, kd,l
	    	integer, dimension (:,:), allocatable 	:: ped
    		integer, dimension (:), allocatable 	:: navg, ii
		real(8), dimension (:), allocatable 	:: f, avginb, avginbn
    		real(8) 				:: minb, minbold
		real(8), dimension (:,:), allocatable 	:: A 

	 	allocate (ped(3, nped)) !, A(nped, nped))

! Read pedigree file including birth year of individulas
 	
	l=-1
	open(12,file=ped_file,status='old')
	open(32, file ='Error.log',status="old",Access = 'append')
    	
		
		do 
			l=l+1
			if(l > nped)then
				stop " ERR: NPED too small in the parameter file"
			end if 

 			read(12,*,end=9)iid,sid,did, year
				if (year < 1000 .OR. year > 9999) then
					write(*,*) 'Birth year should be between 1000 and 9999'
					write(32,*) 'Birth year should be between 1000 and 9999. &
							Check the birth year values in the pedigre file'
					stop " ERR: Check the birth year value"
				end if

 				if (sid <=0 .and. did <=0) then
  					sid = - year
					did = - year
				elseif (sid > 0 .and. did <=0) then
					did = - year
				elseif (sid <= 0 .and. did > 0) then
					sid = - year
				end if
	 		ped(1:3,iid)=(/sid,did,year/)
		end do

	9 continue
	close(12)

! Detemine the maximum number of birth year of individulas in the pedigree file
	ii=maxloc(ped(3,:))
    	k = ped(3, ii(1))
      ! k = max(ped(3,:), 0)
    
     allocate (f(l), avginb(k), avginbn(k), navg(k),A(l, l))
         
! Calculate inbreeding coefficients and average inbreeding per year
	minbold = 0
	avginbn = 0

		do
			f = -1
	            	avginb = avginbn
            		avginbn = 0
            		navg = 0

			do i = 1, l
				if (f(i) == -1) f(i) = inbrec (i, f, k, avginb,ped, l)
				if (ped(1,i) > 0 .and. ped(2,i) > 0) then
					avginbn(ped(3,i)) = avginbn(ped(3,i)) + f(i)
					navg(ped(3,i)) = navg(ped(3,i)) + 1
				endif
			enddo
        
			minb=sum(f)/l
            
			where(navg/=0)
				avginbn=avginbn/navg
			end where
            
			if (abs(minbold-minb) <1.e-6) exit
            
			minbold=minb
    
		enddo

	avginb = avginbn

! Construct relationship matrix for candidate indivduals
A=0

	do i=1,l
    	do j= i,l
        	ks=ped(1,j)
            	kd=ped(2,j)
            
	! CASE one (both parents known)

            if (ks > 0 .and. kd > 0) then
              if(i==j) then 
                A(i,i)=1 + f(i)
                else 
                  A(i,j)=0.5*(A(i,ks)+ A(i,kd))
                  A(j,i)=A(i,j)
                end if
            end if

	!CASE 2 (sire is unkown)

            if (ks <= 0 .and. kd >0) then
              if (i==j) then
                A(i,i)=1 + f(i)
                else
                  A(j,i)=0.5*(A(i,kd) + 2*avginb(abs(ks)))
                  A(j,i)=A(i,j)
              end if
           end if

	! CASE 3 (dam is unknown)

           if (kd <=0 .and. ks >0) then
             if (i==j) then
               A(i,i)=1 + f(i)
               else
                 A(i,j)=0.5*(A(i,ks) + 2*avginb(abs(kd)))
                 A(j,i)=A(i,j)
               end if
            end if

	! CASE 4 (both sire and dam are unknown)

           if (kd <=0 .and. ks <= 0) then
             if (i==j) then
               A(i,i)=1 + f(i)
               else
                 A(i,j)=0.5*( (2*avginb(abs(kd))) + (2*avginb(abs(ks))))
                 A(j,i)=A(i,j)
               end if
            end if

         end do
       end do

Acand=A(cand, cand)

!Write out inbreeding coefficients and average inbreeding per year

	open(1003, file='inbcoef.out')
	open (1004, file='Aveinb.out')
 
		do i=1,l
   			write(1003, 103) i, f(i)
		103 format(i6, 2x, f10.8)
		end do

		do i=1,k
			write (1004, 104) avginb (i)
		104 format(f10.8)
		end do

	close (1003)
	close (1004)
	close (32)

end subroutine vanraden_s
!===================================================================================
! end of subroutine vanraden_s							=======
!===================================================================================


!========================================================================
!This function calculates inbreeding coefficient				===
!========================================================================

real function inbrec(an, f, k, avginb,ped, nped)
		! Returns inbreeding coefficient for animal = an
		! f(an) = 0.5 * cffa(s,d)
		! negative s or d means UPG code
        
	integer, intent(in)	:: nped, k
	integer, dimension (3,nped), intent(in):: ped
    	real(8), dimension (nped), intent(in) :: f
 	real(8), dimension (k), intent(in) :: avginb
	integer:: an, s, d
    
	s = ped(1,an); d = ped(2,an)
           !print*, s , d
		if (s <= 0 .or. d <= 0) then
     
			inbrec = avginb(abs(min(s,d)))
		else
        	inbrec = 0.5 * (cffa(s,d, f, k,avginb,ped, nped))
		endif
if(mod(an,100)==0) print*,'calculation of inbreeding at anim', an
        
end function inbrec
!========================================================================

!========================================================================
!This function calculates inbreeding coefficient				===
!========================================================================

recursive real(8) function cffa(a1,a2, f, k, avginb,ped, nped) result (rel)
		! Returns relationship between a1 and a2
    
    	integer, intent(in)	:: nped, k
   	integer, dimension (3,nped), intent(in):: ped
    	real(8), dimension (nped), intent(in) :: f
 	real(8), dimension (k), intent(in) :: avginb
	integer:: a1, a2
      
		if (a1 <= 0 .or. a2 <= 0) then
			rel = 2 * (avginb(abs(min(a1,a2))))
		elseif (a1 == a2) then
			rel = 1 + f(a1)
		else
			if (a1 < a2) then
				rel = 0.5 * (cffa(a1,ped(1,a2),f, k,avginb, ped, nped) + cffa(a1,ped(2,a2), f, k,avginb,ped, nped))
			else
				rel = 0.5 * (cffa(a2,ped(1,a1), f, k,avginb,ped, nped) + cffa(a2,ped(2,a1), f, k,avginb,ped, nped))
			endif
		endif
end function cffa

!========================================================================

end module vanrnden_mod
