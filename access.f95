module access_mod

!this are a collection of accossory functions for the main prrogram

contains

!==============================================================
! This subroutine solve equations with out involving inverting 
!	of large size matrices
!===============================================================
subroutine solv (Amat, rhs, startsol, n, in_iter, restrict)
    
implicit none

       	integer, intent(in)				:: n
	integer, intent(out)				:: in_iter
	integer, dimension(2), intent(in)		:: restrict
       	real(8), dimension (n,n), intent (in) 		:: Amat
       	real(8), dimension (n), intent (inout)		:: startsol
       	real(8), dimension (n), intent (in)		:: rhs

       	integer						:: iter,i, n_iter
       	real(8)						:: sse, sst, dtmp, tmp, tmp_value
       	real(8), dimension (n)				:: sol

	iter = 0
       	sse=1
       	sol = startsol
	!startsol = 0

	if((restrict(1) .NE. 0) .OR. (restrict(2) .NE. 0)) then
		n_iter=2000
		tmp_value=0
	else
		n_iter=2000
		tmp_value=0
	end if

	do while (sse .GT. 1E-3 .and. iter .LE. 2000)
        	iter=iter+1
          	sse=0
          	sst=0
          	startsol= sol

          	do i=1,n
            		dtmp=(rhs(i) - dot_product(Amat(i,:), sol)) / (Amat(i,i))
              		tmp=sol(i) + dtmp

                	if(tmp .GT. 0.) then
                		sol(i)=tmp
                	else
                		sol(i)= 0
                	end if
               	
			sse=sse + (startsol(i) - sol(i))**2
             		sst=sst + sol(i)**2
           	end do

           	if (sst .GT. 0) then
               		sse=sse/sst
           	end if 
        end do

	in_iter = iter
       	startsol=sol
         
	return

    end subroutine solv

!==============================================================
! End of subroutine 'solv'
!==============================================================

!==============================================================
! This subroutine solve equations with out involving inverting 
!	of large size matrices
!===============================================================
 	subroutine solve (Amat, rhs, startsol, n, in_iter)
    
    	implicit none

        integer, intent(in)					:: n
	integer, intent(out)					:: in_iter
        real(8), dimension (n,n), intent (in) 			:: Amat
        real(8), dimension (n), 	intent (inout)		:: startsol
        real(8), dimension (n), intent (in)			:: rhs

        integer							:: iter,i, yy, y, j, keep
        integer, dimension(n)					:: look
        real(8)							:: sse, sst, dtmp, tmp
        real(8), dimension (n)					:: sol

		look =0
			do i= 1, n
  				look(i)=i
			end do 

		iter = 0
        !thold1= 1E-6
        sse=1
        sol = startsol
		!yy=n
        keep=n
        
       	do while (sse .GT. 1E-6 .and. iter .LE. 2000)
          	iter=iter+1
	       	sse=0
          	sst=0
          	yy=keep
          	y=size(sol(look(1:yy)), 1)
                               
            startsol(1:yy)= sol(look(1:yy))     
        	sol(1:yy)= sol(look(1:yy))
            keep=0
            
 			do i=1,y
            	dtmp=(rhs(look(i)) - dot_product(Amat(look(i),look(1:yy)), sol(1:yy))) / (Amat(look(i),look(i)))
                tmp=sol(i) + dtmp
				
				if(sol(i) .GE. -0.01) then
   					keep = keep +1
    				look(keep)=look(i)
   				end if
                        
                if(tmp .GT. 0) then
                  sol(i)=tmp
                else
                  sol(i)= 0        
                end if
               	
				sse=sse + (startsol(i) - sol(i))**2
                sst=sst + sol(i)**2
            end do

             if (sst .GT. 0) then
               sse=sse/sst
             end if 
        end do 
         
		in_iter = iter
!         startsol=sol
        startsol = 0
        sol(1:yy)= sol(look(1:yy))
         
		do i=1,yy    
			startsol(look(i))= sol(i)
		end do
         
		return
    end subroutine solve

!==============================================================
! End of subroutine 'solve'
!==============================================================

!============================================================
!	This subroutine calcualtes relationship matrix 
!
!============================================================
subroutine A_mat (Acand, ncand, cand, n, ped_file)
implicit none

! define input parameters for the subroutines
	INTEGER, INTENT(IN)						::n, ncand
	real(8), dimension (ncand, ncand), intent(out)		:: Acand
	integer, dimension (ncand), intent(in) 			:: cand
	character(len=100), intent(in) 				:: ped_file

! define inside sburoutine parameters
	real(8), DIMENSION(:,:), allocatable			:: A
	INTEGER 							:: i, j,kd, ks,l, m
	integer 							:: iid, sid, did 
	integer, dimension (2,n)		 			:: ped


! read a pedigree file
	l=-1
	open(12,file=ped_file,status='old')
   	
		do 
			l=l+1
			if(l > n)then
				stop " ERR: NPED too small in the parameter file"
			end if

 			read(12,*,end=9)iid,sid,did
	 		ped(1:2,iid)=(/sid,did/)
		end do

		9 continue
	close(12)
!print*, 'l =',l
!stop
	m=l
	allocate(A(m,m))

! get an identity matrix
	do i=1,m
   		A(i,i)=1
	end do  

! the program begin here
	do i=1,m
    	do j= i,m
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
    
end subroutine A_mat
!=============================================================


!===========================================================
! This subroutine reads in full-saved relation matrix, A, from a file 
! 
!===========================================================

subroutine readin_s (A, ncand, a_file_name)

	implicit none

	integer, intent(in) :: ncand
	real(8), dimension (ncand,ncand) , intent(out)		:: A
	character(len=100), intent(in) 				:: a_file_name

	integer :: i,j

	open(10, file=a_file_name, STATUS='OLD')
		do 11 i=1,ncand
    		read(10,*)(A(i,j),j=1,ncand)
   		11 continue
	close(10)

end subroutine readin_s

!***********************************************************
! end of subroutine 'readin 
!************************************************************

!===========================================================
! This subroutine reads in half-saved relation matrix, A, from a file 
! 
!===========================================================

subroutine readin_half (A, ncand, a_file_name)

	implicit none

	integer, intent(in) :: ncand
	real(8), dimension (ncand,ncand) , intent(out)		:: A
	character(len=100), intent(in) 				:: a_file_name

	real(8), dimension(ncand*(ncand+1)/2)			:: A_half
	integer :: i,j

	j=ncand*(ncand+1)/2
	open(12, file=a_file_name, STATUS='OLD')
		do 11 i=1,j
    		read(12,*)A_half(i)
   		11 continue
	close(12)

	CALL half2full( A_half, A, ncand)

end subroutine readin_half

!***********************************************************
! end of subroutine readin_half 
!************************************************************

! ====================================
  subroutine half2full( ww, hh, mfit)
! ====================================

  implicit none
  integer, intent(in)                              :: mfit
  real(8), dimension(mfit, mfit), intent(out)      :: hh
  real(8), dimension(mfit*(mfit+1)/2), intent(in)  :: ww
  integer                                          :: i, j, ij

  ij = 0
  do i = 1, mfit
     do j = i, mfit
        ij = ij + 1; hh(j,i) = ww(ij); hh(i,j) = ww(ij)
     end do
  end do

  return
  end subroutine half2full
!***************************************************
! End of subroutine half2full
!****************************************************

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

!============================================================
!	A subroutine to for identity matrix
!============================================================

	subroutine idnt(m, mat)

		implicit none

		integer, intent (IN)					:: m
 		integer, dimension (m,m), intent(out)		:: mat
        	integer						:: i

		mat=0

		do i=1, m
			mat(i,i)= 1
		end do
			

	end subroutine idnt
!============================================================

!============================================================
!	A function to calculate a norm of a vector 
!============================================================

	function norm(vec, m)

		implicit none

		real(8)							:: norm, a, b
		integer, intent (IN)					:: m
 		real(8), dimension (:), intent(IN)			:: vec
        	integer						:: i

			a=0
			b=0

			do i=1,m
			  a= vec(i)**2
			  b= b + a
			end do

			norm=b

	end function norm
!============================================================

!=============================================================
! this subroutine calculates inverse of a matrix
!	(not my subroutine)
!=============================================================

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
        IMPLICIT NONE
        !Declarations
        INTEGER, INTENT(IN) 			:: n
        INTEGER, INTENT(OUT) 			:: errorflag  !Return error status. -1 for error, 0 for normal
        real(8), INTENT(IN), DIMENSION(n,n) 	:: matrix  !Input matrix
        real(8), INTENT(OUT), DIMENSION(n,n) 	:: inverse !Inverted matrix
       
        LOGICAL 					:: FLAG = .TRUE.
        INTEGER 					:: i, j, k, l
        real(8) 					:: m
        real(8), DIMENSION(n,2*n) 			:: augmatrix !augmented matrix
       
        !Augment input matrix with an identity matrix
        DO i = 1, n
                DO j = 1, 2*n
                        IF (j <= n ) THEN
                                augmatrix(i,j) = matrix(i,j)
                        ELSE IF ((i+n) == j) THEN
                                augmatrix(i,j) = 1
                        Else
                                augmatrix(i,j) = 0
                        ENDIF
                END DO
        END DO
       
        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
                IF (augmatrix(k,k) == 0) THEN
                        FLAG = .FALSE.
                        DO i = k+1, n
                                IF (augmatrix(i,k) /= 0) THEN
                                        DO j = 1,2*n
                                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                                        END DO
                                        FLAG = .TRUE.
                                        EXIT
                                ENDIF
                                IF (FLAG .EQV. .FALSE.) THEN
                                        PRINT*, "Matrix is non - invertible"
                                        inverse = 0
                                        errorflag = -1
                                        return
                                ENDIF
                        END DO
                ENDIF
                DO j = k+1, n                       
                        m = augmatrix(j,k)/augmatrix(k,k)
                        DO i = k, 2*n
                                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
                        END DO
                END DO
        END DO
       
        !Test for invertibility
        DO i = 1, n
                IF (augmatrix(i,i) == 0) THEN
                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO
       
        !Make diagonal elements as 1
        DO i = 1 , n
                m = augmatrix(i,i)
                DO j = i , (2 * n)                               
                           augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO
       
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
                m = augmatrix(i,k+1)
                        DO j = k, (2*n)
                                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                        END DO
                END DO
        END DO                               
       
        !store answer
        DO i =1, n
                DO j = 1, n
                        inverse(i,j) = augmatrix(i,j+n)
                END DO
        END DO
        errorflag = 0
END SUBROUTINE FINDinv
!==============================================================================


!=============================================================================
!This function tokenize strings into single words. It is not writen by me
!
!=============================================================================
FUNCTION strtok (source_string, delimiters)

IMPLICIT NONE

!     @(#) Tokenize a string in a similar manner to C routine strtok(3c). 
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in DELIMITERS.
!
!             then, if the returned value is not equal to CHAR(0), keep calling until it is
!             with SOURCE_STRING set to CHAR(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning CHAR(0). 
!
!     Input:  source_string =   Source string to tokenize. 
!             delimiters    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!     Output: strtok()
!
!     LIMITATIONS:
!     can not be called with a different string until current string is totally processed, even from different procedures
!     input string length limited to set size
!     function returns fixed 255 character length
!     length of returned string not given

!     PARAMETERS:
      CHARACTER(len=*),intent(in)  :: source_string
      CHARACTER(len=*),intent(in)  :: delimiters
      CHARACTER(LEN=255) :: strtok

!     SAVED VALUES:
      CHARACTER(len=255),save :: saved_string
      INTEGER,save :: isaved_start  ! points to beginning of unprocessed data
      INTEGER,save :: isource_len   ! length of original input string

!     LOCAL VALUES:
      INTEGER :: ibegin        ! beginning of token to return
      INTEGER :: ifinish       ! end of token to return

      ! initialize stored copy of input string and pointer into input string on first call
      IF (source_string(1:1) .NE. CHAR(0)) THEN
          isaved_start = 1                 ! beginning of unprocessed data
          saved_string = source_string     ! save input string from first call in series
          isource_len = LEN(saved_string)  ! length of input string from first call
      ENDIF

      ibegin = isaved_start

      DO
         IF ( (ibegin .LE. isource_len) .AND. (INDEX(delimiters,saved_string(ibegin:ibegin)) .NE. 0)) THEN
             ibegin = ibegin + 1
         ELSE
             EXIT
         ENDIF
      ENDDO

      IF (ibegin .GT. isource_len) THEN
          strtok = CHAR(0)
          RETURN
      ENDIF

      ifinish = ibegin

      DO
         IF ((ifinish .LE. isource_len) .AND.  (INDEX(delimiters,saved_string(ifinish:ifinish)) .EQ. 0)) THEN
             ifinish = ifinish + 1
         ELSE
             EXIT
         ENDIF
      ENDDO

      !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
      strtok = saved_string(ibegin:ifinish-1)
      isaved_start = ifinish

END FUNCTION strtok
!=================================================================================================

!================================================================================================
!
!==============================================================================================
SUBROUTINE DELIM(LINE0,ARRAY,N,ICOUNT,IBEGIN,ITERM,ILEN,DLIM)


IMPLICIT NONE
!     @(#) parse a string and store tokens into an array
!
!     given a line of structure " par1 par2 par3 ... parn "
!     store each par(n) into a separate variable in array.
!
!     IF ARRAY(1).eq.'#NULL#' do not store into string array  (KLUDGE))
!
!     also count number of elements of array initialized, and
!     return beginning and ending positions for each element.
!     also return position of last non-blank character (even if more
!     than n elements were found).
!
!     no quoting of delimiter is allowed
!     no checking for more than n parameters, if any more they are ignored
!
!     input line limited to 1024 characters
!
      CHARACTER*(*)   ::  LINE0, DLIM*(*)
      INTEGER, PARAMETER ::MAXLEN = 1024
      CHARACTER*(MAXLEN) :: LINE
      INTEGER	:: N, IFOUND, I10, IEND,ISTART, IARRAY, ICOL, IDLIM!, MAXLEN
      CHARACTER ::ARRAY(N)*(*)
      INTEGER ::ICOUNT, IBEGIN(N),ITERM(N),ILEN
      LOGICAL LSTORE
      ICOUNT=0
      ILEN=LEN_TRIM(LINE0)
      IF(ILEN.GT.MAXLEN)THEN
         write(*,*)'*delim* input line too long'
      ENDIF
      LINE=LINE0

      IDLIM=LEN(DLIM)
      IF(IDLIM.GT.5)THEN
!        dlim a lot of blanks on some machines if dlim is a big string
         IDLIM=LEN_TRIM(DLIM)
!        blank string
         IF(IDLIM.EQ.0)IDLIM=1
      ENDIF

!     command was totally blank
      IF(ILEN.EQ.0)RETURN
!
!     there is at least one non-blank character in the command
!     ilen is the column position of the last non-blank character
!     find next non-delimiter
      icol=1

!     special flag to not store into character array
      IF(ARRAY(1).EQ.'#NULL#')THEN
         LSTORE=.FALSE.
      ELSE
         LSTORE=.TRUE.
      ENDIF

!     store into each array element until done or too many words
      DO 100 IARRAY=1,N,1
200      CONTINUE
!        if current character is not a delimiter
         IF(INDEX(DLIM(1:IDLIM),LINE(ICOL:ICOL)).EQ.0)THEN
!          start new token on the non-delimiter character
           ISTART=ICOL
           IBEGIN(IARRAY)=ICOL
!          assume no delimiters so put past end of line
           IEND=ILEN-ISTART+1+1

           DO 10 I10=1,IDLIM
              IFOUND=INDEX(LINE(ISTART:ILEN),DLIM(I10:I10))
              IF(IFOUND.GT.0)THEN
                IEND=MIN(IEND,IFOUND)
              ENDIF
10         CONTINUE

!          no remaining delimiters
           IF(IEND.LE.0)THEN
             ITERM(IARRAY)=ILEN
             IF(LSTORE)ARRAY(IARRAY)=LINE(ISTART:ILEN)
             ICOUNT=IARRAY
             RETURN
           ELSE
             IEND=IEND+ISTART-2
             ITERM(IARRAY)=IEND
             IF(LSTORE)ARRAY(IARRAY)=LINE(ISTART:IEND)
           ENDIF
           ICOL=IEND+2
         ELSE
           ICOL=ICOL+1
           GOTO 200
         ENDIF
!        last character in line was a delimiter, so no text left
!        (should not happen where blank=delimiter)
         IF(ICOL.GT.ILEN)THEN
           ICOUNT=IARRAY
           RETURN
         ENDIF
100   CONTINUE
!     more than n elements
      ICOUNT=N
      RETURN
      END SUBROUTINE
 !==========================================================================

 
      FUNCTION TRUNCATE (A)
      IMPLICIT NONE

      real(8) :: A
      real(8) ::TRUNCATE

          TRUNCATE=FLOAT (INT(A + 0.5))

       END  FUNCTION TRUNCATE

end module
