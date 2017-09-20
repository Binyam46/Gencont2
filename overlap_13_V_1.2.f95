module opt_iter_mod

use access_mod

contains

!*****************************************************************************************************
!
!
!*****************************************************************************************************
Subroutine oc_iter (locate, c,c_set, K2, ebv, a_c, Q, A, niter, n_run, ones_1, ncand, nsex, & 
			xxx, Error_c, thold, round, restrict, islow)
  
	Implicit none

	integer, intent(in)						:: ncand, nsex
	integer, intent(inout)						:: xxx, round, Error_c
    	integer, intent(in)						:: niter
	integer, intent(inout)						:: n_run
	integer, intent(in)						:: restrict(2), islow
    	integer, dimension(ncand), intent (in)				:: locate
	real(8), dimension(nsex), intent (in)				:: ones_1
    	real(8), dimension(ncand), intent(out)				:: c, c_set
 	real(8), dimension(ncand, nsex), intent(in)			:: Q
    	real(8), dimension(ncand), intent(in)				:: ebv,  a_c
    	real(8), dimension(ncand, ncand), intent(in)			:: A
    	real(8), intent(in)						:: thold, k2


       	integer								:: i, j, ErrorFlag, x, in_iter, y, x_alt,nsex_1
	integer								:: ii, jj, mini 
	integer, dimension(nsex)					:: rem_animals, r_indx
	integer(4) now(3)
	integer, dimension(ncand, ncand)				:: ID
       	real(8), dimension(xxx, xxx)					:: P
       	real(8)								:: apa, rqr, aqr, DOM, NOM
	integer, dimension(ncand)					:: locate_in, locate_alt
	integer, dimension(xxx)						:: ac_index
       	real(8)								:: dlamda0, deltk, ggain, cAc, lamda0
       	real(8) 							:: deltK_old , sse_1, det, sst, threshold
	real(8), dimension (:,:), allocatable 				:: Q_n, tQ, QQ, QQII, QQII_1
       	real(8), dimension (:), allocatable 				:: Qlamda, ebv_n, lamda, aac,sse, iaac,sse_2
       	real(8), dimension (:), allocatable				:: temp_c, deltC, startsol, Q3, Q2, startsol_a
       	real(8), dimension (:), allocatable 				:: ACCA, RR, deltlamda,delr, r_0
      
       	allocate (ebv_n(xxx),aac(xxx), temp_c(ncand), startsol(xxx))
      	allocate (ACCA(xxx), RR(nsex), deltC(xxx))
	allocate (Qlamda(xxx), iaac(xxx))
	allocate (Q3(xxx), Q2(xxx), startsol_a(xxx))

    	threshold=0.0; ebv_n=0.
	rem_animals=0; nsex_1=0; mini=0; aac=0.
    	locate_in=locate
	ebv_n=ebv(locate_in(1:xxx))
	aac=a_c(locate_in(1:xxx))
	y=xxx
!print*, locate_in(1:xxx)
!print*, 'xxx =', xxx
!stop
!print*, 'the ebv is ', ebv_n
!print*,  'the aac is ', aac
!print'(30i4)',Q
!print*,'before nsex =', nsex
!stop
!print*,'row =',size(Q,1)
!print*,'column',count(Q(locate_in(1:xxx), 1)==1.)

! Check remaining individuals in the run	

	do i=1,nsex
		rem_animals(i) = count(Q(locate_in(1:xxx), i)==1.)

		if( (rem_animals(i) .gt. 0) ) then ! .and. (ones_1(i) .gt. 1E-6)) then
			nsex_1=nsex_1+1
		end if
	end do

!print*, xxx
!print*,'rem_animals =', rem_animals
!print*,'after nsex =', nsex_1
!stop

	allocate (Q_n(xxx, nsex_1), tQ(nsex_1, xxx), QQ(nsex_1, nsex_1), QQII(nsex_1, nsex_1),deltlamda(nsex_1))
	allocate (delr(nsex_1),lamda(nsex_1), r_0(nsex_1), sse(nsex_1), QQII_1(nsex_1,nsex_1), sse_2(3))

Q_n=0.; QQ=0.; tQ=0.; QQII=0.; QQII_1=0.

! Check and correct the ones and Q
	j=0
	do i=1,nsex
		
		if((rem_animals(i) > 0) ) then ! .and. (ones_1(i) .gt. 1E-5)) then
			j=j+1
			Q_n(:,j)=Q(1:xxx, i)
			r_0(j) = ones_1(i)
		end if
	end do

!print*, 'inside nsex is', nsex_1	
!print*, 'inside r is', r_0
!print*, 'inside k2 = ', k2
!print*,'Here 2'
!print*,'Here 2',Q_n
now=0
!stop

!===== Calculation of lamda0 =======!

    	tQ =transpose(Q_n)
	QQ= matmul(tQ, Q_n)
    
	x=size(QQ,1)
    	if(x /=2) then
    		CALL FindInv(QQ, QQII, x, ErrorFlag)
    	else if (QQ(1,1)*QQ(2,2) > 0) then 
      		det=QQ(1,1)*QQ(2,2) - QQ(1,2)**2
       		QQ=reshape((/QQ(2,2), -QQ(1,2), -QQ(1,2), QQ(1,1)/), (/2,2/))
        	QQII=QQ/det
    	else if(QQ(1,1) > 0) then
        	QQII=0
        	QQII(1,1)=1./QQ(1,1)
    	else if(QQ(2,2) > 0) then
      		QQII=0
        	QQII(2,2)=1./QQ(2,2)
    	else 
     	 	print*,'ERROR : QQ=0'
    	end if
!print*,size(QQII,1)
!print*,size(QQII,2)
!print*,size(QQII_1,1)
!print*,size(QQII_1,2)

	QQII_1=QQII
    	P=matmul(matmul(Q_n,QQII),tQ)
    	Q3=matmul(P, ebv_n)
    	Q2= matmul(P, aac)


	NOM= dot_product(ebv_n, ebv_n) - dot_product(ebv_n,Q3)
    	apa= dot_product(aac, aac) - dot_product(aac,Q2)
	!apa=dot_product(aac,Q2)
    	rqr= dot_product(r_0, matmul(QQII, r_0))
    	aqr = dot_product(aac,matmul(matmul(Q_n,QQII),r_0))
	
	!print*, 'EBV', ebv_n,  dot_product(ebv_n,Q3)
	!print*, 'NOM', NOM
	!print*, 'apa', apa
	!print*, 'rqr', rqr	

    	DOM= k2 + apa - rqr -2*(aqr)
!print*, 'DOM', DOM

!*** Minimisation option ******

	if(DOM <= 0 .OR. NOM <= 0 .OR. all(ebv_n==0.0)) then

		mini=1
		lamda0=10
	else
		lamda0=sqrt(NOM/(16*DOM))
	end if


!	if (DOM <= 0) then
! 			stop 'Error constraint cannot be achieved'
!	end if
        
!	lamda0=sqrt(NOM/(16*DOM))


!print*, 'mini', mini	
!print*, 'lamda0', lamda0
!stop
 
! ==== The iteration begins here ====
 
	lamda=0
	deltK=0
	round=0
	c=0
	c_set=0
 
now=0
!open(54, file='lamda0.out')
!open(64, file='dlamda.out')
!open(32, file='sse.out')
!open(74, file='CAC.out')

do i=1, niter  !! For main iteration
  	round=round+1
	locate_alt(1:xxx)=locate_in(1:xxx)
	x_alt = xxx

!if (round ==150) then
! CALL itime(now)
!print*, c
!stop
!end if
   !print*, 'round', round
!if(mod(round,100)==0) print*,'round', round
	!temp_c=c_set 

  	!if(islow==1) n_run=niter
    	if (round .GT. n_run+100) then

			x=size(c(1:xxx),1)
  			xxx=0
			n_run = n_run +100
   			do j= 1, x
 				if(c(j) .GT. 0.0) then
   			 		xxx=xxx+1
    					locate_in(xxx)=locate_in(j)
				end if

			end do
			
			if(xxx==0) then
				xxx=x_alt
				locate_in(1:xxx) = locate_alt(1:xxx)
			end if

		ebv_n(1:xxx)=ebv(locate_in(1:xxx))
		aac(1:xxx)=a_c(locate_in(1:xxx))
		
		deallocate (Q_n, tQ, QQ, QQII_1)
       	allocate (Q_n(xxx, nsex_1), tQ(nsex_1, xxx), QQ(nsex_1, nsex_1), QQII_1(nsex_1, nsex_1))
	
		jj=0
		do ii=1, nsex
			if((rem_animals(ii) .GT. 0)) then ! .and. (ones_1(ii) .gt. 1E-5)) then
				jj=jj+1
				Q_n(:,jj)=Q(locate_in(1:xxx), ii)
			end if
		end do

	!	if(rem_animals(1) == 0 .and. rem_animals(2) .GT. 0) then
	!		Q_n(:,1)=Q(locate_in(1:xxx), 2)
	!		!QQII_1=QQII(1,1)!*(y/xxx)
	!	else if (rem_animals(1) .GT. 0 .and. rem_animals(2) == 0) then
	!		Q_n(:,1)=Q(locate_in(1:xxx), 1)
	!		!QQII_1=QQII(2,2)!*(y/xxx)
	!	else 
	!		Q_n=Q(locate_in(1:xxx), :)
	!		!QQII_1=QQII!*(y/xxx)
	!	end if


!write(74, *) round, xxx, sum(Q_n(:,1)) 
		
    		tQ =transpose(Q_n)
		QQ= matmul(tQ, Q_n)
		
		x=size(QQ,1)
    		if(x /=2) then
    			CALL FindInv(QQ, QQII, x, ErrorFlag)
    		else if (QQ(1,1)*QQ(2,2) > 0) then 
      			det=QQ(1,1)*QQ(2,2) - QQ(1,2)**2
       			QQ=reshape((/QQ(2,2), -QQ(1,2), -QQ(1,2), QQ(1,1)/), (/2,2/))
        		QQII_1=QQ/det
    		else if(QQ(1,1) > 0) then
        		QQII_1=0
        		QQII_1(1,1)=1./QQ(1,1)
    		else if(QQ(2,2) > 0) then
      			QQII_1=0
        		QQII_1(2,2)=1./QQ(2,2)
    		else 
     	 		print*,'ERROR : QQ=0'
    		end if
!do ii=1, xxx
!	print*, xxx, c_set(locate_in(ii)), c(ii)
!end do 
!stop

	end if  

!print*, 'lamda =', lamda
	!lamda = matmul(QQII,matmul(tQ,ebv_n(1:xxx) - lamda0*aac(1:xxx))- r_0)
    	!lamda=matmul(QQII,(matmul(tQ,(ebv_n-(2*lamda0*aac)))-(2*lamda0*r_0)))

!if(round > 499) then
	!print*, size(ACCA,1)
	!print*, xxx
!end if
	if (mini== 1) then
		ACCA(1:xxx) = matmul(Q_n,lamda)/(2*lamda0) + aac(1:xxx)
	else
		ACCA(1:xxx) = ((ebv_n(1:xxx) - matmul(Q_n,lamda))/(2*lamda0))- aac(1:xxx)
		!ACCA = (ebv_n(1:xxx) - (2*lamda0*aac(1:xxx))-matmul(Q_n,lamda))/(2*lamda0)
		!ACCA = ebv_n(1:xxx)/(lamda0) - aac(1:xxx)-matmul(Q_n,lamda)/(2*lamda0)
	end if

	startsol(1:xxx) = c_set(locate_in(1:xxx)) !c(1:xxx) !0
	temp_c(1:xxx) = c_set(locate_in(1:xxx)) !c(1:xxx)
	x=size(c(1:xxx),1)

!CALL itime(now)
!print*, 'matmul(Q_n,lamda)', matmul(Q_n,lamda)

!print*, 'ACCA =', ACCA
!stop

   
	CALL solv (A(locate_in(1:xxx),locate_in(1:xxx)), ACCA(1:xxx), startsol(1:xxx), x, in_iter, restrict)
   	c(1:xxx)=startsol(1:xxx)

!CALL itime(now)
!print*, now

!if (round ==2) then
!stop
!end if
	c_set=0
	do j=1,xxx    
		c_set(locate_in(j))=startsol(j)
	end do

	cAc= dot_product(c(1:xxx), matmul(A(locate_in(1:xxx),locate_in(1:xxx)),c(1:xxx) )) &
			+ 2*dot_product(c(1:xxx), aac(1:xxx))
	
	if (mini ==1) then

	    	delr = (r_0 - (matmul(tQ,c(1:xxx))))
       		!if(any(delr >= 2.0)) then
      		!	delr=1.0
    		!end if
    		deltlamda = (matmul(QQII_1,delr))*(2*lamda0)
	else
	    	delr = -(r_0 - (matmul(tQ,c(1:xxx))))
	    	!delr = (r_0 - (matmul(tQ,c(1:xxx))))

       		!if(any(delr >= 2.0)) then
      		!	delr=1.0
    		!end if
    		deltlamda = (matmul(QQII_1,delr))*(2*lamda0)
	end if

	lamda = lamda + deltlamda
   
       	!if(any(deltlamda == 0.0)) then
      	!	deltlamda=1.0
    	!end if

	
!print*, 'lamda', lamda
!stop
    
	startsol_a(1:xxx) = 0 !iaac(1:xxx)
	!iaac=0
	x=size(aac(1:xxx),1)
  
	CALL solv (A(locate_in(1:xxx),locate_in(1:xxx)), aac(1:xxx), startsol_a(1:xxx), x, in_iter, restrict)
!   	iaac(1:xxx)=startsol_a(1:xxx)
!print*, 'aac =', aac
!print '(30F4.1)', iaac
!print*, A(locate_in(1:xxx),locate_in(1:xxx))
!stop


	deltK_old = deltK

	Qlamda (1:xxx)= matmul(Q_n,lamda)
    	deltK = K2 - dot_product(c(1:xxx),matmul(A(locate_in(1:xxx),locate_in(1:xxx)),c(1:xxx))) &
			- 2*dot_product(c(1:xxx), aac(1:xxx))

	if (mini == 1) then
	      	dlamda0 = -((lamda0**2)*deltK)/dot_product((c(1:xxx) + iaac(1:xxx)), Qlamda(1:xxx) )
	else
	      	!dlamda0 = -((lamda0**2)*deltK)/dot_product((c(1:xxx) + iaac(1:xxx)), ebv_n(1:xxx) - Qlamda(1:xxx) )
		dlamda0 = ((lamda0**2)*deltK)/dot_product((aac(1:xxx) - c(1:xxx)), ebv_n(1:xxx) - Qlamda(1:xxx) - 2*lamda0*aac(1:xxx))

	end if

	lamda = lamda + deltlamda
  	lamda0=lamda0 + dlamda0

       if(lamda0 <= 0.0) then
      		lamda0=0.0001
	!elseif (lamda0 > 2000) then
	!	lamda0= 1999
    	end if
	
!write(54,*) lamda0
!write(64,*) deltlamda
!write(74,54) cAc
!54 format(f12.7)


       
	sse_1=dot_product(temp_c (1:xxx)-c(1:xxx), temp_c(1:xxx)-c(1:xxx))
	sst = dot_product(c(1:xxx),c(1:xxx))

      	if (sst .NE. 0) then
               	sse_2(1) = sse_1/sst
	else
		sse_2(1) = 1.0
      	end if 

	sse_2(2)=deltk !- deltK_old
	sse_2(3)=maxval(delr)

!if(round==niter) print*, sse_2

!write(32,222) sse_2(1), sse_2(2), sse_2(3)
!222 format(3f12.7)  
       
     	sse_2=abs(sse_2)

	!if((restrict(1) .GT. 0) .OR. (restrict(2) .GT. 0)) then

if(mini==1) then

        if ((sse_2(1) .LT. thold) .and. (sse_2(2) .LT. 1E-3) .and. (sse_2(3) .LT. 1E-4) ) then
      
          	Error_c=0
          	!print*, 'Converged after', round, ' itrations'
          	exit 
           
        else if (niter == round .AND. ((sse_2(1) .GE. thold) .OR. (sse_2(2) .LT. 1E-4) .OR. (sse_2(3) .LT. 1E-4))) then
          	Error_c=1
          	!print*, 'The solutions after',round, ' itrations are provided in the output file' 
          	!print*
        end if 
else         
	if ((sse_2(1) .LT. thold) .and. (sse_2(2) .LT. 1E-4) .and. (sse_2(3) .LT. 1E-3) ) then
  	!if ((sse_2(1) .LT. thold)) then
          	Error_c=0
          	!print*, 'Converged after', round, ' itrations'
          	exit 
           
        else if (niter == round .AND. ((sse_2(1) .GE. thold) .OR. (sse_2(2) .LT. 1E-4) .OR. (sse_2(3) .LT. 1E-3))) then
          	Error_c=1
          	!print*, 'The solutions after',round, ' itrations are provided in the output file' 
          	!print*
        end if 
end if

if(niter < 10001) Error_c=0
end do	!!!! For main iteration

!close(54)
!close(64)
!close(32)
!close(74)
	c_set=0
	do i=1,xxx    
		c_set(locate_in(i))=startsol(i)
!		if(c_set(locate_in(i))==0.000001) c_set(locate_in(i))=0
	end do
	
	c(1:xxx)=startsol(1:xxx)
!print*, 'xxx is =', xxx
!print'(17i6)', locate(1:xxx)
!print'(8F6.3)', startsol(1:xxx)

!print*, locate_in(1:5)
!print*
!print*, 'the rema r is ', r_0 - matmul(tQ,c(1:xxx))
!print*, 'delta k is', deltk
!print*, 'round =', round
!print*, 'delr', delr
!print*, 'the c*ac is', dot_product(c(1:xxx),matmul(A(locate_in(1:xxx),locate_in(1:xxx)),c(1:xxx)))  &
!	+ 2*dot_product(c(1:xxx), aac(1:xxx)) 
!print'(A10, F6.3)', 'sum(c) is', sum(c(1:xxx))
!print*
end subroutine oc_iter

end module  opt_iter_mod
