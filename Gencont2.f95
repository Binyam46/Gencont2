!=======================================================================================
!
! This program calculate optimum contricbution of animals 
! constraining on average relation among candidates	
!
!======================================================================================

include 'access.f95'
!include 'inbreeding_m.f95'
include 'vanraden_in.f95'
include 'overlap_13.f95'
!include 'inbreeding_m2.f95'
include 'relation.f95'

program optim_cont

use access_mod
use vanrnden_mod
use inbreeding_mod
!use inbreeding_mod2
use opt_iter_mod


Implicit none

	integer			:: i,j,l, Error_c, round,kk
	integer			:: ninp, inp(20), k_cand, ncand_, iout(20), nout
	integer			:: row
	real(8)			:: gen_gain, PaveA, covA
    	real(8)			:: PaveA_selct

	integer, dimension (:), allocatable 		:: sex, cand,locate,cand_group, loc_post
    	integer, dimension (:), allocatable 		:: anim, cand_id
    	character(20), dimension (:), allocatable 	:: name
	real(8), dimension (:,:), allocatable 		:: Q, Acand, A_young, A_anim
	real(8), dimension (:), allocatable 		:: aveA, c, ebv, c_young,w_young

type input
	character(20) 	:: name
    	integer		:: id
    	integer		:: sex
       	integer		:: age
       	integer		:: avail
    	real(8)		:: c_max
       	real(8)		:: c_min
       	real(8) 	:: c_prev
    	real(8)		:: ebv
	real(8)		:: avg_rel
end type input

type(input), dimension (:), allocatable		:: data_info
    
  	character(120)				:: str_in, line, line1
  	character (5) 				:: dlm =' '
  	integer, parameter 			:: fh = 15
  	integer 				:: ios = 0
  	integer 				:: line_1 = 0
  	Integer, parameter 			:: NT=13
	Integer 				:: Icount, Ilen
	Integer, dimension (NT)			:: Ibegin, Iterm
    	character(len=20)			:: str_out(NT)
    	integer(4) today(3), now(3)
    
! Control file variables
  	character(len=100) 			:: ped_file, data_file , Parafile,  outfile, a_file_name
  	integer 				:: nsex, ncand, nped, niter, nprog, minimize, nanim
    	integer					:: readin, vanraden, henderson, inbreedingm, overlap, a_half
	integer					:: restrict(2), select_all_females, n_females, multi_stage
  	integer					:: nages, islow, n_young, young_select, age_ch
	integer, dimension (20)			:: ngroup
	real(8) 				:: k, thold, deltaF,Mrelat, temp_val
	!real(8)				:: cmax_male, cmin_male, cmax_female, cmin_female, temp_val
	!real(8)				:: cmax_young, cmin_young
	real(8), dimension (:), allocatable	:: sum_cprev
	real(8), dimension (10)			:: cmax_group, cmin_group, s_group  

! Some varaibles for overlapping generation
	
	!character, dimension(:,:), allocatable		:: char_mat
	integer						:: iprev, jj, ii
	integer, dimension (:), allocatable 		:: list, ages, diff_ages
    	integer, dimension (:,:), allocatable 		:: nbar, Ntemp, n12
	real(8), dimension (:,:), allocatable 		:: abar, A12, Ctemp, char_mat
	real(8), dimension (:), allocatable 		:: CprevA12, ACprev, cbar, xbar 
	real(8), dimension (:), allocatable 		:: cmaxx, cminn, cmin, cmax, sum_cprev_age
	real(8)						:: CprevACprev, fact, aa

	integer						:: niter_w, i_w, i_age, i_sex, i_s
	real(8), dimension (:), allocatable 		:: w, avg_rel, cbar_1, A12c2, s, s_, w0
	real(8)						:: sumw, c2a22c2, cc
    	real(8)						:: o_sse, o_sst, avg_young
	real(8)						:: a11_young,a12_young
	integer						:: n11_young, n12_young


      CALL idate(today)
      CALL itime(now)
      CALL getarg(1, Parafile)
      CALL getarg(2, outfile)

	open(32, file ='Error.log',status="unknown")
	write(32, 1010 )  today(2), today(1), today(3)
    	write(32, 1011) now
 1010 format ( '				Date: ', i2.2, '/', i2.2, '/', i4.4)
 1011 format ( '				Time: ', i2.2, ':', i2.2, ':', i2.2)

	if(Parafile == '') then
		write(32, *) 'Gencont2 could not find a parameter file'
		write(32, *) 'Please refer to the user manual'
		stop 'Gencont2 could not find a parameter file'
	end if

	if(outfile == '') then
		outfile = 'gencont2.out'
	end if 

	open(31, file = outfile)

write(*,*)' .-----------------------------------------------------------------------------.'
write(*,*)'|                                                                               |'
write(*,*)'|    .--.    --   .    .   .--.   .--.   .    .   -----        .--.             |'
write(*,*)'|   :       |     |\   |  :      :    :  |\   |     |             :             |'
write(*,*)'|   |   _   :--   : \  :  |      |    |  : \  :     :     __    .-`             |'
write(*,*)'|   :   :   |     |  \ |  :      :    :  |  \ |     |          |                |' 
write(*,*)'|    `--`    --   :   \:   `--`   `--`   :   \:     :          :__.             |'
write(*,*)'|                                                                               |'
write(*,*)'|      (A program for calculation of Optimal Genetic Contributions)             |' 
write(*,*)'|                              Version 1.0                                      |'
write(*,*)'|                                                                               |'
write(*,*)'|                  Dept of Animal and Aguacultural Sciences                     |' 
write(*,*)'|                    Norwegian University of Life Sciences                      |'
write(*,*)'|                                    Norway                                     |'
write(*,*)'|                                                                               |'
write(*,*)'|                     Binyam Dagnachew and Theo Meuwissen                       |'
write(*,*)'|                         (Last update: October, 2015)                          |'           
write(*,*)' `-----------------------------------------------------------------------------`'


write(31,*)' .-----------------------------------------------------------------------------.'
write(31,*)'|                                                                               |'
write(31,*)'|    .--.    --   .    .   .--.   .--.   .    .   -----        .--.             |'
write(31,*)'|   :       |     |\   |  :      :    :  |\   |     |             :             |'
write(31,*)'|   |   _   :--   : \  :  |      |    |  : \  :     :     __    .-`             |'
write(31,*)'|   :   :   |     |  \ |  :      :    :  |  \ |     |          |                |' 
write(31,*)'|    `--`    --   :   \:   `--`   `--`   :   \:     :          :__.             |'
write(31,*)'|                                                                               |'
write(31,*)'|      (A program for calculation of Optimal Genetic Contributions)             |' 
write(31,*)'|                              Version 1.0                                      |'
write(31,*)'|                                                                               |'
write(31,*)'|                  Dept of Animal and Aguacultural Sciences                     |' 
write(31,*)'|                    Norwegian University of Life Sciences                      |'
write(31,*)'|                                    Norway                                     |'
write(31,*)'|                                                                               |'
write(31,*)'|                     Binyam Dagnachew and Theo Meuwissen                       |'
write(31,*)'|                         (Last update: October, 2015)                          |'           
write(31,*)' `-----------------------------------------------------------------------------`'

write(*,*)	
write(*,*) 'GENCONT-2 is running ...'
write(*,*) 'The results will be writen out to ',outfile
write(*,*)

	write(31,*)
    	write(31,*)
	write(31,*) '============================================================================='
   	write(31,*)
   	write(31,*) 'GENCONT starts @@ '
	write(31, 1000 )  today(2), today(1), today(3)
    	write(31, 1001) now
 1000 format ( '				Date: ', i2.2, '/', i2.2, '/', i4.4)
 1001 format ( '				Time: ', i2.2, ':', i2.2, ':', i2.2)
 	write(31,*)
 	write(31,*)'=============================================================================='
	
    	write(31,*)
	write(31,*)'Input parameters are read from "',Parafile,'"'
	write(31,*)
  

  	open(fh,file = Parafile)
	nprog=0; ngroup=0
  	minimize = 0
  	Mrelat = 0; deltaF=0; islow=0
	readin = 0
 	vanraden = 0
	henderson = 0
	inbreedingm = 0
	!cmax_male =1.0; cmin_male=0.0; cmax_female=1.0; cmin_female=0.0
	cmax_group =1.0; cmin_group= 0.0; s_group=0.0
	restrict = 0; a_half=0; n_young=0
	select_all_females = 0
	overlap = 0; ncand_ = 0; nped = 0; multi_stage = 0
	k_cand=0  !cmax_young =1.0; cmin_young = 0.0
	niter=100000 
	thold=1E-8
	n11_young=0;n12_young=0
	a11_young=0.0; a12_young=0.0
	
  do while (ios == 0)
    read(fh, '(A)', iostat=ios) str_in
    write(31,'(a)')trim(str_in)    
    if (ios == 0) then
        line_1 = line_1 + 1

    if (str_in(1:1) == '!'.OR. str_in(1:1) == '%') cycle
    if (str_in(1:1) == ' ') cycle
    
	CALL DELIM(str_in,str_out,NT,ICOUNT,IBEGIN,ITERM,ILEN,DLM)

     select case (trim(str_out(1)))
       case ('thold')
           read(str_out(2), *, iostat=ios) thold
       case ('ngroup')
           read(str_out(2), *, iostat=ios) nsex
       case ('nped')
           read(str_out(2), *, iostat=ios) nped
       case ('ncand')
           read(str_out(2), *, iostat=ios) ncand_
       case ('niter')
	   read(str_out(2), *, iostat=ios) niter
       case ('nprog')
           read(str_out(2), *, iostat=ios) nprog
	case('per_group')
		do i=1,nsex
           		if(ICOUNT .GE. i+1) read(str_out(i+1), *, iostat=ios) s_group(i)
		end do			
       case ('cmax')
		do i=1,nsex
           		if(ICOUNT .GE. i+1 .and.cmax_group(i) == 1.) read(str_out(i+1), *, iostat=ios) cmax_group(i)
		end do
		restrict(1)= 1
       case ('cmin')
		do i=1,nsex
           		if(ICOUNT .GE. i+1 .and.cmin_group(i) == 0.) read(str_out(i+1), *, iostat=ios) cmin_group(i)
		end do
		restrict(2)= 1
	case ('n_group')
		do i=1,nsex
			if(str_out(i+1) .EQ. 'OPT' .OR. str_out(i+1) .EQ. 'opt') then
				if(cmax_group(i) == 1.) cmax_group(i)= 1.
				if(cmin_group(i) == 0.) cmin_group(i)= 0.
			else
				read(str_out(i+1), *, iostat=ios) temp_val
				ngroup(i) = temp_val
				cmax_group(i)= 1.0/temp_val
				cmin_group(i)= (1.0/temp_val)
			end if
		end do
		restrict= 1
	case ('constraint')
            	if (str_out(2)=='relat') then
			read(str_out(3), *, iostat=ios) Mrelat
		end if
        	if (str_out(2)=='deltaf')read(str_out(3), *, iostat=ios) deltaF
            	if (str_out(2)=='minim') minimize=1
	case ('overlap', 'OVERLAP')
		if (str_out(2)=='yes') overlap=1
	case ('multistage', 'MULTISTAGE')
		!if (str_out(2)=='NO') multi_stage=0
		if (str_out(2)=='post_run') multi_stage=1
		if (str_out(2)=='run') multi_stage=2
    	case ('slow', 'SLOW')
		if (str_out(2)=='yes') islow=1
     	case ('Amatrix')
        	if (str_out(2)=='vanraden') vanraden =1
            	if (str_out(2)=='readin') then
			readin =1
			read(str_out(3), *, iostat=ios) a_file_name
			if(str_out(4)=='half')a_half=1				
		endif
            	if (str_out(2)=='henderson') henderson =1
		if (str_out(2)=='inbreeding') inbreedingm =1 
      	case ('ped_file')
           read(str_out(2), *, iostat=ios) ped_file
       	case ('data_file')
           read(str_out(2), *, iostat=ios) data_file
		ninp=icount-2
     		inp=0
     		do i=3,icount
!print*, str_out(i)(1:7)
       		if(str_out(i)(1:4)=='name' .or. str_out(i)(1:4)=='NAME')then
           			inp(i-2)=1
      	 		else if( str_out(i)(1:2)=='id' .or. str_out(i)(1:2)=='ID')then
           			inp(i-2)=2
       		else if( str_out(i)(1:3)=='ebv' .or. str_out(i)(1:3)=='EBV')then
           			inp(i-2)=3
       		else if( str_out(i)(1:5)=='avail' .or. str_out(i)(1:6)=='AVAIL')then
           			inp(i-2)=4
       		else if( str_out(i)(1:6)=='c_prev' .or. str_out(i)(1:6)=='C_PREV')then
           			inp(i-2)=5
       		else if( str_out(i)(1:4)=='c_max' .or. str_out(i)(1:4)=='C_MAX')then
           			inp(i-2)=6
				restrict(1)= 1
       		else if( str_out(i)(1:4)=='c_min' .or. str_out(i)(1:4)=='C_MIN')then
           			inp(i-2)=7
				restrict(2) = 1
       		else if( str_out(i)(1:3)=='age' .or. str_out(i)(1:3)=='AGE')then
           			inp(i-2)=8
       		else if( str_out(i)(1:3)=='sex' .or. str_out(i)(1:3)=='SEX')then
           			inp(i-2)=9
       		end if
    		 end do
      	case default
           print *, 'Skipping invalid INPUT at line', line_1
      end select
    end if
  end do

!print*, 'Reading file finshed'
write(31,*)'======================END READING PARAMETER FILE================================'

! Check for compulsory fields

	if(data_file == '') then
		write(32, *) 'data_file is missing'
		stop 'data_file is missing'
	end if
	if(ped_file == '') then
		write(32, *) 'a pedigree file is missing'
		stop 'ped_file is missing' 
	end if
	!if (vanraden == 1) then

	if(ncand_ == 0) then
		write(32, *) 'ncand is missing'
		stop 'ncand is missing'
	end if
	if(nped == 0) then
		write(32, *) 'nped is missing'
		stop 'nped is missing'
	end if 
	if(multi_stage > 0 .and. overlap .eq. 0) then
		write(32,*) 'Overlap is missing from the parameter file'
		stop 'Overlap is missing from the parameter file'
	end if
	if(multi_stage == 1 .and. all(ngroup ==0)) then
		write(32,*) 'The number of selected in completed selection stage should be provided'
		stop 'The number of selected in completed selection stage should be provided'
	end if
	if(overlap ==1 .and. count(inp==4)/=1) then 
		write(32, *)'Avail is missing from the data_file'
		stop 'Avail is missing from the data_file'
	end if
	if(overlap ==1 .and. count(inp==8)/=1) then
		write(32, *)'Age is missing from the data_file'
		stop 'Age is missing from the data_file'
	end if
	if(nsex .eq. 2 .and. count(inp==9)/=1) then
		write(32, *)'sex is missing from the data_file'
		stop 'sex is missing from the data_file'
	end if
	if(count(inp==3) /= 1) write(31,*) 'NOTE: Relationships will be minimized'

allocate (cand(ncand_), name(ncand_), sex(ncand_))
allocate (data_info(ncand_))

!print*, 'multistage = ', multi_stage
!print*,'nsex = ', nsex
!stop

!! ### Read input data  ### !!
	l=0

	open(11,file=data_file, status='old')
	data_info%age=1; data_info%sex=1; data_info%avail=1; data_info%c_min=0.0; data_info%c_max=0.0; 
	data_info%ebv=0.0; data_info%c_prev=0.0; data_info%avg_rel=0.0

	do 
 		read(11,'(a)',end=9) str_in
 		CALL DELIM(str_in,str_out,NT,ICOUNT,IBEGIN,ITERM,ILEN,DLM)

 		if(ICOUNT < ninp)then
			write(32,*) ' INPUTFILE does not contain all declared fields'
			stop ' INPUTFILE does not contain all declared fields'
		end if

 		l=l+1      
		if(l > ncand_)then
			write(32,*) 'The NCAND specified in the parameter file is too small'
			stop " ERR: NCAND too small in 'para.inp'"
		end if

 		do i=1,ninp
   			select case(inp(i))
     				case(1)
       					data_info(l)%name=str_out(i) 
     				case(2)
       					read(str_out(i),*)data_info(l)%id 
     				case(3)
       					if(minimize==0)read(str_out(i),*)data_info(l)%ebv 
     				case(4)
       					read(str_out(i),*)data_info(l)%avail
       					if(data_info(l)%avail==1)then
         					k_cand=k_cand+1; cand(k_cand)=l
       					end if 
     				case(5)
       					read(str_out(i),*)data_info(l)%c_prev 
     				case(6)
       					read(str_out(i),*)data_info(l)%c_max 
     				case(7)
       					read(str_out(i),*)data_info(l)%c_min 
     				case(8)
       					read(str_out(i),*)data_info(l)%age 
     				case(9)
       					read(str_out(i),*)data_info(l)%sex 
   			end select
 		end do       
	end do
	9 continue
	close(11)


! Determin number of available candidates
	nanim = l
	ncand = sum(data_info(1:l)%avail)
print*,'number of candidates =', ncand

	allocate (c(ncand), c_young(ncand), ebv(ncand), cand_id(ncand), aveA(ncand), anim(nanim))
	allocate (locate(ncand), Acand(ncand, ncand), A_anim(nanim, nanim), loc_post(ncand), cand_group(nsex))

! Select identification for candidates and animals
	anim = data_info(1:nanim)%id
	if (k_cand == 0) cand(1:ncand) = (/(i,i=1, ncand)/)
	cand_id = data_info(cand(1:ncand))%id
 
! Count number of candidates per group in the datafile

	do i=1,nsex
		cand_group(i)= count(data_info(cand(1:ncand))%sex == i)
	end do

	if(nsex .eq. 1) data_info%sex = 1

	allocate (Q(ncand,nsex), sum_cprev(nsex))
	allocate (ages(nsex), diff_ages(nsex), s(nsex), s_(nsex), w0(nsex))


! Creat age*sex groups

	if (overlap == 0) data_info(1:nanim)%age = 1
	
	j=0; age_ch=0; jj=0
	do i=1,nsex

		ages(i) = maxval(data_info(1:nanim)%age, mask=data_info(1:nanim)%sex == i)

		where (data_info(1:nanim)%sex == i+1)
			data_info(1:nanim)%age = data_info(1:nanim)%age + ages(i)
		end where

		do jj=jj+1,ages(i)
			age_ch=count(data_info(1:nanim)%age == jj)
!print*, ' age_ch',  jj, age_ch
	
			if(age_ch == 0) then
				where (data_info(1:nanim)%age .GT. jj)
					data_info(1:nanim)%age = data_info(1:nanim)%age-1
				end where
				ages(i)=ages(i)-1
			end if
		end do
		
		diff_ages(i)=ages(i)-j
		j=ages(i); jj=ages(i)
	end do

	nages = maxval(data_info(1:nanim)%age)
	if(multi_stage == 0 .or. nsex < 3) ages(2)=nages

!print*, 'diff_ages = ', diff_ages
!print*, 'ages = ', ages
!print*, 'nages = ', nages
!stop
        
	allocate (nbar(nages, nages), abar(nages, nages), cbar(nages), xbar(nages), char_mat(maxval(diff_ages), nsex))
 	allocate (CprevA12(nages), ACprev(ncand), A12(ncand, nages), n12(ncand, nages))
	allocate (w(nages), cbar_1(nages), avg_rel(0:ncand), A12c2(ncand), w_young(nages))
	allocate (cmaxx(ncand), cminn(ncand), cmin(ncand), cmax(ncand), sum_cprev_age(nages))


	!if (overlap == 1) niter=100000

!***======== Check for previous contributions for each age*sex group ========***

	sum_cprev_age = 0
	do i=1, nanim
		if(data_info(i)%c_prev .GT. 0) &
		sum_cprev_age(data_info(i)%age)= sum_cprev_age(data_info(i)%age) + data_info(i)%c_prev 
	end do

	do i=1,nsex
		sum_cprev(i) = sum(data_info%c_prev, mask= data_info%sex ==i) 
	end do

	do i=1,nsex
		if(sum_cprev (i) > 1.0) then
			write(32, *) 'ERROR: Sum of of previous contribution > 1 for group',i
			stop 'ERROR: Sum of Cprev > 1'
		end if
	end do

!*** ======== Define maximum and minimum contribution of individuals =========****

	!!!!!===== First collect the cmax ======!!!!!

do i=1,nsex	
	if(cmax_group(i) .gt. 0) then
		where((data_info(1:nanim)%c_max .GE. 1.0 .OR. data_info(1:nanim)%c_max .LE. 0.0) .AND. data_info(1:nanim)%sex ==i)
			data_info(1:nanim)%c_max = cmax_group(i)
		end where
	end if
end do 

	cmax(1:ncand)=data_info(cand(1:ncand))%c_max
    
	!!!!!===== Second collect the cmin ======!!!!!

do i=1,nsex
	if(cmin_group(i) .GT. 0.) then
		where((data_info(1:nanim)%c_min .GE. 0.) .AND. data_info(1:nanim)%sex ==i)
			data_info(1:nanim)%c_min = cmin_group(i)
		end where
	end if
end do

	cmin(1:ncand)=data_info(cand(1:ncand))%c_min

!print*, cmin_group(1:nsex)
!print*, data_info(1:nanim)%c_min !cmin
!print*,
!print*, cmin
!stop

!*** ======== Calculate (obtain) relationship matrix =========***

	if (vanraden == 1) then 
   		CALL vanraden_s (nped, nanim, A_anim, anim, ped_file)
   	elseif (readin == 1) then
		if (a_half == 1) then
			CALL readin_half (A_anim, nanim, a_file_name)
		else
    			CALL readin_s (A_anim, nanim, a_file_name)
		end if
	elseif ((nped .GT. 70000) .OR. (inbreedingm == 1)) then
		CALL inbreeding_s (nped, nanim, A_anim, anim, ped_file)
    	else 
		CALL A_mat (A_anim, nanim, anim, nped, ped_file)
	end if 
!print*, 'DONE'
!stop

!*** ======== Make list of previously contributed animals ======== ****

	iprev = count(data_info(:)%c_prev  > 0)
	if(iprev > 0)then
  		allocate(Ctemp(nanim,nages),Ntemp(nanim,nages),list(nanim))
 		 Ctemp=0; Ntemp=0
  		list(1:iprev)=pack((/(i,i=1,nanim)/),mask=data_info(:)%c_prev > 0)
	end if

	CprevACprev=0; ACprev=0; CprevA12=0
	a12=0; n12=0; abar=0; nbar=0
	data_info%avg_rel=0.0
	ii=0 !ii=cand_number

	do i=1,nanim
		if(data_info(i)%avail==1)ii=ii+1
	  	jj=0
	  	do j=1,i
	    		aa=A_anim(i,j)

	    		if(data_info(j)%avail==1)then
	      			jj=jj+1
	    		end if

	    		if(i/=j)then
	      			abar(data_info(i)%age,data_info(j)%age)=abar(data_info(i)%age,data_info(j)%age)+aa
	      			abar(data_info(j)%age,data_info(i)%age)=abar(data_info(j)%age,data_info(i)%age)+aa
	      			nbar(data_info(i)%age,data_info(j)%age)=nbar(data_info(i)%age,data_info(j)%age)+1.
	      			nbar(data_info(j)%age,data_info(i)%age)=nbar(data_info(j)%age,data_info(i)%age)+1.
	    		end if
	
	    		if(data_info(i)%avail==1)then
	      			a12(ii,data_info(j)%age)=a12(ii,data_info(j)%age)+aa
	      			n12(ii,data_info(j)%age)=n12(ii,data_info(j)%age)+1
	      			ACprev(ii)=ACprev(ii)+aa*data_info(j)%c_prev
	    		end if
	
	    		if(data_info(j)%avail==1 .and. i/=j)then
	      			a12(jj,data_info(i)%age)=a12(jj,data_info(i)%age)+aa
	      			n12(jj,data_info(i)%age)=n12(jj,data_info(i)%age)+1
	      			ACprev(jj)=ACprev(jj)+aa*data_info(i)%c_prev
	    		end if

	    		if(iprev > 0)then
	     			if(any(list(1:iprev)==i))then   
	      				Ctemp(i,data_info(j)%age)=Ctemp(i,data_info(j)%age)+aa
	      				Ntemp(i,data_info(j)%age)=Ntemp(i,data_info(j)%age)+1
	     			end if

	     			if(any(list(1:iprev)==j) .and. i/=j)then
	      				Ctemp(j,data_info(i)%age)=Ctemp(j,data_info(i)%age)+aa
	      				Ntemp(j,data_info(i)%age)=Ntemp(j,data_info(i)%age)+1
	     			end if
	    		end if

	    		fact=1.; if(i/=j)fact=2.
	    		CprevACprev=CprevACprev+data_info(j)%c_prev*aa*data_info(i)%c_prev*fact
	    		if(i/=j)data_info(i)%avg_rel=data_info(i)%avg_rel+aa
	    		if(i/=j)data_info(j)%avg_rel=data_info(j)%avg_rel+aa
	  	end do
	end do

	data_info%avg_rel = data_info%avg_rel/(nanim-1)
	where(nbar > 0)abar = abar/nbar
	where(n12 > 0)a12 = a12/n12

	if(iprev > 0)then
	  	where(Ntemp > 0) Ctemp = Ctemp/Ntemp
	  	CprevA12 = matmul(data_info(1:nanim)%c_prev,Ctemp)
	  	deallocate(Ctemp,Ntemp,list)
	end if


!***========= Get prior for contribution of ageclasses (cbar) ========= ****

	Acand=A_anim(cand(1:ncand), cand(1:ncand))
	deallocate(A_anim)

	cbar=0
	if(all(cbar <= 0))then
		
  		do i=1,nages
    			cbar(i) = count(data_info(cand(1:ncand))%age==i)
		end do
!print*,'cbar =',  cbar
		j=0
		do i=1,nsex			
  			cbar(j+1:ages(i)) = cbar(j+1:ages(i))/sum(cbar(j+1:ages(i)))/nsex
			j=ages(i)
		end do
	end if
!print*, ages
!print*,'cbar =',  cbar
!print*, 'sum(cbar)', sum(cbar(1:4))
!stop
!***=========== Collect constraints ============= ****

    	xbar=0   !contribution that abgeclass is still going to make

  	jj=0; ii=0	
	do j=1, nsex
		
		do i=jj+1,ages(j)
    			xbar(i)=sum(cbar(i:ages(j)))
  		end do
   		xbar(jj+1:ages(j))=xbar(jj+1:ages(j))/sum(xbar(jj+1:ages(j)))
		jj=ages(j)
	end do

!print*, xbar
!print*, 'sum(xbar) =', sum(xbar)
!stop 
  		
	PaveA=dot_product(cbar,matmul(abar,cbar)) 

	if(Mrelat==0)then

  		k = PaveA + deltaF*(2.- PaveA)

    		if (minimize == 1) then
      			k = PaveA
		end if
    	else
      		k = Mrelat
	end if

!print*,'current PaveA =', PaveA
!print*, 'the constraint is =', k
!stop

!***============ Construct design matrix for sex of the candidates ===========***
	Q=0

	if (nsex>1) then 

		do i=1,nsex
			where(data_info(cand(1:ncand))%sex == i)
				Q(:,i)=1		
			end where
		end do
	else
		Q=1
	end if

!do i=100,105
!	print*, Q(i,1), Q(i,2)
!end do
!stop
!print*, size(Q,1), size(Q,2)
!stop

!***======== determing contribution per group =======****

	if (any(s_group .GT. 0.)) then
		s=s_group(1:nsex)
	else
		s=1./nsex
	end if

!print*, 's', s
!stop

!***============ Get weight for each age*sex group ==============****

! ============ At first get Weights from Cbar ================

	niter_w=1
	do i_w=1,niter_w
!print*, 'the weight round is =', i_w

		w0=1./nsex
		sumw=0
!print*, 'before w0 = ', w0
!print*, 'nages =', nages
!print*, 'ages =', ages
!print*, 'multi_satge 1', multi_stage

	j=0
	do i=1,nsex
  		i_sex=i
		do i_age=j+1,ages(i)
  			w0(i_sex)=w0(i_sex)-cbar(i_age)
  			w(i_age)=w0(i_sex) 
			sumw=sumw+w(i_age)
			j=j+1

!print*, 'i_sex =', i_sex
!print*, 'i_age =', i_age
!print*, 'j =', j
		end do
	end do

		sumw=sumw +1.  !add contributions of age class 0
		w=w/sumw

	do i=1,nsex
		s_(i)=1./nsex/sumw  !sum contributions of age class 0 males
	end do

!print*,'w = ', w
!print*,'sumw = ', sumw
!stop
	! Calculate c2a22c2
  		c2a22c2=dot_product(W,matmul(Abar,W))  &   
       			+ CprevACprev*s_(1)**2 + 2.*dot_product(CprevA12,W)*s_(1)
  		A12c2=matmul(a12,W) + ACprev*s_(1)
!print*,
!print*, 's_ before = ', s_
!print*,

	! Correct s_ for previousc onctributions of each sex group
		s_=s_*(1.-sum_cprev)

	! Correct cmax and cmin for s  
  		cmaxx=cmax*s_(1)
	  	cminn=cmin*s_(1)
!print*,'cmax', cmaxx(1:5)
!stop

!print*, 's_ after= ', s_
!print*, 'sum prev =', sum_cprev
!print*, 'c2a22c2 = ', c2a22c2
!stop

	! Call for the main subroutine whcih deals wiht the optimization problem
 
        CALL overl_iter (Acand, Q, k, s_, ncand, nsex, data_info(cand(1:ncand))%ebv,data_info(cand(1:ncand))%sex, c, niter, &
		thold, Error_c, round,cmaxx, cminn, c2a22c2, a12c2, islow, restrict, ngroup, cand_group)

	! Get new cbar (cbar_1)
	
	j=0
	do i=1,nsex
  		cbar_1(j+1:ages(i))=sum_cprev_age(j+1:ages(i))*s(i)
		j=ages(i)
	end do
	
  	do i=1,ncand
		do j=1,nsex
			if(Q(i,2) > .0)i_s=j
		end do
				
		i_age = data_info(cand(i))%age
		cc =(1.-sum_cprev(i_s))*s(i_s)*c(i)/s_(i_s)
		cbar_1(i_age)=cbar_1(i_age) + cc
	end do

	j=0
	do i=1,nsex
		cbar_1(j+1:ages(i))=cbar_1(j+1:ages(i))/sum(cbar_1(j+1:ages(i)))/nsex
		j=ages(i)
	end do

!if (i_w >1 ) then
!print*, 'the cbars ='
!do i=1,nages
!print*,cbar(i)
!end do
!stop
!end if

! Check convergence
  		o_sse=dot_product(cbar-cbar_1,cbar-cbar_1)
  		o_sst=dot_product(cbar_1,cbar_1)
	  	o_sse=o_sse/o_sst
	  	cbar=cbar_1
  		if(o_sse<1.e-3)exit
	end do  !for weight iteration 
!print*
!print*,'DONE GENCONT!'

!print*, 'w0',w0

!do i=1,nages
!print*,w(i)
!end do

!stop

!print*, 'the cbars ='
!do i=1,nages
!print*,cbar(i)
!end do

!***======== Get average relationship of solution ========***

	avg_rel(1:)=matmul(Acand,c)
	avg_rel(0)=dot_product(c,avg_rel(1:)) + 2.*dot_product(c,a12c2) + c2a22c2
	avg_rel(1:)=avg_rel(1:) + a12c2
!print*, avg_rel(0)
!print*, 'c2a22c2 = ', c2a22c2
!print*, 'a12c2 = ', a12c2


! *** ===== Calculate increased average relationship due to selection at stage 1 =====***

IF (multi_stage == 1) THEN
	open(34, file ='Pre_run.out')

	young_select = count(data_info(cand(1:ncand))%sex==2 .AND. c .gt. 0.)
	allocate(A_young(young_select,young_select))
!print*, 'n_selected',young_select

	c_young=0
	do i=1,ncand
		if(data_info(cand(i))%sex == 2) then
			c_young(i)=c(i)
		end if
	end do

avg_young = dot_product(c_young,matmul(Acand,c_young)) + 2.*dot_product(c_young,a12c2)

!print*, 'avg_young',avg_young
!print*, 'avg_young/avg_rel(0)', (avg_young/avg_rel(0))*100



!w_young=0; w_young(nages) =w(nages)
!print*,'w_young', dot_product(w_young,matmul(Abar,w_young)) 

!	avg_young=dot_product(c_young,matmul(Acand,c_young)) + (2.*dot_product(c_young,a12c2)) &
!			+ dot_product(w_young,matmul(Abar,w_young)) 
!print*, 'avg_young', avg_young
!print*, 'avg_young/avg_rel(0)', (avg_young/avg_rel(0))*100

	write(34,*) 'number of selected at stage 1 = ', young_select
	write(34,*)'Relationship incease due to stage 1'
	write(34, 34) avg_young
	34 format(f15.11)
	close(34)

END IF ! for the multi_stage selection pre_run

!print*,'before sum', sum(c)
!print*, 'sum_cprev',sum_cprev 
!print*, 's = ',s
!print*, 's_ =',s_
!stop

! ========== Make sure the contributions sum to s ============================***
	
  	do i=1,ncand
		do j=1,nsex
			if(Q(i,2) > .0)i_s=j
		end do
				
		c(i)=(1.-sum_cprev(i_s))*s(i_s)*c(i)/s_(i_s)
	end do

covA=sum(c,1)
!print*, 'before sum(c) = ', covA
!	do i=1, ncand
!		c(i)=c(i)/covA
!	end do

!print*,'after sum', sum(c)

!***======== Check the status of convergence and give report ========***

        if (Error_c == 0) then
         	write(31,*)
		write(31,*) 'Converged after', round, ' itrations'
            	write(31,*)
        else          
         	write(31,*)
		write(31,*) ' WARNING!!! Do not converge with the current criteria'
         	write(31,*) '	Reduce "thold" or increase "niter"'
            	write(31,*)
         	write(31,*) ' The solutions after',niter, ' itrations are provided here' 
         	write(31,*)
        end if

		write(31,'(A44,f13.10)') ' Population Average Relationship (current) = ', PaveA
   		write(31,'(A44,f13.10)') ' Constraint = ', k


!*** ======== Summerize results and Write solutions to an output file ======== ****!

! Calculate the average relationship among selected parents
	PaveA_selct = dot_product(c , matmul(Acand,c))
        
write(31,*)

do i=1,nsex
	write(31,'(A31,i2,A3,i7)') 'Number of candidates in group',i,' = ',cand_group(i)
end do
write(31,*)


!***======== writing outputs for the age-groups ========== ***

write(31,*)
write(31,*) 'Average r/ship per group'
write(31,*)

j=1; line1=''
	write(line1(j:j+13),'(a13)') "Age group"
		j=j+13

	do i=1, nsex
		write(line1(j:j+12),'(A10, i2)') "Group",i
			j=j+12
	end do
write(31,'(a)') trim(line1)
!print*, diff_ages
!stop	

j=1; char_mat=-1.1;
	do i=1, nsex
		do kk=1,diff_ages(i)
			char_mat(kk, i) = abar(j,j)
			j=j+1
		end do	
		
	end do

	do i=1,maxval(diff_ages)
		j=1; line1=''	
		write(line1(j:j+13),'(i13)') i
		j=j+13
		
		do kk=1,nsex
			if(char_mat(i,kk) .eq. -1.1) then
				write(line1(j:j+12),'(A12)') '------'
				j=j+12
			else
				write(line1(j:j+12),'(f12.7)') char_mat(i,kk)
				j=j+12
			end if
		end do

	write(31,'(a)') trim(line1)
	end do
		
!***========= Writing out solutions ================================================

! Calculate genetic merit of the selected candidates
    	gen_gain=dot_product(c,data_info(cand(1:ncand))%ebv) &
        			+ dot_product(data_info(1:nanim)%c_prev/nsex,data_info(1:nanim)%ebv)

write(31,*)
write(31,*)
write(31,'(A44,f13.10)') 'Population Average Relationship (solution) = ', avg_rel(0)
write(31,*)
write(31,'(A44,f13.6)') 'Genetic merit of the parents = ', gen_gain
write(31,*)
write(31,*)
 
!## Determine the number of outputs

	nout=0
	if(any(inp==1)) then
		iout(1)=1		! For name
		nout=nout+1
	else 
		iout(1)=-1		! For ID
		nout=nout+1
	end if
		iout(2)=2		! For contribution (c)
		nout=nout+1
	if(any(data_info(1:nanim)%c_prev > 0)) then
		iout(3)=3		! For previous contribution (Cprev)
		nout=nout+1
	end if
	if(any(data_info(1:nanim)%c_min > 0)) then
		nout=nout+1		! For c minimum constraint
		iout(nout) = 4
	end if
	if(any(data_info(1:nanim)%c_max > 0 .AND. data_info(1:nanim)%c_max < 1.)) then
		nout=nout+1		! For c max constraint
		iout(nout) = 5
	end if
	if(any(inp==3)) then
		nout=nout+1		! For EBV
		iout(nout) = 6
	end if
	nout=nout+1; iout(nout) = 7 	! For average r/ship
	if(any(inp==8)) then
		nout=nout+1
		iout(nout)=8		! For age
	end if
	if(any(inp==4)) then
		nout=nout+1
		iout(nout)=9		! For avail
	end if

j=1; line1=''

do i=1, nout
	if(abs(iout(i))==1) then
		write(line1(j:j+20), '(a)') ' Name  '
			j=j+20
	end if
	if(abs(iout(i))==2) then
		write(line1(j:j+10), '(a)') ' %_progeny'
			j=j+10
	end if
	if(abs(iout(i))==3) then
		write(line1(j:j+10), '(a)') ' %_Cprev'
			j=j+10
	end if
	if(abs(iout(i))==4) then
		write(line1(j:j+10), '(a)') ' %_Cmin'
			j=j+10
	end if
	if(abs(iout(i))==5) then
		write(line1(j:j+10), '(a)') ' %_Cmax'
			j=j+10
	end if
	if(abs(iout(i))==6) then
		write(line1(j:j+10), '(a10)') '  EBV '
			j=j+10
	end if
	if(abs(iout(i))==7) then
		write(line1(j:j+10), '(a)') ' Ave_relat'
			j=j+10
	end if
	if(abs(iout(i))==8) then
		write(line1(j:j+10), '(a)') '  Age '
			j=j+10
	end if
	if(abs(iout(i))==9) then
		write(line1(j:j+10), '(a)') '  Avail '
			j=j+10
	end if
end do

do kk=1,nsex

write(31,*)

write(31,'(A31,i2,A3,i5)') 'Number of animals in Group',kk, '=',count(data_info(1:nanim)%sex==kk)
write(31,'(A31,i2,A3,i5)') 'Number of candidates in Group',kk, '=',count(data_info(cand(1:ncand))%sex==kk)
write(31,'(A31,i2,A3,i5)') 'Number of selected from Group',kk, '=',count(data_info(cand(1:ncand))%sex==kk .AND. c > 0)

write(31,*)
write(31,*) '____________________________________________________________________'
write(31,'(a)') trim(line1)
write(31,*) '____________________________________________________________________'
write(31,*)
	jj=0
	do ii=1, nanim
		if(data_info(ii)%avail == 1) jj=jj+1
		if(data_info(ii)%avail == 1 .OR. data_info(ii)%c_prev > 0) then
			if(data_info(ii)%sex==kk) then
				line=''; j=1
				do i=1, nout
					if(iout(i)==1) then
						write(line(j:j+20), '(a)') data_info(ii)%name
						j=j+20
					else if(iout(i)==-1) then
						write(line(j:j+20), '(i8)') data_info(ii)%id
						j=j+20
					end if
					if(iout(i)==2) then
						if(data_info(ii)%avail == 1) then
							write(line(j:j+10), '(f10.3)')c(jj)*100.*nsex
							j=j+10
						else
							write(line(j:j+10), '(f10.3)') 0.0
							j=j+10
						end if
					end if
					if(iout(i)==3) then
						write(line(j:j+10), '(f10.3)') data_info(ii)%c_prev*100.
						j=j+10
					end if
					if(iout(i)==4) then
						write(line(j:j+10), '(f10.3)') data_info(ii)%c_min*100.
						j=j+10
					end if
					if(iout(i)==5) then
						write(line(j:j+10), '(f10.3)') data_info(ii)%c_max*100.
						j=j+10
					end if
					if(iout(i)==6) then
						write(line(j:j+10), '(f10.3)') data_info(ii)%ebv
						j=j+10
					end if
	
					if(iout(i)==7) then
						if(data_info(ii)%avail == 1) then
							write(line(j:j+10), '(f10.3)')avg_rel(jj)
							j=j+10
						else
							write(line(j:j+10), '(f10.3)') 0.0
							j=j+10
						end if
					end if
		
					if(iout(i)==8) then
						write(line(j:j+10), '(i5)') data_info(ii)%age
						j=j+10
					end if
					if(iout(i)==9) then
						write(line(j:j+10), '(i5)') data_info(ii)%avail
						j=j+10
					end if
				end do
					write(31, '(a)') trim(line)
			end if
		end if 
	end do
write(31,*) '____________________________________________________________________'
write(31,*)
end do

!print*,'sum of c', sum(c(1:ncand)) 

print*

do i=1,nsex

	print*, 'Number of selected in group',i,' = ',count(data_info(cand(1:ncand))%sex ==i .AND. c > 0)
!print*
end do

print*
print*, 'The genetic merit of the parent = ', gen_gain
print*
print*, '***', ' Solutions are written in the output file ==> ', outfile

	today=0
	now=0
      		call idate(today)
      		call itime(now)
      
write(31,*) '==========================================================================='
write(31,*)
write(31,*) 'GENCONT2 ends @@ '
write(31, 1002)  today(2), today(1), today(3)
write(31, 1003) now
 	1002 format ( '				Date: ', i2.2, '/', i2.2, '/', i4.4)
 	1003 format ( '				Time: ', i2.2, ':', i2.2, ':', i2.2)
write(31,*)
write(31,*)'============================================================================'
    
close (31)
close (32)

!write output for plotting
	!open(33, file='con_ebv.out')
	!	do row=1,ncand
	!		write(33,23) data_info(row)%id, c(row), data_info(row)%ebv , c(row)*100.*nsex , data_info(row)%sex 
	!	end do
	!	23 format(i8,f12.7, 2f12.4, i3)
!
	!close(33)

contains

!=============================================================================================
!
!	This subroutine calculate optimum contricbution of animals 				===
! 		constraining on average relation among candidates				===
!
!=============================================================================================

Subroutine overl_iter (A, Q, k, ones, ncand, nsex, ebv, sex, c, niter, thold, Error_c, round, cmax, cmin, & 
				c2a22c2, a12c2, islow, restrict, ngroup, cand_group)

	Implicit none

! Define input parameters for the subroutines

	integer, intent (in)					:: ncand, islow, niter, restrict(2), ngroup(20), nsex
	integer, intent (inout)					:: round, Error_c
	integer, intent(in)					:: cand_group(nsex), sex(ncand)
   	real(8), dimension(ncand), intent(out)			:: c
	real(8), intent(in)					:: k, thold, c2a22c2
	real(8), dimension(ncand), intent(in)			:: ebv, a12c2
    	real(8), dimension(nsex), intent(in)			:: ones
	real(8), dimension (ncand, nsex),intent(in)		:: Q
	real(8), dimension (ncand, ncand),intent(in) 		:: A
	real(8), dimension (ncand), intent(inout), optional	:: cmax, cmin

! Define related parameters 

	integer								:: i, j, x, xx,xxx, n_run, nin, nout, y 
	integer								:: select_all_group(20), ii, jj, ibig, iter_c, nsex_1
	real(8)								:: k1,threshold, big, zero, c2a22 
	integer, dimension (ncand)	 				:: locate
	integer, dimension(ncand)					:: loc_post, rem_post
   	integer, dimension (:), allocatable				:: in, out
   	real(8), dimension (:), allocatable 				:: a_c, c_set, ones_r , c_store

! Allocate respective dimensions for the array

    	allocate (ones_r(nsex))
	allocate (a_c(ncand), c_set(ncand), c_store(ncand))
	allocate (in(ncand), out(ncand))

	n_run = niter 
	n_run = 500 
	!if(overlap==1) n_run = niter 
	x=ncand; y=0; zero=0.
	locate(1:x)=(/(i,i=1,x)/)  !all animals are in solution
	ones_r=ones; nsex_1=nsex
	a_c=0.; k1=0.; c2a22=0.
	select_all_group=0

! Check if one or more groups of candidates are all selected 

!print*, 'cand_group', cand_group
!print*, 'ngroup', ngroup(1:nsex)

ii=0; jj=0; c=0.

	do i=1,nsex

		if (ngroup(i) .EQ. cand_group(i)) then
			select_all_group(i) = 1
		end if
	end do 

IF (any(select_all_group .GT. 0) ) THEN
	do i=1,nsex
		if (ngroup(i) .GT. cand_group(i)) then
			write(32, *)'# of required candidates is more than availabel candidates'
			stop '# of required candidates is more than availabel candidates'
		end if		
		
		if(select_all_group(i) == 1) then 

			do j=1, ncand
				if (sex(j)==i) then
					jj=jj+1
					loc_post(jj)=locate(j)
					!c(j)=cmax(j)
					c(j)=ones_r(i)/cand_group(i)
				end if
			end do
		else
			do j=1, ncand
				if (sex(j)==i) then
					ii=ii+1
					rem_post(ii)=locate(j)
				end if
			end do

		end if
	x=ii
	Y=jj
	end do

	locate(1:x)=rem_post(1:x)

!print*, loc_post(1:y)
!print*, c
!print*, rem_post(1:x)
!print*, 'x is =', x
!print*, 'y is =', y

!stop	

	nsex_1=nsex-sum(select_all_group)
	a_c(rem_post(1:x))=matmul(A(rem_post(1:x), loc_post(1:y)), c(loc_post(1:y)))
	c2a22= dot_product(c(loc_post(1:y)), matmul(A(loc_post(1:y), loc_post(1:y)), c(loc_post(1:y)))) &
			+ 2*dot_product(c(loc_post(1:y)), a12c2(loc_post(1:y)))


!print*, '1 a_c is = ', a_c
!print*, 'c2a22c2 =', c2a22c2
!print*, locate(1:x)

END IF

!print*, 'nsex_1 =', nsex_1
!stop

! Impose cmax & Cmin restriction (Cmax first)


	where(cmax<=0)cmax=1.
	if(restrict(2) /= 0)then
  		if(restrict(1) /= 0)where(cmax<abs(cmin))cmax=abs(cmin)
	end if

! Initiate some variables 
	xx = x
	k1 = k-c2a22c2 -c2a22
	nout = 0; nin=xx; icount=0
	a_c(locate(1:xx)) = a_c(locate(1:xx)) + a12c2(locate(1:xx))
	c_set = 0
	c_store = c
!print*, 'the const k1 = ', k1
!print*, '2 a_c is = ', a_c
!stop

!#####################################################
!***** The main loop begins here **********	######	
!#####################################################

do iter_c=1,100000
	xxx=xx

	!print*, 'the xx is ', xxx
	!print*, 'itreation is ', iter_c

	if(iter_c .GT. 1) n_run=500

  	CALL oc_iter (locate, c,c_set, k1, ebv,  a_c, Q, A, niter, n_run, ones_r, ncand, nsex_1, &
			xxx, Error_c, thold, round, restrict, islow)

	
!*** If there is restriction on cmax and cmin, Exit ***

	if(restrict(1) == 0 .and. restrict(2) == 0) then
	 	!c =c_set
		c_store= c_set
		nin=0
		exit
	endif

	nin=xx
!if(iter_c == 5 .or. iter_c == 6) then
!	print*, c(1:xxx)
!	print*, locate(1:xxx)
!end if

!print*, 'Here 1'

!*** If there is restriction, deal with cmax first ***

	IF(restrict(1) == 1 .or. restrict(2) /= 0) THEN
  		nin=0
  		threshold=0.
  		if(islow==1)threshold=max(maxval(c_set(locate(1:xx)) - cmax(locate(1:xx))),zero)
  	
		do i=1,xx
    			ii=locate(i)
   			if(c_set(ii)-cmax(ii)>=threshold)then    !fix solution to cmax
				!print*, locate(ii)
				!print*, c_set(ii)
				!print*, c_set(ii)-cmax(ii)
				!STOP
      				c_store(ii)=cmax(ii)
		! Adjust KK:
      				k1=k1-A(ii,ii)*cmax(ii)**2 -2*cmax(ii)*a_c(ii)
		! Adjust a_c:
      				a_c(locate(1:xx))=a_c(locate(1:xx)) + cmax(ii)*A(locate(1:xx),ii)
		! Adjust ones:
      				ones_r= ones_r - Q(ii,:)*cmax(ii)
    			else
      				nin=nin+1
      				in(nin)=i
    			end if
  		end do
!print*, 'Here 2'
print*, 'the xx is ', xx
print*, 'the nin is ', nin

		if(nin==xx .and. restrict(2)/=0)then

    			if(any(cmin < 0))then  !fix their solutions to abs(cmin(i))
      			nin=0
print*, 'I am Here'
      			do i=1,xx
        			ii=locate(i)
        			if(cmin(ii) < 0 .and. c(i) < -cmin(ii))then    !fix solution to -cmin
          				c_store(ii)=-cmin(ii)
		! Adjust K1:
          				k1=k1-A(ii,ii)*c_store(ii)**2 -2*c_store(ii)*a_c(ii)
		! Adjust a_c:
          				a_c(locate(1:xx))=a_c(locate(1:xx)) + c_store(ii)*A(locate(1:xx),ii)
		! Adjust ones_r :
          				ones_r = ones_r - Q(ii,:)*c_store(ii)
         	 			cmax(ii)=max(cmax(ii) + cmin(ii),zero)  !adjust cmax
          
          				if(cmax(ii) > 0.0)then                 !animal remains in IND list
            					nin=nin+1
            					in(nin)=i
          				end if
        			else
            				nin=nin+1
            				in(nin)=i
        			end if                  
      			end do
		if(nin ==1) exit
    			end if
  		end if

		!if(nin==xx .and. ngroup_in(
	END IF
!print*, 'the nin is ', nin
!print*
!print*, in(1:nin)
!print*
!print*, 'the const is ', k1
!print*, 'the ones is ', ones_r

!*** If cmax is done, then deal with cmin ***

	IF(nin==xx .and. restrict(2) /= 0)THEN
  
  		if(nin==xx)then    !set solution with biggest CMINi-Ci to zero (1_by_1)
    			nin=0
    			big=1.e-6; ibig=0                   !1.e-6 leaves some margin

    			do i=1,xx        !exclude all zero solutions
    				if(c(i) > 0)then
       				nin=nin+1
        				in(nin)=i
        				if(big < cmin(locate(i))- c(i) .and. cmin(locate(i)) > 0)then
          					big = cmin(locate(i))- c(i)
          					ibig=i   
        				end if       
      				end if
   			end do 
print*, 'I am Here'
    			if(ibig==0)exit    
    			c_store(locate(ibig))=0             !set smallest solution to zero & exclude animal
    			nin=nin-1
    			in(ibig:nin)=in(ibig+1:nin+1)
  		end if
	ELSE
  		if(nin==xx)exit
	END IF

!if(iter_c .ge. 7) then
print*, nin
!stop
!print*, c
!print*
!print*, c_store
!endif

	if(nin==0)exit
	locate(1:nin)=locate(in(1:nin))
	xx=nin
if(any(ones_r <= 0)) exit
end do 

	if(nin > 0) c_store(locate(1:xx))=c_store(locate(1:xx))+ c_set(locate(1:xx))
	c = c_store
 	!open(53, file='test_c_1.out')
	!	do row=1,ncand
	!		write(53,53) c(row)
	!	end do
	!	53 format(f12.7)
!
	!close(53)
!print*, 'the c*ac is', dot_product(c,matmul(A,c))
!stop

end subroutine overl_iter

end program optim_cont
