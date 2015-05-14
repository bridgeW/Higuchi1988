MODULE higuchi1988_subs
CONTAINS
    !**********************************************
    !****SUBROUTINE check_in***********************
    !**********************************************
	SUBROUTINE check_in(one_day_data,x,status)
	!  Purpose:
	!    Check the input file when Opening,Reading &Ending
	!  References:
	!   Stephen J. Chapman,Fortran 95/2003 for Scientists and Engineers(Third Edition)    
	!	    McGraw-Hill Primis

    USE  higuchi1988_type
    
	IMPLICIT NONE
	CHARACTER(len=*),INTENT(IN)::one_day_data
	INTEGER(I1B),INTENT(OUT)::status        !0--Open success
	                                        !>0 Open or READ error
	                                        !=-1 READ Ending of FIle
    REAL(DP),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::x
    
    INTEGER(I4B)::i     !Loop Index	                                       
	INTEGER(I4B)::nvals	!Count the number of input-values
	REAL(DP)::hx,hy,z,h,f                !D,H,Z component of magnetic field
	
	REAL(DP)::temp=0.				!Temperary varible
	
    nvals=0
	!//Open the Input File
	OPEN(LU,FILE=one_day_data,STATUS='OLD',IOSTAT=status,ACTION='READ')

	!//Check to see the OPEN failed
	!//Check open error on opening Input file
	openIF:IF(status==0)THEN
	  !//Opening success
	  countLoop:DO
	    READ(LU,*,IOSTAT=status)
	  	IF(status/=0)EXIT	!//EXIT if not Valid
		nvals=nvals+1
	  ENDDO countLoop
	  
	  !//The readloop has terminated. Was it because of a READ error
	  !// or because of the end of file?
	  readIF:IF(status>0)THEN 
	    WRITE(*,2)nvals+1
		2 FORMAT(1X,'Reading Error Occured at INPUT Reading Line',I6/)
		WRITE(*,*)'INPUT Reading Status=',status
		WRITE(*,*)"Maybe some Non-real or Non-integer data In Input File."
		RETURN      !Return to main Program
	  ELSE readIF
	    WRITE(*,*)'Read file Status=',status
	    WRITE(*,3)nvals
	    3 FORMAT(1X,'End of INPUT File Reached. There Are ',I6,&
					' Values in the INPUT File.'/)
	  ENDIF readIF
	ELSE openIF
	  ! Opening failed
	  WRITE(*,4)status
	  4 FORMAT(1X,'Error Opening File: IOSTAT=',I6/)
	  WRITE(*,*)'File:',one_day_data,'Does not exist!'
	  RETURN    !Return to main Program
	ENDIF openIF
	
	ALLOCATE(x(nvals),STAT=status)
	REWIND(LU)
	readLoop:DO i=1,nvals
        
        !**********************************************************
        !> the input data only one column
        !**********************************************************
        READ(LU,*,IOSTAT=status)temp
        
        
	    !**********************************************************
	    !> NOTE:if the input file have more than one column******
	    !> for an example: the IAGA data which has D,H,Z and time
        !**********************************************************
	    ! READ(LU,*,IOSTAT=status)hx,hy,z,temp,f 
	    !**********************************************************
        
	  IF(status/=0)EXIT	!//EXIT if not Valid
	  x(i)=temp
	ENDDO readLoop
	
	CLOSE(LU)
		
	END SUBROUTINE check_in
	
	!**********************************************
	!****SUBROUTINE higuchi1988********************
    !**********************************************
	SUBROUTINE higuchi1988_v5(one_day_data,x,output_lmk,slope,rcoef,status)
	!  Purpose:
	!    Calculating the fractal dimension
	!  References:
	!	1.Higuchi,T.,1988.Approach to an irregular time 
	!	   series on the basis of fractal theory. 
	!	   Physica D 31.277-283
	!	2.Burlaga,L.F., Klein,L.W.,1986.Fractal structure of
	!	   the interplanetary magnetic field. J.Geophys.Res.91,374-350
	!	3.Gotoh,K.,Hayakawa,M.,Smirnova,N.A.,etal.2004. Fractal 
	!	   analysis of seimogenic ULF emissions. Physis and Chemistry of 
	!	   the Earth.29,419-424

    USE  higuchi1988_type
    
	IMPLICIT NONE
	CHARACTER(len=*),INTENT(IN)::one_day_data
	CHARACTER(len=*),INTENT(OUT)::output_lmk
	REAL(DP),INTENT(OUT)::slope
	REAL(DP),INTENT(OUT)::rcoef     !Relative coefficiency of Least-Square
	INTEGER(I1B),INTENT(OUT)::status            !status=0 for success
	REAL(DP),ALLOCATABLE,DIMENSION(:),INTENT(IN)::x
		                                !When INTENT is IN,x is not allowed to modified,
	                                !When INTENT is INOUT,x is allowed to modified.
	INTEGER(I4B)::nvals             !Count numbers of data
	INTEGER(I4B)::ndata
	
	REAL(DP)::Lmk=0.				!Length of curves at (m,k)
			
	REAL(DP),ALLOCATABLE,DIMENSION(:)::log2_Lmk    !Log2(Lmk)
	                                !!Dimension varies with the the loop of m
	REAL(DP)::std_log2Lmk           !Standard deviation of Log2(Lmk)
	INTEGER(I4B)::k,m,i		        !Control-parameters of loops
	INTEGER(I4B),ALLOCATABLE,DIMENSION(:)::interval		!Intervals for averaging
	                                         !Dimension varies with the numbers of Inputdata
	
	REAL(DP)::temp=0.				!Temperary varible
	
	INTEGER(I2B),PARAMETER::nstarts=11		!Inorder to make INT(2**(nstarts-1)/4.)=5
										    ! So nstarts=11
	INTEGER(I2B)::nends			
	
	!//FOR Least-square
	REAL(DP)::intercept=0.
	REAL(DP)::sum_k				!Sum of Log2(Interval)
	REAL(DP)::sum_L              !Sum of Log2(Average of Lmk)
	REAL(DP)::sum_kL             !Sum of Log2(Interval)*Log2(Average of Lmk)
	REAL(DP)::sum_k2             !Sum of (Log2(Interval))**2
	REAL(DP)::sum_L2             !Sum of (Log2(Average of Lmk))**2
	REAL(DP)::ave_k              !Average of Log2(Interval)
	REAL(DP)::ave_L              !Average of Log2(Average of Lmk)
		
	!Initialize output file's name
	output_lmk='LMK'//TRIM(one_day_data)
    
    !Get the size of x(INput dummy argument)
    nvals=SIZE(x)
    WRITE(*,*)'Sub higuchi1988: Size(x)=',nvals
	
	!//ALLOCATE memory to interval & Initialize it
	nends=4*INT(LOG10(nvals/10.)/LOG2)+1      !Inorder to make INT(2**(nends-1)/4)<=(nvals/10)
										    ! So nends=29
    ALLOCATE(interval(nends-nstarts+5),STAT=status)
	interval=(/1,2,3,4,(INT(2**((i-1)/4.)),i=nstarts,nends)/)
	WRITE(*,*)"SUB_higuchi1988_v*: ALLOCATE interval's Size=",SIZE(interval)
	
	
	!//Open output file to store Lmk.
	OPEN(LU,FILE=output_lmk,STATUS='REPLACE',IOSTAT=status)
	WRITE(*,*)'SUB: higuchi1988',output_lmk,'Output File Status=',status
	WRITE(LU,5)
	5 FORMAT(1X,'%',T5,'K',T15,'CURVE LENGTH',10X,'Standard Deviation'/&
			 1X,'%','======',6X,'===========',10X,'===================')

	
	!//Open a file to save each step of m
	OPEN(LU+1,FILE='Each_K_to_Lm.txt',STATUS='REPLACE',ACTION='READWRITE',IOSTAT=status)
	WRITE(*,*)"Open FILE Each_K_to_Lm.txt's status=",status
	WRITE(LU+1,6)
	6 FORMAT(1X,'%Interval Time',1X,'Initial Time',3X,'log2LMK'/,&
			1X,'%======',6X,'======',12X,'======')

	!//Calculate Lmk
	WRITE(*,*)"Calculating the Curve Length......"

    !Initialize Least-square's Parameters ot sum.
    sum_k=0.
	sum_L=0.
	sum_kL=0.             
	sum_k2=0.             
	sum_L2=0.            
	ave_k=0.              
	ave_L=0. 
    
    WRITE(*,*)'SIZE(interval)=',SIZE(interval),interval
	LOOP_k: DO k=1,SIZE(interval)
	  
	  !ALLOCATE memory to log2_Lmk
	  ALLOCATE(log2_Lmk(interval(k)),STAT=status)
	  allocated_okIF:IF(ALLOCATED(log2_Lmk).AND.ALLOCATED(x))THEN
	    ndata=0       !Fix on where the data in log2_Lmk
	  
	    !Initialize Average's Parameters ot sum.
        Lmk=0.
         
	    LOOP_m: DO m=1,interval(k)
	      !Initialize each-step-Lmk's Parameters ot sum.
	      temp=0.
	    
	  	  LOOP_i: DO i=1,INT((nvals-m)/interval(k))
	  	  
	  	    temp=temp+ABS(x(m+i*interval(k))-x(m+(i-1)*interval(k)))*	&
					(nvals-1)/(INT((nvals-m)/interval(k))*interval(k))/interval(k)
	
		  ENDDO LOOP_i
		
		  !Fix on where the data in log2_Lmk
		  ndata=ndata+1
		
		  !//TO Calculate the standard deviation of Log2(Lmk)
		  log2_Lmk(ndata)=LOG10(temp)/LOG2
		
		  !//OUTPUT each step of m (m=1,2,...,interval(k)) 
		  WRITE(LU+1,7)interval(k),ndata,log2_Lmk(ndata)
		  7 FORMAT(1X,I6,6X,I6,6X,F13.5)

		  !Sum each step of m to calculate the average of Lmk
		  Lmk=Lmk+temp
		
	    ENDDO LOOP_m

        m=size(log2_Lmk)
      
        !//Average of Lmk
	    Lmk=Lmk/m	  
	  
	    !//Calculate standard deviation of Log2(Lmk)
	    CALL std_dev(log2_Lmk,m,std_log2Lmk,status)
	    WRITE(*,*)'SUB:std_dev,status=',status
	    DEALLOCATE(log2_Lmk,STAT=status)
	    WRITE(*,*)'SUB: higuchi1988,Dellocate log2_Lmk status=',status 

	    !!//OUTPUT the Average of Lmk & standard deviation of Log2(Lmk)
	    WRITE(LU,8)interval(k),Lmk,std_log2Lmk
	    8 FORMAT(1X,I6,6X,G15.3,8X,G15.5)

	    !//TO Get the Slope of the Curve: 
	    sum_k=sum_k+LOG10(REAL(interval(k)))/LOG2
	    sum_k2=sum_k2+(LOG10(REAL(interval(k)))/LOG2)**2
	    sum_L=sum_L+LOG10(Lmk)/LOG2
	    sum_L2=sum_L2+(LOG10(Lmk)/LOG2)**2
	    sum_kL=sum_kL+(LOG10(REAL(interval(k)))/LOG2)*(LOG10(Lmk)/LOG2)
      ELSE allocated_okIF
        WRITE(*,*)'In Sub higuchi1988: the array(x or log2_Lmk) is not allocated'
        RETURN
      ENDIF allocated_okIF
	ENDDO LOOP_k
	
	CLOSE(LU+1)
	WRITE(*,9)
	9 FORMAT(1X,'Calculating Curve Length Finished !' )
	
	nvals=REAL(SIZE(interval))      !Total numbers of values in interval
	ave_k=sum_k/nvals
	ave_L=sum_L/nvals
	slope=(sum_kL-sum_k*ave_L)/(sum_k2-sum_k*ave_k)
	intercept=ave_L-slope*ave_k
	rcoef=(nvals*sum_kL-sum_k*sum_L)/   &
	    SQRT((nvals*sum_k2-sum_k**2)*(nvals*sum_L2-sum_L**2))
	
	
	WRITE(LU,10)slope,intercept,rcoef
	10 FORMAT(1X,'%======',6X,'===========',12X,'==========='/,&
			 1X,'% SLOPE',6X,' INTERCEPT',12X,'Relative Coefficiency',/,&
			 1X,F10.4,6X,F10.4,12X,F10.4)
	CLOSE(LU)
	
	DEALLOCATE(interval,STAT=status)
	WRITE(*,*)'SUB: higuchi1988,Dellocate interval status=',status

	WRITE(*,*)one_day_data,'SLOPE=',slope
	WRITE(*,*)'Relative Coefficiency=',rcoef
			
	END SUBROUTINE higuchi1988_v5
	
	!**********************************************
	!****SUBROUTINE std_dev************************
    !**********************************************
    SUBROUTINE std_dev(a,k,sd,error)
    !
    !  Purpose:
    !    To calculate the standard deviation of an array.
    USE  higuchi1988_type
    
    IMPLICIT NONE
  
    INTEGER(I4B),INTENT(IN)::k              !No. of vals in array a.
    REAL(DP),DIMENSION(k),INTENT(IN)::a   !Input data
    REAL(DP),INTENT(OUT)::sd        !Standard deviation
    INTEGER(I1B),INTENT(OUT)::error        !Flag:0--No error
                                    !     1--std invalid
                                    !     2--Error
  
    INTEGER(I4B)::i                !Loop index
    REAL(DP)::sum_x               !Sum of Input values
    REAL(DP)::sum_x2              !Sum of input values squared
  
    !Initialize the sums to zero.
    sum_x=0.
    sum_x2=0.
  
    !Accumulate sums.
    DO i=1,k
      sum_x=sum_x+a(i)
      sum_x2=sum_x2+a(i)**2
    ENDDO
    
    !Check to see if we have enough input data
    IF(k>=2)THEN  !Enough data
      !Calculate Standard deviation
      sd=sqrt((k*sum_x2-sum_x**2)/(k*(k-1)))
      error=0
    ELSEIF(k==1)THEN  !No valid std_dev
      sd=0.
      error=1
    ELSE
      sd=0.
      error=2
    ENDIF
  
    END SUBROUTINE std_dev
    
    function num2str(intNum)
    integer,intent(in) :: intNum
    character(len=:),allocatable:: num2str
    character(len=255):: tmpstr
    write(tmpstr,*)intNum
    allocate(character(len=len_trim(adjustl(tmpstr)))::num2str)
    num2str=trim(adjustl(tmpstr))
    end function num2str
    
    function num2strN(intNum,n)
    integer,intent(in) :: intNum,n
    character(len=:),allocatable:: num2strN
    character(len=255):: tmpstr
    write(tmpstr,1)intNum
1   format(i<n>.<n>) 
    allocate(character(len=len_trim(adjustl(tmpstr)))::num2strN)
    num2strN = trim(adjustl(tmpstr))
    end function num2strN
 
END MODULE higuchi1988_subs
