	PROGRAM FractalAnalysis
	!  Purpose:
	!    Calculating the fractal dimension
	!	 Step1: Get the length of the curve,L(m,k) 
	!		m:initial time
	!		k:interval time
	!	  L(m,k)=sum(i=1,[(N-m)/k],|X(m+i*k)-X(m+(i-1)*k)|/k)
	!   
	!    Step2: Get the average of L(m,k)& its Standard_dev
	!       ave_Lmk=sum(m=1,k,L(m,k))/k
	!       std_dev_Lmk=sqrt(k*L(m,k)**2-sum(L(m,k))**2))/(k*(k-1))
	!      
	!  References:
	!	1.Higuchi,T.,1988.Approach to an irregular time 
	!	   series on the basis of fractal theory. 
	!	   Physica D 31.277-283
	!	2.Burlaga,L.F., Klein,L.W.,1986.Fractal structure of
	!	   the interplanetary magnetic field. J.Geophys.Res.91,374-350
	!	3.Gotoh,K.,Hayakawa,M.,Smirnova,N.A.,etal.2004. Fractal 
	!	   analysis of seimogenic ULF emissions. Physis and Chemistry of 
	!	   the Earth.29,419-424

	!  Record of revisions:
	!	Date	Programmer	Description of change
	!	====	==========	=====================
	! 2012/1/5    Q.Wang		Original code
	! 2012/1/8	  Q.Wang		Revised
	! 2012/1/9	  Q.Wang		Revised
	! 01/14/12    Q.Wang        Revised
    
    USE higuchi1988_type
    USE higuchi1988_subs
    use dflib
    
  	IMPLICIT NONE
  	
    !> for DFA
	CHARACTER(len=30)::filename,lmk
	CHARACTER(len=4)::nof  !No. of file
	INTEGER(I4B)::nf       !No. of file
	REAL(DP),SAVE,ALLOCATABLE,DIMENSION(:) :: x
    REAL(DP),SAVE,allocatable,DIMENSION(:) :: section1year
    integer,parameter:: Nsubsec=180
	INTEGER(I1B)::fileindex
	INTEGER(I4B)::i, j, Nsec!Loop Index
	REAL(DP)::slope
	REAL(DP)::rcoef
	INTEGER(I1B)::status
    
    !> for cmd.exe
    character(len=255) :: cmd
    logical :: res

	!> input file's name
    WRITE(*,1)
1   FORMAT("Please Enter Other File's filename: ",$)
    READ(*,*)filename
    !filename='JIHminSr160'

	  
	!  Open a file to output everyday's slope
	OPEN(100,FILE=trim(filename)//'slopes'//num2str(Nsubsec)//'.txt',STATUS='REPLACE',ACTION='WRITE')
	
	  !  Open File & Read data
	  CALL check_in(filename,x,status)
      
	
	  IF(ALLOCATED(x))THEN
          allocate(section1year(Nsubsec))
          Nsec=floor((size(x)-size(section1year)+10)/10.0)
         write(*,*)'Total section: ', Nsec
          do i=1,floor((size(x)-size(section1year)+10)/10.0)
              
              section1year(1:Nsubsec) = x(1 + (i-1)*10 : NsubSec + (i-1)*10 )
              
              !>> output subsection's data into files
              open(90,file=trim(filename)//'sec'//num2strN(i,len_trim(num2str(Nsec)))//'.dat',status='replace',action='write')
                    do j=1,Nsubsec
                        write(90,*)section1year(j)
                    enddo
              close(90)
              !> mkdir subdataBJIminSr160
              cmd='mkdir subdata_'//trim(filename)//num2str(Nsubsec)
              res = systemqq(cmd)      
              !> move LMKBJIminSr160section??? into LMK_BJIminSr160
              cmd='move '//trim(filename)//'sec*.dat '//'subdata_'//trim(filename)//num2str(Nsubsec)
              res = systemqq(cmd)
      
              !>> DFA process
              CALL higuchi1988_v5(trim(filename)//'section'//num2str(i),section1year,lmk,slope,rcoef,status)
              WRITE(*,*)'STATUS=',status
              WRITE(100,10)ABS(slope),abs(rcoef)
10            FORMAT(1X,2F16.5)
              
              !> write xTickDays file for GMT5 plotting
              !if(mod(i,9)==0) write(10,*)i, 'a', (i-1)*10+1
              
          enddo
          deallocate(section1year)
          close(100)

          IF(status/=0)STOP
        
      ELSE
          WRITE(*,*)'Sub chek_in: the array is not allocated.'        
      ENDIF
      
      
      !> mkdir LMK_BJIminSr160
      cmd='mkdir LMK_'//trim(filename)//num2str(Nsubsec)
      res = systemqq(cmd)      
      !> move LMKBJIminSr160section??? into LMK_BJIminSr160
      cmd='move LMK'//trim(filename)//'section* '//'LMK_'//trim(filename)//num2str(Nsubsec)
      res = systemqq(cmd)
      
      
      
      
   
	END PROGRAM FractalAnalysis
