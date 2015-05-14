# Higuchi1988
Calculate the Fractal Dimension D. It has a relation with Detrended Fluctuation Analysis of 3-D.

  Purpose:  Calculating the fractal dimension D
  
  Step1: Get the length of the curve,L(m,k) 
        m: initial time
        k: interval time
        L(m,k)=sum(i=1,[(N-m)/k],|X(m+i*k)-X(m+(i-1)*k)|/k)
	   
	Step2: Get the average of L(m,k)& its Standard_dev
	       ave_Lmk=sum(m=1,k,L(m,k))/k
	       std_dev_Lmk=sqrt(k*L(m,k)**2-sum(L(m,k))**2))/(k*(k-1))
	      
References:

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
