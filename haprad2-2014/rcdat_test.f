      program rcdat
      implicit none
      include 'haprad_consts.inc'
	integer i,j
	double precision Ebeam,x,q2,z,t,phi,nu
      double precision sib,sig,delta,tail,tai2
      double precision mn,DSIampS,rc,s
      double precision phi_min,phi_max,t_min,t_max,z_min,z_max,m2th
      double precision x_min,x_max
      data mn/0.93956536d0/
	double precision xmas(16),q2mas(16),zmas(16),ptmas(16),mxmas(16),del(16,4)
	data xmas/0.2645,0.2645,0.31,0.37,0.27,0.27,0.27,0.27,0.18,0.24,0.31,0.37,0.26,0.25,0.23,0.20/
	data q2mas/1.913,1.913,1.6,2.0,1.46,1.44,1.44,1.43,1.1,1.4,1.7,2.0,1.44,1.41,1.37,1.26/
	data zmas/0.6481,0.6481,0.61,0.61,0.54,0.61,0.69,0.77,0.58,0.57,0.56,0.56,0.54,0.61,0.69,0.77/
	data ptmas/0.3046,0.8,0.41,0.39,0.43,0.43,0.42,0.36,0.43,0.36,0.32,0.28,0.38,0.34,0.30,0.24/
	data mxmas/1.596,1.596,1.1,1.1,1.1,1.1,1.1,1.1,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4/
	
	amhh=0.1395675d0  ! mass of detected hadron
	amhu=0.93956536d0 ! mass of undetected hadron for exclusive radiative tail 
	m2th=(amp+ampi)**2 ! minimal invariant mass square of undetected hadron
c	m2th=1.4d0**2     ! another limit on  minimal invariant mass square of undetected hadron 
       Ebeam=5.498d0         ! energy of lepton beam
       s=2d0*amp*Ebeam
	open(12,file='test.dat')
	write(12,'(6a12)')'Energy ','x      ','q2     ','z      ','pt      ','phi     '
	do j=1,96
	phi=0.0327249+dble(j-1)*0.0654498 ! must be consistent with line above (96)
	do i=2,2 !1,16
	print*,j,i
      call fhaprad(Ebeam,xmas(i),q2mas(i),zmas(i),ptmas(i),phi,mxmas(i),sib,sig,delta,tail)
	write(*,'(a5,12g12.4)') 'aaaaa',xmas(i),q2mas(i),zmas(i),ptmas(i),phi,sib,sig/sib,tail/sib,delta
	write(12,'(11g12.4)')Ebeam, xmas(i),q2mas(i),zmas(i),ptmas(i),phi,sib,sig/sib,tail/sib,delta
	enddo
	enddo
      end
