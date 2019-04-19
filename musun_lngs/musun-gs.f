*-----------------------------------------------------------------------*
*									*
*	Program: 		MUSUN-GS				*
*	Author: 		Vitaly A. Kudryavtsev			*
*	Institution: 	Department of Physics and Astronomy             *
*                       The University of Sheffield                     *
*	Date: 		March, 2003					*
*									*
*-----------------------------------------------------------------------*
c
c	The code samples single atmospheric muons at the Gran Sasso
c       underground laboratory (taking into account the slant depth distribution)
c	in the energy range E1-E2, zenith angle range theta1-theta2 (0-90 degrees)
c	and azimuthal angle range phi1-phi2 (0-360 degrees).
c       At present only the following ranges of parameters are supported:
c       E1 = 0, E2 = 10^6 GeV, theta1 = 0, theta2 = 90 deg, phi1 = 0, phi2 = 360 deg.
c
c	Program uses muon energy spectra at various depths and zenith
c	angles obtained with MUSIC code for muon propagation and Gaisser's
c	formula for muon spectrum at sea level
c	(T.K.Gaisser, Cosmic Rays and Particle Physics, Cambridge
c	University Press, 1990) modified for large zenith angles and
c	prompt muon flux with parameters from the best fit to the LVD data. 
c
c	Note: the muon spectrum at sea level does not take into account
c	the change of the primary spectrum slope at and above the knee
c	region (3*10^15-10^16 eV).
c
c	Program uses the tables of muon energy spectra at various
c       zenith and azimuthal angles at Gran Sasso
c       calculated with the muon propagation code MUSIC and the
c       angular distribution of muon intensities at LNGS (Hall A) from
c       the best fit to the LVD experimental data.
c
c       Coordinate system is connected to the LVD detector.
c       This is done to allow the simulation of muons on a surface of
c       a rectangular parallelepiped, the planes of which are parallel
c       to the walls of the laboratory halls (and also parallel to the
c       planes of the LVD detector).
c       The zero point is the far right corner of the Hall A as seen from
c       the main road connecting the halls.
c       X-axis is oriented along the short side of the Hall A (to the south-west).
c       Y-axis is oriented along the long side of the Hall A (to the south-east).
c       Z-axis is running up.
c       See Astropart. Phys., vol. 2, p. 103 (1994) for more information.
c       Note that I use here azimuthal angle range from 0 to 360 deg 
c       (unlike the system described in the paper: from -180 to +180 deg)
c       To convert phi in the LVD system to the geographical azimuth, the 
c       following equation should be used: phi_geogr = phi_lvd + 141.6 deg,
c       where geographical azimuth is calculated from north to west.
c       


***************************************************************************
	subroutine initialization
     *  (theta1,theta2,phi1,phi2,igflag,s_hor,s_ver1,s_ver2,FI)
        parameter (pi=3.141592654)
	real*4 em1(121)
	common/sam/spmu(121,62,51),fnmu(32401),idepth(360,91),
	1    fmu(360,91),e1,e2,the1,the2,ph1,ph2
	data fmu/32760*0./,fnmu/32401*0./,spmu/382602*0./
	
	open(unit=3,file='muint-gs.dat',
	1    form='formatted',status='old')
	read(3,104)fmu
	close(3)
 104	format(91(36(10f8.3/)/))
	open(unit=2,file='musp-gs.dat',
	1    form='unformatted',status='old')
	read(2)spmu
	close(2)
	open(unit=1,file='depth-gs.dat',
	1    form='formatted',status='old')
	read(1,103)idepth
 103	format(36(10I7/))
	close(1)
	
	do i=2,121
	em1(i)=0.05*(i-1)
	end do
	em1(1)=alog10(0.105658)

	the1 = theta1
	the2 = theta2
	ph1 = phi1
	ph2 = phi2
	c1=cos(pi/180.*theta1)
	c2=cos(pi/180.*theta2)
	ph1=pi/180.*phi1
	ph2=pi/180.*phi2
	dph=ph2-ph1

	ipc = 1
	theta = theta1
	dc = 1.
	sc = 0.
 1	theta = theta + dc/2.
	theta0 = pi / 180. * theta
	cc = cos(theta0)
	ash = s_hor * cc
	asv01 = s_ver1 * sqrt(1. - cc * cc)
	asv02 = s_ver2 * sqrt(1. - cc * cc)
	ic1 = theta + 0.999
	ic2 = ic1 + 1
	if(ic2.gt.91) ic2 = 91
	if(ic1.lt.1) ic1 = 1
	phi = phi1
	dp = 1.
 2	phi = phi + dp/2.
	phi0 = pi / 180. * phi
	asv1 = asv01 * abs(cos(phi0))
	asv2 = asv02 * abs(sin(phi0))
	asv0 = ash + asv1 + asv2
	fl = 1.
	if(igflag.eq.1) fl = asv0
	ip1 = phi + 0.999
	ip2 = ip1 + 1
	if(ip2.gt.360) ip2 = 1
	if(ip1.lt.1) ip1 = 360
	sp1 = 0.
	do ii = 1, 4
	   iic = (ii - 1) / 2
	   iip = ii - iic * 2 - 1
	   if(fmu(ip1+iip,ic1+iic).lt.0.) then
	      sp1 = sp1 + 10.**fmu(ip1+iip,ic1+iic) / 4.
	   end if
c	   print *,ip1,iip,ic1,iic,fmu(ip1+iip,ic1+iic)
	end do
c	print *,sp1,fl,dp,dc,pi
	sc = sc + sp1 * fl * dp * pi / 180. * 
	1	sin(theta0) * dc * pi / 180.
	ipc = ipc + 1
c	print *,ipc,sc
	fnmu(ipc) = sc
	phi = phi + dp / 2.
	if(phi.lt.phi2-dp/2.) go to 2
	theta = theta + dc / 2.
	if(theta.lt.theta2-dc/2.) go to 1
	FI = sc

	do ipc1 = 1, ipc
	   fnmu(ipc1) = fnmu(ipc1) / fnmu(ipc)
	end do
c	print *, fnmu(1), fnmu(ipc)

	return
	end

 
************************************************************************
	subroutine sampling(E,theta,phi,dep)
	common/sam/spmu(121,62,51),fnmu(32401),idepth(360,91),
	1    fmu(360,91),e1,e2,theta1,theta2,phi1,phi2
        parameter (pi=3.141592654)

	call ranlux(yfl,1)
	i=1
 1	i=i+1
	if(yfl.gt.fnmu(i)) go to 1
c	i = 2
	ic=(i-2)/360
	ip=i-2-ic*360
c	print *,i,ic,ip
	call ranlux(yfl,1)
	theta = theta1 + 1. * (ic + yfl)
 	call ranlux(yfl,1)
	phi = phi1 + 1. * (ip + yfl)
c	print *, theta, phi
	dep = idepth(ip+1,ic+1) * 2.71
c	theta = 0.
c	phi = 0.
c	dep = 3000.
c	print *, dep,theta,phi

	call ranlux(yfl,1)
	ic1 = cos(pi / 180. * theta) / 0.02 + 1.
	if(ic1.lt.1) ic1 = 1
	if(ic1.gt.51) ic1 = 51
	ip1 = dep / 200. - 13 
	if(ip1.lt.1) ip1 = 1
	if(ip1.gt.62) ip1 = 62
	
	k=1
 11	k=k+1
	if(yfl.gt.spmu(k,ip1,ic1)) go to 11
	En1 = 0.05 * (k-2)
	En2 = 0.05 * (k-1)
	call ranlux(yfl,1)
	E = 10.**(En1 + (En2 - En1) * yfl)
c	if(k.eq.2) E = 10.**(alog10(0.105658) + 
c	1    (0.05 - alog10(0.105658)) * yfl)
c	E = 10.**En1
c	print *,k,ic1,ip1,en1,en2,E
	
	return
	end

