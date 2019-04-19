*-----------------------------------------------------------------------*
*									*
*	Program: 		MUSUN-GS				*
*	Author: 		Vitaly A. Kudryavtsev			*
*	Institution: 	Department of Physics and Astronomy             *
*                       The University of Sheffield                     *
*	Date: 		March, 2003, updated 2012, March 2014		*
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
c       Z-axis is pointing up.
c       See Astropart. Phys., vol. 2, p. 103 (1994) for more information.
c       Note that I use here azimuthal angle range from 0 to 360 deg 
c       (unlike the system described in the paper: from -180 to +180 deg)
c       To convert phi in the LVD system to the geographical azimuth, the 
c       following equation should be used: phi_geogr = phi_lvd + 141.6 deg,
c       where geographical azimuth is calculated from north to west.
c       


*************************************************************************
c
c       This is a test program which reads spectra and muon intensities
c       at the LNGS site (subroutine initialization) and samples a number of muons
c       according to their spectra and angular distribution (subroutine sampling).
c       Input and parameters are described below.
c       Subroutine initialization returns the absolute muon intensity for
c       the specified ranges of energies and angles, which can be used
c       for normalisation. The units are described in the format statement.

	integer*4 ipar,ierr,isl
	real*4 depth,E1,E2,theta1,theta2,phi1,phi2,rc

c	Input parameters:

	E1=1.		!Muon energy range (GeV): E1-E2
	E2=1.e6		!(E2>=E1)
	theta1=0.	!Range of zenith angles (degrees): theta1-theta2
	theta2=90.	!(theta2>=theta1)
	phi1=0.		!Range of azimuthal angles (degrees): phi1-phi2
	phi2=360.	!(phi2>=phi1)


	iranlux=101
	call rluxgo(3,iranlux,0,0)

c       Muons are sampled on the surface of a sphere with a unit area
c       perpendicular to the muon flux.
c       Zenith and azimuthal angles are sampled from the slant depth
c       distribution (muon intensity distribution) without taking
c       into account the detector or cavern geometries.
c       This is the case for Borexino simulations.

c	igflag = 0
c        s_hor=1.
c        s_ver1=1.
c        s_ver2=1.
	
c       Alternative way is to sample muons on a surface of a parallelepiped.
c       If you want muons to be sampled on the surface of a parallelepiped
c       decomment the following 4 lines and specify the size (the areas of
c       the planes) of a parallelepiped.

        igflag = 1

	zh1 = -250.	  ! You have to change all these dimensions if you change
        zh2 = +250.	  ! the size of the cuboid. The dimensions are given in cm.
        xv1 = -250.      ! The example shown here is for the lab+rock having the size of
        xv2 = +250.	  ! 40 x 20 x 35 m^3 in x, y and z, respectively. 
        yv1 = -250.	  ! The centre of the coordinate system is the centre of the
        yv2 = +250.	  ! lab. 
			  ! The volume of the cuboid for muon sampling extends 
	                  ! to the rock on all sides.

	s_hor=(xv2-xv1)*(yv2-yv1) ! area of the horizontal plane of the parallelepiped
        s_ver1=(yv2-yv1)*(zh2-zh1)   ! area of the vertical plane of the parallelepiped, 
                                     ! perpendicular to x-axis
        s_ver2=(xv2-xv1)*(zh2-zh1)   ! area of the vertical plane of the parallelepiped,
                                     ! perpendicular to y-axis

	call initialization
     *  (theta1,theta2,phi1,phi2,igflag,s_hor,s_ver1,s_ver2,FI)

	write(6,1004)E1,E2,theta1,theta2,phi1,phi2,FI
 1004	format(2x,'Material - Gran Sasso rock'/
     *  2x,'Density = 2.71 g/cm^3'/
     *  2x,'Parameters for muon spectrum are from the LVD best fit'/
     *  2x,'Muon energy range = ',e11.4,' - ',e11.4,' GeV'/
     *  2x,'Zenith angle range = ',f6.2,' - ',f6.2,' degrees'/
     *  2x,'Azimuthal angle range = ',f6.2,' - ',f6.2,' degrees'/
     *  2x,'Global intensity = ', e11.4,2x,
     *  '(cm^2 s)^(-1) or s^(-1)
     *  (for sampling muons on the surface of the cuboid)')
	
	nmumax=1000000
	se=0.
	st=0.
	sp=0.
	sd=0.
	sse = 0.
	sst = 0.
	ssp = 0.
	ssd = 0.
	i = 1

        open(unit=1,file='muons_gs_01.dat',
     *       form='formatted',status='unknown')

	do nmu=1,nmumax
	   call sampling(E,theta,phi,dep)
	theta1 = theta
	phi1 = phi
        theta = theta * 3.141592654 / 180.

	if(phi.gt.360.) phi = phi - 360.
	if(phi.lt.0) phi = phi + 360.
        phi = phi * 3.141592654 / 180.

*       Muon sign (particle ID - GEANT3 definition)

	call ranlux(yfl,1)
        id_part = 10                ! positive muon
        if(yfl.lt.1./2.38) id_part = 11   !mu+/mu- ration is 1.38; change this if necessary.

*       Direction cosines (momenta projections - normalised to the unit vector)

      cz = -cos(theta)
*       The minus sign above is for z-axis pointing up, so the z-momentum 
*       is pointing down. 
      cx = -sin(theta)*cos(phi)
      cy = -sin(theta)*sin(phi)

*       Muon coordinates

      sh1 = s_hor * cos(theta)
      sv1 = s_ver1 * sin(theta) * abs(cos(phi))
      sv2 = s_ver2 * sin(theta) * abs(sin(phi))
      ss = sh1 + sv1 + sv2
	call ranlux(yfl1,1)
        if(yfl1.le.sh1/ss) then
           z0 = zh2
           call ranlux(yfl,1)
           x0 = (xv2 - xv1) * yfl + xv1
           call ranlux(yfl,1)
           y0 = (yv2 - yv1) * yfl + yv1
        else if(yfl1.le.(sh1+sv1)/ss) then
           if(cx.ge.0.) x0 = xv1
           if(cx.lt.0.) x0 = xv2
           call ranlux(yfl,1)
           z0 = (zh2 - zh1) * yfl + zh1
           call ranlux(yfl,1)
           y0 = (yv2 - yv1) * yfl + yv1
        else
           if(cy.ge.0.) y0 = yv1
           if(cy.lt.0.) y0 = yv2
           call ranlux(yfl,1)
           z0 = (zh2 - zh1) * yfl + zh1
           call ranlux(yfl,1)
           x0 = (xv2 - xv1) * yfl + xv1
        end if

c       Parameters written to the file muons_gs_test.dat
c       nmu - muon number
c       id_part - muon id (positive or negative according to GEANT3 definitions)
c       E - muon energy in GeV assuming ultrarelativistic muons: E_total = E_kinetic
c       x0, y0, z0 - muon coordinates on the surface of parallelepiped
c                    size: 40 m * 20 m * 35 m with the centre in the middle of the lab;
c                    x-axis and y-axis are pointing in the directions described by the
c                    LVD coordinate system, see above.
c                    z-axis is directed upwards.
c       cx, cy, cz - direction cosines.

        write(1,1)nmu,id_part,E,x0,y0,z0,cx,cy,cz
 1      format(i9,i3,f10.1,3f8.1,3f10.6)

	se=se+E/nmumax*10.
	st=st+theta1/nmumax*10.
	sp=sp+phi1/nmumax*10.
	sd=sd+dep/nmumax*10.
	if(nmu.eq.i*nmumax/10) then
	  i = i + 1
	  print *, se, st, sp, sd
	  sse = sse + se
	  sst = sst + st
	  ssp = ssp + sp
	  ssd = ssd + sd
	  se = 0.
	  st = 0.
	  sp = 0.
	  sd = 0.
	end if
	end do

	close(1)

	write(6,1005)nmumax,sse/10.,sst/10.,ssp/10.,ssd/10.
 1005	format(2x,'Number of muons = ',i12/
     *  2x,'Mean muon energy = ',e11.4,' GeV'/
     *  2x,'Mean zenith angle = ',f6.2,' degrees'/
     *  2x,'Mean azimuthal angle = ',f6.2,' degrees'/
     *  2x,'Mean slant depth = ',f10.1,' m w.e.'/)

	stop
	end
