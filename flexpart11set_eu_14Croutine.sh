#!/bin/bash 
### create pathnames, RELEASES, COMMAND file for flexpart
### script should be in parent folder of options
### run : bash flexpart11set.sh 11 inputs
### zhendong.wu@nateko.lu.se

HOME_PATH=$(pwd)
OPTION_PATH=$HOME_PATH/options

LIFETIME=10 # in days

echo ------ make pathnames, RELEASES, COMMAND, OUTGRID, OUTGRID_NEST, AGECLASSES file ------
echo

# ======================================
# arguments
# ======================================
simulationidproject=$1 # simulations id and project name, e.g. htm150C14eu
input=$2 # meteo input path
output=$3 # flexpart output path
id=$4 # species id, NOTE one id for now 
start=$5 # start date, e.g. YYYYMMDD
end=$6 # end date, e.g. YYYYMMDD
step=$7 # time step in hour
lon=$8 # longitude of release box -180 < LON1 <180
lat=$9 # latitude of release box, -90 < LAT1 < 90
z=${10} # height of release
particles=${11} # Total number of particles to be released

output=${output}${start}/
mkdir -p ${output}
WORK_PATH=${HOME_PATH}/${simulationidproject}/${start}
mkdir -p ${WORK_PATH}
# ======================================
# backup orignal files
# ======================================
cp -r ${OPTION_PATH} ${WORK_PATH}
# cp -n ${OPTION_PATH}/RELEASES ${OPTION_PATH}/RELEASES_backup
# cp -n ${OPTION_PATH}/COMMAND ${OPTION_PATH}/COMMAND_backup
# cp -n ${OPTION_PATH}/OUTGRID ${OPTION_PATH}/OUTGRID_backup
# cp -n ${OPTION_PATH}/OUTGRID_NEST ${OPTION_PATH}/OUTGRID_NEST_backup

# ======================================
# make pathnames file
# ======================================
cat <<EOM >${WORK_PATH}/pathnames
${WORK_PATH}/options/
${output}
${input}
${input}AVAILABLE
EOM
echo pathnames is done, under ${WORK_PATH}/
echo

# ======================================
# make AGECLASSES file
# ======================================
cat <<EOM >${WORK_PATH}/options/AGECLASSES
************************************************
*                                              *
*Lagrangian particle dispersion model FLEXPART *
*         Please select your options           *
*                                              *
*This file determines the ageclasses to be used*
*                                              *
*Ages are given in seconds. The first class    *
*starts at age zero and goes up to the first   *
*age specified. The last age gives the maximum *
*time a particle is carried in the simulation. *
*                                              *
************************************************
&NAGE
 NAGECLASS=          1, ! Number of ageclasses for the age spectra calculation
 /

&AGECLASS
 LAGE=    $(( $LIFETIME * 24 * 60 * 60 )),      ! Maximum age of particles in seconds for each ageclass
 /
EOM

# ======================================
# make RELEASES file
# ======================================
cat <<EOM >${WORK_PATH}/options/RELEASES
***************************************************************************************************************
*                                                                                                             *
*                                                                                                             *
*                                                                                                             *
*   Input file for the Lagrangian particle dispersion model FLEXPART                                          *
*                        Please select your options                                                           *
*                                                                                                             *
*                                                                                                             *
*                                                                                                             *
***************************************************************************************************************
&RELEASES_CTRL
 NSPEC      =           1, ! Total number of species
 SPECNUM_REL=          ${id}, ! Species numbers in directory SPECIES
 /
EOM

start_release=$(date -u -d "$start")
end_release=$(date -u -d "$end -$step hours")

i=0
releaseno=0
while [ "$thedate" != "$start_release" ]
do
# 	echo $thedate
	IDATE1=$(date -u -d "$end -$i hours" +%Y%m%d)
	ITIME1=$(date -u -d "$end -$i hours" +%H%M%S)

	IDATE2=$(date -u -d "$end -$i hours" +%Y%m%d)
	ITIME2=$(date -u -d "$end -$i hours" +%H%M%S)

	# Release end date, same as IDATE1
	i=$(( i - $step ))

 	thedate=$(date -u -d "$end -$i hours")
	
cat <<EOM >>${WORK_PATH}/options/RELEASES
&RELEASE                   ! For each release 
IDATE1  =       $IDATE1, ! Release start date, YYYYMMDD: YYYY=year, MM=month, DD=day
ITIME1  =         $ITIME1, ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds
IDATE2  =       $IDATE2, ! Release end date, same as IDATE1
ITIME2  =         $ITIME2, ! Release end time, same as ITIME1
LON1    =          $lon, ! Left longitude of release box -180 < LON1 <180
LON2    =          $lon, ! Right longitude of release box, same as LON1
LAT1    =         $lat, ! Lower latitude of release box, -90 < LAT1 < 90
LAT2    =         $lat, ! Upper latitude of release box same format as LAT1 
Z1      =         $z, ! Lower height of release box meters/hPa above reference level
Z2      =         $z, ! Upper height of release box meters/hPa above reference level
ZKIND   =              1, ! Reference level 1=above ground, 2=above sea level, 3 for pressure in hPa
MASS    =       1.0000E0, ! Total mass emitted, only relevant for fwd simulations
PARTS   =          $particles, ! Total number of particles to be released
COMMENT =    "RELEASE $releaseno $IDATE1:$ITIME1", ! Comment, written in the outputfile
/
EOM
	
done
echo RELEASES file is done, under ${WORK_PATH}/options/. i.e. $((- $i / $step  )) releases from $start to $end. 
echo



# ======================================
# make COMMAND file
# ======================================
IBDATE=$(date -u -d "$start -$LIFETIME days" +%Y%m%d)
IBTIME=$(date -u -d "$start -$LIFETIME days" +%H%M%S)
IEDATE=$(date -u -d "$end" +%Y%m%d)
IETIME=$(date -u -d "$end" +%H%M%S)
OUTSTEP=$(( $step * 60 * 60 )) # in seconds

cat <<EOM >${WORK_PATH}/options/COMMAND
***************************************************************************************************************
*                                                                                                             *
*      Input file for the Lagrangian particle dispersion model FLEXPART                                       *
*                           Please select your options                                                        *
*                                                                                                             *
***************************************************************************************************************
&COMMAND
 LDIRECT=               -1, ! Simulation direction in time   ; 1 (forward) or -1 (backward)
 IBDATE=         $IBDATE, ! Start date of the simulation   ; YYYYMMDD: YYYY=year, MM=month, DD=day  
 IBTIME=           $IBTIME, ! Start time of the simulation   ; HHMISS: HH=hours, MI=min, SS=sec; UTC
 IEDATE=         $IEDATE, ! End date of the simulation     ; same format as IBDATE 
 IETIME=           $IETIME, ! End  time of the simulation    ; same format as IBTIME
 LOUTSTEP=           $OUTSTEP, ! Interval of model output; average concentrations calculated every LOUTSTEP (s)  
 LOUTAVER=           $OUTSTEP, ! Interval of output averaging (s)
 LOUTSAMPLE=          900, ! Interval of output sampling  (s), higher stat. accuracy with shorter intervals
 LOUTRESTART=       -1, ! Interval of writing restart files (s), switched off when set to -1
 LRECOUTSTEP=        3600, ! Interval of model output at receptors (s)
 LRECOUTAVER=        3600, ! Interval of receptor output averaging (s)
 LRECOUTSAMPLE=      1200, ! Interval of receptor output sampling (s)
 LSYNCTIME=           900, ! All processes are synchronized to this time interval (s)
 CTL=          -5.0000000, ! CTL>1, ABL time step = (Lagrangian timescale (TL))/CTL, uses LSYNCTIME if CTL<0
 IFINE=                 4, ! Reduction for time step in vertical transport, used only if CTL>1 
 IOUT=                  13, ! Gridded output type: [0]off [1]mass 2]pptv 3]1&2 4]plume 5]1&4, +8 for NetCDF output     
 IPOUT=                 1, ! Particle position output: 0]off 1]every output 2]only at end
 LSUBGRID=              1, ! Increase of ABL heights due to sub-grid scale orographic variations;[0]off 1]on 
 LCONVECTION=           1, ! Switch for convection parameterization;0]off [1]on
 LTURBULENCE=           1, ! Switch for turbulence parameterisation;0]off [1]on
 LTURBULENCE_MESO=      0, ! Switch for mesoscale turbulence parameterisation;0]off (recommended) [1]on
 LAGESPECTRA=           1, ! Switch for calculation of age spectra (needs AGECLASSES);[0]off 1]on
 IPIN=                  0, ! Warm start from particle dump; [0]no 1]from restart.bin file 2]from previous partoutput file 3]self made initial conditions 4]restart.bin and self made initial conditions
 IOUTPUTFOREACHRELEASE= 1, ! Separate output fields for each location in the RELEASE file; [0]no 1]yes 
 IFLUX=                 0, ! Output of mass fluxes through output grid box boundaries
 MDOMAINFILL=           0, ! Switch for domain-filling, if limited-area particles generated at boundary
 IND_SOURCE=            1, ! Unit to be used at the source; [1]mass 2]mass mixing ratio 
 IND_RECEPTOR=          2, ! Unit to be used at the receptor; [0]no receptor [1]mass 2]mass mixing ratio 3]wet depo. 4]dry depo.
 MQUASILAG=             0, ! Quasi-Lagrangian mode to track individual numbered particles 
 NESTED_OUTPUT=         0, ! Output also for a nested domain 
 LNETCDFOUT=            1, ! Gridded netcdf output: [0]no [1]yes
 LINIT_COND=            0, ! Output sensitivity to initial conditions (bkw mode only) [0]off 1]conc 2]mmr 
 SFC_ONLY=              0, ! Output only for the lowest model layer, used w/ LINIT_COND=1 or 2
 CBLFLAG=               0, ! Skewed, not Gaussian turbulence in the convective ABL, need large CTL and IFINE
 OHFIELDS_PATH= "../../flexin/", ! Default path for OH file
 NXSHIFT=             359, ! Shift of the global meteorological data. Default 359 for ECMWF and 0 for GFS if not given
 MAXTHREADGRID=         1, ! Set maximum number of threads for doing grid computations. Recommended to set this no higher than 16. High numbers create more overhead and a larger memory footprint, 1=no parallelisation on grid.
 MAXFILESIZE=       10000, ! Maximum output of each partoutput NetCDF-4 file in Mb before a new one is created
 LOGVERTINTERP=         0, ! Flag to set all vertical interpolation to logarithmic instead of linear
 LCMOUTPUT=             0, ! Switch for the Linear Chemistry Module; [0] off [1] on
 /
EOM

echo COMMAND file is done, under ${WORK_PATH}/options/.  
echo


# ======================================
# make OUTGRID file
# ======================================
cat <<EOM >${WORK_PATH}/options/OUTGRID
!*******************************************************************************
!                                                                              *
!      Input file for the Lagrangian particle dispersion model FLEXPART        *
!                       Please specify your output grid                        *
!                                                                              *
! OUTLON0    = GEOGRAPHYICAL LONGITUDE OF LOWER LEFT CORNER OF OUTPUT GRID     *
! OUTLAT0    = GEOGRAPHYICAL LATITUDE OF LOWER LEFT CORNER OF OUTPUT GRID      *
! NUMXGRID   = NUMBER OF GRID POINTS IN X DIRECTION (= No. of cells + 1)       *
! NUMYGRID   = NUMBER OF GRID POINTS IN Y DIRECTION (= No. of cells + 1)       *
! DXOUT      = GRID DISTANCE IN X DIRECTION                                    *
! DYOUN      = GRID DISTANCE IN Y DIRECTION                                    *
! OUTHEIGHTS = HEIGHT OF LEVELS (UPPER BOUNDARY)                               *
!*******************************************************************************
&OUTGRID
 OUTLON0=    -15.00,
 OUTLAT0=     33.00,
 NUMXGRID=       500,
 NUMYGRID=       400,
 DXOUT=        0.1,
 DYOUT=        0.1,
 OUTHEIGHTS=  50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0,
 /
EOM

echo OUTGRID file is done, under ${WORK_PATH}/options/.
echo



# cd ${WORK_PATH}
# exec bash
