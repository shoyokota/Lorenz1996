#!/bin/sh
set -ex

I_CMP_NR=1
I_CMP_DA=1
I_RUN_NR=1
I_RUN_DA=1
I_RUN_VE=1
I_RUN_BE=1

#FC=gfortran
#MKL=llapack
FC=ifort
MKL=mkl

if [ ${I_CMP_NR} -eq 1 ]
then  
    ${FC} \
	prm.f90 io.f90 random.f90 lorenz1996.f90 obsope.f90 \
	nature_run.f90 -o nature_run.out
    rm -f *.mod
    echo 'finish compiling nature_run'
fi

if [ ${I_CMP_DA} -eq 1 ]
then  
    ${FC} \
	prm.f90 io.f90 random.f90 lorenz1996.f90 obsope.f90 \
	matrix.f90 minimizer.f90 \
	kf.f90 envar.f90 etlmvar.f90 fdvar.f90 ensrf.f90 \
	da_run.f90 \
	-${MKL} \
	-o da_run.out
    rm -f *.mod
    echo 'finish compiling da_run'
fi

if [ ${I_RUN_NR} -eq 1 ]
then  
    cat <<EOF > namelist.txt
&NAMPRM
nx=40, forcing=8.0, dt=0.01, ft_nature_pre=73.0, ft_nature=20.0,
nslot=1, dt_cycle=0.05,
nobs=40, obsthin=1, obsid=1, obserr=1.0,
!nobs=20, obsthin=2, obsid=1, obserr=1.0,
&END
&NAMFILE
truthfile='TRUTH.txt',
obsfile='OBS.txt',
&END
EOF
    ./nature_run.out < namelist.txt > log_nature.txt || echo 'abort nature_run'
    echo 'finish running nature_run'
fi

if [ ${I_RUN_DA} -eq 1 ]
then  
    for NAME in EnSRF 3DVar 4DVar En3DVar En4DVar ETLMVar HEn3DVar HEn4DVar HETLMVar
    do
	case ${NAME} in
	    "KF")       I_DETDA=1 ; I_ENSDA=0 ; HYB_BETAB=1.0 ; HYB_BETAE=0.0 ;;
	    "EnSRF")    I_DETDA=0 ; I_ENSDA=5 ; HYB_BETAB=0.0 ; HYB_BETAE=1.0 ;;
	    "3DVar")    I_DETDA=2 ; I_ENSDA=0 ; HYB_BETAB=1.0 ; HYB_BETAE=0.0 ;;
	    "4DVar")    I_DETDA=4 ; I_ENSDA=0 ; HYB_BETAB=1.0 ; HYB_BETAE=0.0 ;;
	    "En3DVar")  I_DETDA=2 ; I_ENSDA=5 ; HYB_BETAB=0.0 ; HYB_BETAE=1.0 ;;
	    "En4DVar")  I_DETDA=4 ; I_ENSDA=5 ; HYB_BETAB=0.0 ; HYB_BETAE=1.0 ;;
	    "ETLMVar")  I_DETDA=3 ; I_ENSDA=5 ; HYB_BETAB=1.0 ; HYB_BETAE=0.0 ;;
	    "HEn3DVar") I_DETDA=2 ; I_ENSDA=5 ; HYB_BETAB=0.5 ; HYB_BETAE=0.5 ;;
	    "HEn4DVar") I_DETDA=4 ; I_ENSDA=5 ; HYB_BETAB=0.5 ; HYB_BETAE=0.5 ;;
	    "HETLMVar") I_DETDA=3 ; I_ENSDA=5 ; HYB_BETAB=0.5 ; HYB_BETAE=0.5 ;;
	esac
	cat <<EOF > namelist.txt
&NAMPRM
i_detda=${I_DETDA}, i_ensda=${I_ENSDA},
nx=40, forcing=8.0, dt=0.01, ft_nature_pre=73.0, ft_nature=20.0,
nslot=1, aslot=1, niter=20, dt_cycle=0.05, ft_cycle_pre=20.0, ft_cycle=10.0, ft_max=1.0,
nobs=40, obsthin=1, obsid=1, obserr=1.0,
nmem=3, hyb_betab=${HYB_BETAB}, hyb_betae=${HYB_BETAE}, amp_bcli=1.0, sigma_bcli=2.0,
sigma_sloc=2.0, sigma_tloc=0.2, sigma_sloce=1.0, sigma_tloce=0.2, infl=0.0, rtpp=0.8
!nslot=1, aslot=1, niter=20, dt_cycle=0.05, ft_cycle_pre=20.0, ft_cycle=20.0, ft_max=1.0,
!nobs=20, obsthin=2, obsid=1, obserr=1.0,
!nmem=5, hyb_betab=${HYB_BETAB}, hyb_betae=${HYB_BETAE}, amp_bcli=1.0, sigma_bcli=1.0,
!sigma_sloc=0.5, sigma_tloc=0.2, sigma_sloce=0.5, sigma_tloce=0.2, infl=0.0, rtpp=0.8
&END
&NAMFILE
truthfile='TRUTH.txt',
obsfile='OBS.txt',
detdafile='DATDA.txt',
ensdafile='ENSDA.txt',
biasfile='BIAS.txt',
rmsefile='RMSE.txt',
sprdfile='SPRD.txt',
&END
EOF
	rm -f MAT_${NAME}*.txt
	./da_run.out < namelist.txt > log_da_${NAME}.txt || echo 'abort da_run'
	echo 'finish running da_run'
	mv RMSE.txt RMSE_${NAME}.txt
	mv BIAS.txt BIAS_${NAME}.txt
	mv SPRD.txt SPRD_${NAME}.txt
    done
fi

if [ ${I_RUN_VE} -eq 1 ]
then  
    for NAME in RMSE BIAS SPRD
    do
	gnuplot <<EOF
set term postscript eps enhanced color
set out "${NAME}.eps"
#set term png enhanced
#set out "${NAME}.png"
set grid
set xrange [0:10]
#set yrange [0:1]
set xlabel "time"
set ylabel "${NAME}"
plot \
"${NAME}_EnSRF.txt"    w l lt 1 lw 3 lc rgb "gray"       title "EnSRF",\
"${NAME}_3DVar.txt"    w l lt 1 lw 3 lc rgb "khaki"      title "3DVar",\
"${NAME}_4DVar.txt"    w l lt 1 lw 3 lc rgb "skyblue"    title "4DVar",\
"${NAME}_En3DVar.txt"  w l lt 1 lw 3 lc rgb "green"      title "En3DVar",\
"${NAME}_En4DVar.txt"  w l lt 1 lw 3 lc rgb "blue"       title "En4DVar",\
"${NAME}_ETLMVar.txt"  w l lt 1 lw 3 lc rgb "red"        title "ETLMVar"
EOF
#"${NAME}_HEn3DVar.txt" w l lt 1 lw 1 lc rgb "dark-green" title "Hybrid(3DVar+En3DVar)",\
#"${NAME}_HEn4DVar.txt" w l lt 1 lw 1 lc rgb "dark-blue"  title "Hybrid(4DVar+En4DVar)",\
#"${NAME}_HETLMVar.txt" w l lt 1 lw 1 lc rgb "dark-red"   title "Hybrid(ETLMVar+En3DVar)"
	convert -density 300 ${NAME}.eps ${NAME}.png
    done
fi

if [ ${I_RUN_BE} -eq 1 ]
then
    for NAME in 4DVar_Bc 4DVar_Be 4DVar_MBcM 4DVar_MBeM ETLMVar_MBcM ETLMVar_MBeM
    do
	head -41 MAT_${NAME}.txt > tmp.txt
	gnuplot <<EOF
set term postscript eps enhanced color
set out "MAT_${NAME}.eps"
#set term png enhanced
#set out "MAT_${NAME}.png"
set grid
set pm3d map corners2color c4
set palette defined (-1 "blue", -0.5 "skyblue", 0 "white", 0.5 "yellow", 1 "red")
set size ratio -1
set cbrange [-1:1]
splot "tmp.txt" matrix
EOF
	convert -density 300 MAT_${NAME}.eps MAT_${NAME}.png
	rm -f tmp.txt
    done
fi

exit 0
