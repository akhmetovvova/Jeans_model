#!/bin/bash


#i_th	j_R	k_Z	R	Theta	Z	VR	Vtheta	VZ	sigVR2	sigVth2	sigVRVZ	gradVR_R	grad_sigVR2_R	grad_sigVRVZ_Z	Vcru	density
#1	3	2	0.557731	2.597466	-0.768468	-1.862521	-58.018113	-1.656494	26695.413126	16134.378281	474.174166	-0.754943	0.000000	3.897214	139.666762	0.064090

#./plot-gifs_index_2_data res.txt Vtheta_Vrot.gif '($4):($8)' '($4):($16)' 0 20 180 250 R Vrot 1 0 120 theta 120 1 2 10 100

#./plot-gifs_index_2_data res.txt Vtheta_Vrot_khr_f_woVR.gif '($4):($8)' '($4):($16)' 4 20 180 250 'R, kpc' 'Vcir, km/s' 1 30 90 Vth Vcr 120 1 2 10 100

#cut -d"," -f 5-7,10-12,16-18,63-65 RGBGigants_80.csv > RGBGigants_80_cut.csv

./Jeans_model RGBGigants_80_cut.csv > res_80.txt 2>res_80.err

awk '{ fname="res_z_"$3".txt"; print $0 > fname; }'  res_80.txt

for((i=0;i<5;i++));
do
fname=res_z_$i.txt;
echo $fname;
z=$(echo "scale=3; ($i-2)*0.5" | bc -l | xargs printf "%.1f")
echo $z;
./plot-gifs_index_2_data $fname Vtheta_Vrot_z_$z.gif '($4):($8)' '($4):($16)' 4 20 180 250 'R, kpc' 'V, km/s z= '$z'' 1 30 90 Vth Vcr 120 1 2 10 100
plot-1D_3f_error_point Vcir_R_z_$z.eps Vcru_th_160_z_$z.txt '($1):($19)' Vcru_th_180_z_$z.txt '($1):($19)' '(0):(0):(0)'  Vcru_th_200_z_$z.txt '($1):($19)' '(0):(0):(0)' Eilers_2019_Vrot.dat '($1):($2)' '($1):($2):($3)' 4 20 180 250 'R, kpc' 'V, km/s' 2 10 '@{/Symbol q }  = 160'  '@{/Symbol q }  = 180' '@{/Symbol q }  = 200' 'Eilers2019'
done

#Vcru_th_160_z_0.0.txt
#R	sR	Theta	sTheta	Z	sZ	VR	sVR	Vtheta	sVtheta	VZ	sVZ	sigVR2	s_sigVR2	sigVth2	s_sigVth2	sigVRVZ	s_sigVRVZ	Vcru	sVcru	Nbin

#plot-1D_3f_error_point Vcir_R.eps Vcru_th_160_z_0.0.txt '($1):($19)' Vcru_th_180_z_0.0.txt '($1):($19)' '(0):(0):(0)'  Vcru_th_200_z_0.0.txt '($1):($19)' '(0):(0):(0)' Eilers_2019_Vrot.dat '($1):($2)' '($1):($2):($3)' 4 20 180 250 'R, kpc' 'V, km/s' 2 10 '@{/Symbol q }  = 160'  '@{/Symbol q }  = 180' '@{/Symbol q }  = 200' 'Eilers2019'
#plot-1D_3f_error_point cor_pmdec_gmag.eps cor_pix_$namefile.txt '($22):($5)' 3d_bin_by_Gmag_cor_$namefile '($4):($9)' '(0):(0):(0)' 3d_bin_by_Gmag_S_cor_$namefile '($4):($9)' '($4):($9):($15)'  3d_bin_by_Gmag_N_cor_$namefile '($4):($9)' '($4):($9):($15)'  14 21 -0.5 0.5 Gmag "{/Symbol m_d}, mas/yr" 1 0.1 All S N "Corrected {/Symbol m_d} from Gmag"
