/*
 * main.cpp
 *
 *  Created on: Nov 13, 2025
 *      Author: vova
 */


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define rad2deg 180.0/3.14159265358979323846
#define deg2rad 3.14159265358979323846/180.0

//double theta = -(atan2(y,(par.R0-x)));
//	double theta0 = M_PI+theta;
//	double Ro = sqrt(pow(par.R0-x,2)+y*y);
//	double newX,newY,newVx,newVy;

//rotation coordinate system
//						newX = (xx*cos(theta)+yy*sin(theta));
//						newY = (yy*cos(theta)-xx*sin(theta));
//
//						newVx = (cell[i][j][k].s.vx*cos(theta)+cell[i][j][k].s.vy*sin(theta));
//						newVy = (cell[i][j][k].s.vy*cos(theta)-cell[i][j][k].s.vx*sin(theta));


//Vx = w2*z;
//Vz = -w2*x;

struct data
{
	double R,theta,Z;
	double VR,Vtheta,VZ;

	double sigVR2,sigVth2,sigVZ2;
	double sigVRVth,sigVRVZ,sigVthVZ;

	double den;
	double grad_VR_R;
	double grad_sigVR2_R;
	double grad_sigVRVZ_Z;

	double grad_den_R;
	double grad_sigVR2_VR2_R;
	double grad_den_sigVRVZ_Z;

	double Vcru;
	double Vcru_my;

	int ok;
	int nR;
	int nZ;
};

void skipline(FILE * fr)
{
	while(!feof(fr))
		if(fgetc(fr)=='\n')
			break;
}

int get_index(double z, double max, double min, int n)
{
	if(z==max ){
		return n-1;
	}else if (z==min){
		return 0;
	}else if(z<min || z>max){
		//fprintf(stderr,"Invalid calc index!!!\n");
		return -1;
	} else	{return (int)((z-min)/(max-min)*(n-1));}
}

double calc_den(double R, double Z, double vo, double Rexp, double h_z)
{
	return vo*exp(-R/Rexp - fabs(Z)/h_z);
}

data push_obj( double avgZ,double  avgR, double avgTh, double avgVR,double avgVTh,double avgVZ,double sigVRVZ,double sigVR2,double sigVTh2, double h_R, double h_Z)
{
	data obj;

	obj.R = avgR;
	obj.Z = avgZ;
	obj.theta = avgTh;

	obj.VR = avgVR;
	obj.VZ = avgVZ;
	obj.Vtheta = avgVTh;

	obj.sigVR2 = sigVR2;
	obj.sigVth2 = sigVTh2;
	obj.sigVRVZ = sigVRVZ;

	obj.den = calc_den(avgR, avgZ, 1.0, h_R, h_Z);

	obj.ok = 1;

	return obj;
}

double fun_sigVth2 (double x)
{
//=======================            ==========================
//	a0              = -529.955         +/- 5.502        (1.038%)
//	a1              = 13425.3          +/- 17.6         (0.1311%)
//f(x)=a0+a1/x
	//return -529.955  + 13425.3/x;
	double f = (x-20)*7.5;
	return (-529.955  + 13425.3/x) + f*f;
}

double fun_sigVR2_VR2 (double x)
{
//=======================            ==========================
//a0              = 5.97659          +/- 4.765e+10    (7.973e+11%)
//a1              = -0.309958        +/- 0.0008517    (0.2748%)
//a2              = 3.22338          +/- 7.961e+09    (2.47e+11%)
//a3              = 25.7774          +/- 0.07164      (0.2779%)
	return pow(5.97659*exp(-0.3099586*x+3.22338) + 25.7774, 2);
}


double fun_sigVR2 (double x)
{
//	a0              = 5.98654          +/- 5.842e+10    (9.759e+11%)
//	a1              = -0.316276        +/- 0.001045     (0.3303%)
//	a2              = 3.24087          +/- 9.743e+09    (3.006e+11%)
//	a3              = 25.3984          +/- 0.07169      (0.2823%)


	return pow(5.98654*exp(-0.316276*x+3.24087) + 25.3984,2);
}


void calc_gradient_my(data *** obj, int ii, int jj, int kk, int size_i,int size_j, int size_k, int Nbin_i, int Nbin_j, int Nbin_k)
{
	int imin = ii-size_i;
	if(imin<0){imin=0;}
	int imax = ii+size_i;
	if(imax>=Nbin_i){imax=Nbin_i-1;}

	int jmin = jj-size_j;
	if(jmin<0){jmin=0;}
	int jmax = jj+size_j;
	if(jmax>=Nbin_j){jmax=Nbin_j-1;}

	int kmin = kk-size_k;
	if(kmin<0){kmin=0;}
	int kmax = kk+size_k;
	if(kmax>=Nbin_k){kmax=Nbin_k-1;}

	double S_R=0;
	double S_d=0;
	double S_RR=0;
	double S_Rd=0;
	double S_vR =0;
	double S_RvR =0;

	double Zx=0;
	double Zy=0;
	double Zxx=0;
	double Zxy=0;

	double R,d,vR;
	double Z,yZ;
	int nR=0;
	int nZ=0;

	//double Vcru;

	//double R0 = obj[ii][jj][kk].R;

//	by plane in TH_R
	for(int i = imin; i<=imax; ++i)
		for(int j = jmin; j<=jmax; ++j)
		{
			if(obj[i][j][kk].ok==1)
			{
				//fprintf(stderr,"R=%lf den = %lf sigVR2 = %lf\n",obj[i][j][kk].R,obj[i][j][kk].den,obj[i][j][kk].sigVR2);
				R = obj[i][j][kk].R;
				d = obj[i][j][kk].den;
//				vR = obj[i][j][kk].VR;
				vR = obj[i][j][kk].VR*obj[i][j][kk].VR+obj[i][j][kk].sigVR2;

				if(d>0)
				{
					d= log(d);
					S_R += R;
					S_RR+= R*R;
					S_d += d;
					S_Rd+= R*d;
					S_vR += vR;
					S_RvR += R* vR;

					nR++;
				}
			}
		}

	obj[ii][jj][kk].grad_den_R = (nR*S_Rd - S_R*S_d)/(nR*S_RR-S_R*S_R);
	//obj[ii][jj][kk].grad_den_R = 0;
	obj[ii][jj][kk].grad_VR_R = (nR*S_RvR-S_R*S_vR)/(nR*S_RR-S_R*S_R);

	double dR = 0.1;
	obj[ii][jj][kk].grad_sigVR2_VR2_R = ( fun_sigVR2_VR2(obj[ii][jj][kk].R+dR/2.0)-fun_sigVR2_VR2(obj[ii][jj][kk].R-dR/2.0)) / dR;

	obj[ii][jj][kk].nR = nR;

//	double Z0 = obj[ii][jj][kk].Z;

	for(int i = imin; i<=imax; ++i)
		for(int k = kmin; k<=kmax; ++k)
		{
			if( obj[i][jj][k].ok==1)
			{
			Z = obj[i][jj][k].Z;
			yZ = obj[i][jj][k].den*obj[i][jj][k].sigVRVZ;

			Zx += Z;
			Zxx+= Z*Z;
			Zy += yZ;
			Zxy+= Z*yZ;

			nZ++;
			}
		}

	obj[ii][jj][kk].grad_den_sigVRVZ_Z = (nZ*Zxy-Zx*Zy)/(nZ*Zxx-Zx*Zx);
	obj[ii][jj][kk].nZ = nZ;


//	obj[ii][jj][kk].Vcru = obj[ii][jj][kk].Vtheta*obj[ii][jj][kk].Vtheta + obj[ii][jj][kk].sigVth2
//			- (obj[ii][jj][kk].VR*obj[ii][jj][kk].VR + obj[ii][jj][kk].sigVR2)*(1.0 + obj[ii][jj][kk].grad_VR_R) - obj[ii][jj][kk].R * obj[ii][jj][kk].grad_sigVRVZ_Z / obj[ii][jj][kk].den;


	obj[ii][jj][kk].Vcru_my = obj[ii][jj][kk].Vtheta*obj[ii][jj][kk].Vtheta + obj[ii][jj][kk].sigVth2 - (obj[ii][jj][kk].VR*obj[ii][jj][kk].VR + obj[ii][jj][kk].sigVR2)
	- obj[ii][jj][kk].R*obj[ii][jj][kk].grad_sigVR2_VR2_R - obj[ii][jj][kk].R*(obj[ii][jj][kk].VR*obj[ii][jj][kk].VR+obj[ii][jj][kk].sigVR2)* obj[ii][jj][kk].grad_den_R - obj[ii][jj][kk].R * obj[ii][jj][kk].grad_den_sigVRVZ_Z / obj[ii][jj][kk].den	;

	if(obj[ii][jj][kk].Vcru_my>0)
	{
		obj[ii][jj][kk].Vcru_my=sqrt(obj[ii][jj][kk].Vcru_my);
	}else
	{
		obj[ii][jj][kk].Vcru_my=0;
	}

}

void calc_gradient_khrob(data *** obj, int ii, int jj, int kk, int size_i,int size_j, int size_k, int Nbin_i, int Nbin_j, int Nbin_k, double h_R, double h_Z)
{
	if(jj>150)
	{
		size_i=5*size_i;
		size_j=5*size_j;
		size_k=5*size_k;
	}
	//by theta
	int imin = ii-size_i;
	if(imin<0){imin=0;}
	int imax = ii+size_i;
	if(imax>=Nbin_i){imax=Nbin_i-1;}

	//by R
	int jmin = jj-size_j;
	if(jmin<0){jmin=0;}
	int jmax = jj+size_j;
	if(jmax>=Nbin_j){jmax=Nbin_j-1;}

	//by Z
	int kmin = kk-size_k;
	if(kmin<0){kmin=0;}
	int kmax = kk+size_k;
	if(kmax>=Nbin_k){kmax=Nbin_k-1;}

	double S_R=0;
	double S_RR=0;
	double S_sigVR2=0;
	double S_vR=0;
	double S_RsigVR2=0;
	double S_RvR=0;

	//double Vcru;

	double S_Z=0;
	double S_sigVRVZ=0;
	double S_ZZ=0;
	double S_ZsigVRVZ=0;

	double R,sigVR2,vR;
	double Z,sigVRVZ;
	int nR=0;
	int nZ=0;

//	double R0 = obj[ii][jj][kk].R;

//	by plane in TH_R
	for(int i = imin; i<=imax; ++i)
		for(int j = jmin; j<=jmax; ++j)
		{
			if(obj[i][j][kk].ok==1)
			{
				//fprintf(stderr,"R=%lf den = %lf sigVR2 = %lf\n",obj[i][j][kk].R,obj[i][j][kk].den,obj[i][j][kk].sigVR2);
				R = obj[i][j][kk].R;
				sigVR2 = obj[i][j][kk].sigVR2;
				vR = obj[i][j][kk].VR;

				S_R  += R;
				S_RR += R*R;
				S_sigVR2 += sigVR2;
				S_vR += vR;
				S_RsigVR2+= R*sigVR2;
				S_RvR+= R*vR;

				nR++;
				}
		}

	obj[ii][jj][kk].grad_sigVR2_R = (nR*S_RsigVR2-S_R*S_sigVR2)/(nR*S_RR-S_R*S_R);

	obj[ii][jj][kk].nR = nR;

	if(nR>1)
	{
	obj[ii][jj][kk].grad_VR_R = (nR*S_RvR-S_R*S_vR)/(nR*S_RR-S_R*S_R);
	}else
	{obj[ii][jj][kk].grad_VR_R = 0;}

	//	double Z0 = obj[ii][jj][kk].Z;

	for(int i = imin; i<=imax; ++i)
		for(int k = kmin; k<=kmax; ++k)
		{
			if( obj[i][jj][k].ok==1)
			{
			Z = obj[i][jj][k].Z;
			sigVRVZ = obj[i][jj][k].sigVRVZ;
			//sigVRVZ = obj[i][jj][k].VR*obj[i][jj][k].VZ;

			S_Z    += Z;
			S_ZZ   += Z*Z;
			S_sigVRVZ  += sigVRVZ;
			S_ZsigVRVZ += Z*sigVRVZ;

			nZ++;
			}
		}

	obj[ii][jj][kk].nZ = nZ;
	if(nZ>1)
	{
	obj[ii][jj][kk].grad_sigVRVZ_Z = (nZ*S_ZsigVRVZ-S_Z*S_sigVRVZ)/(nZ*S_ZZ-S_Z*S_Z);
	}else
	{obj[ii][jj][kk].grad_sigVRVZ_Z=0;}

//	obj[ii][jj][kk].grad_VR_R = 0;
//	obj[ii][jj][kk].grad_sigVRVZ_Z = 0;

	sigVR2 =  obj[ii][jj][kk].sigVR2;
	double  VR = obj[ii][jj][kk].VR;
	double sigVth2 = obj[ii][jj][kk].sigVth2;
	double dR = 0.1;


	double grad_sigVR2_R = obj[ii][jj][kk].grad_sigVR2_R ;
	double grad_VR_R = obj[ii][jj][kk].grad_VR_R;

	double grad_sigVRVZ_Z = obj[ii][jj][kk].grad_sigVRVZ_Z;
	sigVRVZ  = obj[ii][jj][kk].sigVRVZ;

		if(jj>150)
		{
		grad_sigVR2_R = ( fun_sigVR2(obj[ii][jj][kk].R+dR/2.0)-fun_sigVR2(obj[ii][jj][kk].R-dR/2.0)) / dR;
		sigVR2 =  fun_sigVR2(obj[ii][jj][kk].R);
		sigVth2 = fun_sigVth2(obj[ii][jj][kk].R);

		//grad_sigVR2_R = 0;
		//grad_VR_R = 0;
		grad_sigVRVZ_Z= 0;
		sigVRVZ  = 0;

		}



//	obj[ii][jj][kk].Vcru =  obj[ii][jj][kk].Vtheta*obj[ii][jj][kk].Vtheta + obj[ii][jj][kk].sigVth2  + (obj[ii][jj][kk].VR*obj[ii][jj][kk].VR+ obj[ii][jj][kk].sigVR2) * (obj[ii][jj][kk].R-h_R)/h_R
//							- 2.0*obj[ii][jj][kk].R*obj[ii][jj][kk].VR * obj[ii][jj][kk].grad_VR_R - obj[ii][jj][kk].R * obj[ii][jj][kk].grad_sigVR2_R
//							+ obj[ii][jj][kk].R*obj[ii][jj][kk].Z * obj[ii][jj][kk].sigVRVZ /(h_Z*fabs(obj[ii][jj][kk].Z) ) - obj[ii][jj][kk].R * obj[ii][jj][kk].grad_sigVRVZ_Z ;

	obj[ii][jj][kk].Vcru =  obj[ii][jj][kk].Vtheta*obj[ii][jj][kk].Vtheta + sigVth2  + (VR*VR+ sigVR2) * (obj[ii][jj][kk].R-h_R)/h_R
								- 2.0*obj[ii][jj][kk].R*VR * grad_VR_R - obj[ii][jj][kk].R * grad_sigVR2_R
								+ obj[ii][jj][kk].R*obj[ii][jj][kk].Z * sigVRVZ /(h_Z*fabs(obj[ii][jj][kk].Z) ) - obj[ii][jj][kk].R * grad_sigVRVZ_Z ;

	if(obj[ii][jj][kk].Vcru>0)
	{
		obj[ii][jj][kk].Vcru = sqrt(obj[ii][jj][kk].Vcru);

	}else
	{
		obj[ii][jj][kk].Vcru = 0;
	}
}

int calc_statistic(double x, int n, double * m_x, double * sqr_x, double * var_x )
{
	if(n==1)
	{
		*m_x = x;
		*sqr_x = x*x;
		*var_x = 0;
	}else
	{
	*m_x	+=	(x- *m_x)/ n;
	*sqr_x	+=  (x*x- *sqr_x)/ n;

	//дисперсия равна разности средней арифметической квадратов всех вариант статистической совокупности и квадрата средней самих этих вариант.
	*var_x = sqrt((*sqr_x-pow(*m_x,2))*n/(n-1));
	}
	return 0;
}


int main(int argc, char * argv[])
{
	//int Nz = 5;
	//int NR = 240;
	//int Nth = 120;

	double h_R = 3.0;
	double h_Z = 0.3;

	double zmin=-1.0;
	double zmax=1.0;
	//double dz=0.5;
	//int Nz = (int)((zmax-zmin)/dz);
	int Nz = 5;

	double Rmin=0.0;
	double Rmax=27;
	double dR=0.1;
	int NR = (int)((Rmax-Rmin)/dR);

	double Th_min=119;
	double Th_max=241;
	double dTh=1.0;
	int Nth = (int)((Th_max-Th_min)/dTh);

	fprintf(stderr,"Nz= %d\tNR= %d\tNTh= %d\n",Nz,NR,Nth);


	data *** obj = new  data ** [Nth];
	for(int i=0;i<Nth;++i)
	{
		obj[i] = new data * [NR];
		for(int j=0;j<NR;++j)
			{
				obj[i][j]= new data[Nz];
			}
	}

	fprintf(stderr,"Creat data array done\n");

	char * fname;
	fname = argv[1];

	FILE * f = fopen(fname,"r");

	if(!f)
	{
		fprintf(stderr,"Can not open file %s\n", fname);
	}

	double R,Th,Z;
	double avgR,avgTh,avgZ;
	double avgVR,avgVTh,avgVZ;
	double sigVRVZ,sigVR2,sigVTh2;

	fprintf(stderr,"Load data from %s file\n",fname);

	int n=0;
	int skip_n=0;
	int i,j,k;
	skipline(f);
	int n0=0;
	int n1=0;
	int n2=0;
	int n3=0;
	int n4=0;
	while(!feof(f))
	{
		if(fscanf(f,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&Z,&R,&Th,&avgZ,&avgR,&avgTh,&avgVR,&avgVTh,&avgVZ,&sigVRVZ,&sigVR2,&sigVTh2)==12)
		{
			Th = Th * rad2deg;
			i = get_index(Th,Th_max,Th_min,Nth);
			j = get_index(R,Rmax,Rmin,NR);
			k = get_index(Z,zmax,zmin,Nz);

			if( !(i<0 || j<0 || k<0) && sigVR2>0 && sigVTh2>0)
			{
				obj[i][j][k] = push_obj(avgZ,avgR,avgTh,avgVR,avgVTh,avgVZ,sigVRVZ,sigVR2,sigVTh2, h_R, h_Z);
				if(++n%100000==0)
				{
					fprintf(stderr,"Load %d objects for Th=%d R=%d Z=%d\n",n,i,j,k);
					fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",avgZ,avgR,avgTh,avgVR,avgVTh,avgVZ,sigVRVZ,sigVR2,sigVTh2);
				}
				switch (k)
				{
				case 0: n0++;
				break;
				case 1: n1++;
				break;
				case 2: n2++;
				break;
				case 3: n3++;
				break;
				case 4: n4++;
				break;
				}
			}else
			{
				++skip_n;
				fprintf(stderr,"skip_n n=%d  %lf %lf %lf i=%d j=%d k=%d sigVR2=%lf sigVTh2=%lf\n",n+1,Th,R,Z, i, j, k, sigVR2, sigVTh2 );
			}
		}else
		{
			skipline(f);
		}
	}

	fprintf(stderr,"Load %d  and skip %d objects for n0= %d n1= %d n2= %d n3= %d n4= %d\n",n,skip_n,n0,n1,n2,n3,n4);

	fprintf(stderr,"Start fill gradients\n");

	printf("i_th\tj_R\tk_Z\tR\tTheta\tZ\tVR\tVtheta\tVZ\tsigVR2\tsigVth2\tsigVRVZ\tgradVR_R\tgrad_sigVR2_R\tgrad_sigVRVZ_Z\tVcru_khrob\tdensity\tVcru_my\tgrad_sigVR2_VR2_R\tgrad_den_R\tgrad_den_VRVZ_Z\tNR\tNZ\n");

	//int iz = 2; // for z=0;
	//double Vcru;

	for(int iz=0;iz<Nz;++iz)
	for(int i=0;i<Nth;++i)
	{
		for(int j=0;j<NR;++j)
		{
			calc_gradient_khrob(obj,i,j,iz, 5, 5, 5, Nth, NR, Nz, h_R,h_Z );
			calc_gradient_my(obj,i,j,iz, 5,5,5, Nth, NR, Nz);

			if(obj[i][j][iz].ok==1)// && obj[i][j][iz].Vcru>0 && obj[i][j][iz].nR>5 && obj[i][j][iz].nZ>5)// && obj[i][j][iz].R<14
			{

			printf("%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",i,j,iz,
					obj[i][j][iz].R, obj[i][j][iz].theta, obj[i][j][iz].Z,
					obj[i][j][iz].VR, -obj[i][j][iz].Vtheta, obj[i][j][iz].VZ,
					obj[i][j][iz].sigVR2, obj[i][j][iz].sigVth2,obj[i][j][iz].sigVRVZ,
					obj[i][j][iz].grad_VR_R, obj[i][j][iz].grad_sigVR2_R,obj[i][j][iz].grad_sigVRVZ_Z, obj[i][j][iz].Vcru, obj[i][j][iz].den,obj[i][j][iz].Vcru_my,obj[i][j][iz].grad_sigVR2_VR2_R,obj[i][j][iz].grad_den_R,obj[i][j][iz].grad_den_sigVRVZ_Z, obj[i][j][iz].nR, obj[i][j][iz].nZ );
			}
		}


	}


	//avreging Vcru
	char filez [1024];
	char file [1024];
		for(int iz=0;iz<Nz;++iz)
		{
			sprintf(filez,"Vcru_z_%.1f.txt",(iz-2)*0.5);
			FILE * fz = fopen(filez,"w");
			if(!fz)
			{
				fprintf(stderr,"Can not open file %s to write results\n",filez);
				break;
			}
			fprintf(fz,"R\tsR\tTheta\tsTheta\tZ\tsZ\tVR\tsVR\tVtheta\tsVtheta\tVZ\tsVZ\tsigVR2\ts_sigVR2\tsigVth2\ts_sigVth2\tsigVRVZ\ts_sigVRVZ\tVcru\tsVcru\tNbin\n");




			for(int j=0;j<NR;j++)
			{
				data avg_z;
				data qrt_z;
				data var_z;
				int n_z=0;
				for(int i=30;i<Nth-30;i++)
				{
					if(obj[i][j][iz].Vcru>0)
					{
					n_z++;
					calc_statistic(obj[i][j][iz].R,n_z,&avg_z.R,&qrt_z.R,&var_z.R);
					calc_statistic(obj[i][j][iz].theta,n_z,&avg_z.theta,&qrt_z.theta,&var_z.theta);
					calc_statistic(obj[i][j][iz].Z,n_z,&avg_z.Z,&qrt_z.Z,&var_z.Z);

					calc_statistic(obj[i][j][iz].VR,n_z,&avg_z.VR,&qrt_z.VR,&var_z.VR);
					calc_statistic(obj[i][j][iz].Vtheta,n_z,&avg_z.Vtheta,&qrt_z.Vtheta,&var_z.Vtheta);
					calc_statistic(obj[i][j][iz].VZ,n_z,&avg_z.VZ,&qrt_z.VZ,&var_z.VZ);

					calc_statistic(obj[i][j][iz].sigVR2,n_z,&avg_z.sigVR2,&qrt_z.sigVR2,&var_z.sigVR2);
					calc_statistic(obj[i][j][iz].sigVth2,n_z,&avg_z.sigVth2,&qrt_z.sigVth2,&var_z.sigVth2);
					calc_statistic(obj[i][j][iz].sigVRVZ,n_z,&avg_z.sigVRVZ,&qrt_z.sigVRVZ,&var_z.sigVRVZ);

					calc_statistic(obj[i][j][iz].Vcru,n_z,&avg_z.Vcru,&qrt_z.Vcru,&var_z.Vcru);
					}
				}

			if(n_z>1)
			fprintf(fz,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
					avg_z.R,var_z.R,avg_z.theta,var_z.theta,avg_z.Z,var_z.Z,avg_z.VR,var_z.VR,avg_z.Vtheta,var_z.Vtheta,avg_z.VZ,var_z.VZ,
					avg_z.sigVR2,var_z.sigVR2,avg_z.sigVth2,var_z.sigVth2,avg_z.sigVRVZ,var_z.sigVRVZ,avg_z.Vcru,var_z.Vcru,n_z);
			}


			for(int ii=39;ii<80;ii+=20)
			{
				sprintf(file,"Vcru_th_%d_z_%.1f.txt",ii+1+120,(iz-2)*0.5);
				FILE * f = fopen(file,"w");
				if(!f)
				{
					fprintf(stderr,"Can not open file %s to write results\n",file);
					break;
				}
				fprintf(f,"R\tsR\tTheta\tsTheta\tZ\tsZ\tVR\tsVR\tVtheta\tsVtheta\tVZ\tsVZ\tsigVR2\ts_sigVR2\tsigVth2\ts_sigVth2\tsigVRVZ\ts_sigVRVZ\tVcru\tsVcru\tNbin\n");
				for(int j=0;j<NR;j++)
				{
					data avg;
					data qrt;
					data var;
					int n=0;
					for(int i=ii-10;i<ii+11;i++)
					{
						if(obj[i][j][iz].Vcru>0)
						{
							n++;

							calc_statistic(obj[i][j][iz].R,n,&avg.R,&qrt.R,&var.R);
							calc_statistic(obj[i][j][iz].theta,n,&avg.theta,&qrt.theta,&var.theta);
							calc_statistic(obj[i][j][iz].Z,n,&avg.Z,&qrt.Z,&var.Z);

							calc_statistic(obj[i][j][iz].VR,n,&avg.VR,&qrt.VR,&var.VR);
							calc_statistic(obj[i][j][iz].Vtheta,n,&avg.Vtheta,&qrt.Vtheta,&var.Vtheta);
							calc_statistic(obj[i][j][iz].VZ,n,&avg.VZ,&qrt.VZ,&var.VZ);

							calc_statistic(obj[i][j][iz].sigVR2,n,&avg.sigVR2,&qrt.sigVR2,&var.sigVR2);
							calc_statistic(obj[i][j][iz].sigVth2,n,&avg.sigVth2,&qrt.sigVth2,&var.sigVth2);
							calc_statistic(obj[i][j][iz].sigVRVZ,n,&avg.sigVRVZ,&qrt.sigVRVZ,&var.sigVRVZ);

							calc_statistic(obj[i][j][iz].Vcru,n,&avg.Vcru,&qrt.Vcru,&var.Vcru);

						}
					}

					if(n>1)
					fprintf(f,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
							avg.R,var.R,avg.theta,var.theta,avg.Z,var.Z,avg.VR,var.VR,avg.Vtheta,var.Vtheta,avg.VZ,var.VZ,
							avg.sigVR2,var.sigVR2,avg.sigVth2,var.sigVth2,avg.sigVRVZ,var.sigVRVZ,avg.Vcru,var.Vcru,n);
				}
			}

		}

	fprintf(stderr,"Stop fill gradients\n");


		for(int i=0;i<Nth;++i)
		{
			for(int j=0;j<NR;++j)
				{
					delete [] obj[i][j];
				}
			delete [] obj [i];
		}

		delete [] obj;

	return 0;
}

