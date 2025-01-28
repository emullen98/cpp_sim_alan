/*
Created: Alan Long Jun 20 2019
Edited: Alan Long Mar 02 2021
Edited: Ethan Mullen Jan 28 2025

This code simulates our mean-field model exactly. 
***IMPORTANT!!!***
Ethan modified this script on Jan 28 2025 to take in command line arguments.
I will describe those below.

It takes six inputs: 
(1) time_max (int) is the length of time you wish to simulate (in steps simulation timesteps). 
(2) Area (int) is the system size (in cells). 
(3) consv (double) is the conservation parameter c. 
(4) rate (double) is the strain rate, set to 0 for an adiabatic system. 
(5) w (double) is the spread of arrest stresses with failure stress normalized to 1. 
(6) weakening (double) is the weakening parameter epsilon.

Command line arguments that can be passed in any order:
-t or --time: Simulation timesteps. Default is 1000. Enter as a simple integer, i.e., 100000 instead of 1e6 or 100,000.
-s or --size: Simulation size (number of cells). Default is 1000. Enter as a simple integer.
-r or --rate: Simulation strain (driving) rate. Default value 0.0. Enter as a simple float.
-d or --disorder: Width of arrest stress distribution. Default value is 0.05. Enter as a simple float.
-w or --weakening: Weakening parameter epsilon. Default value is 0.0. Enter as a simple float.

It write two files as an output:
- stress_s={-s}_r={-r}_d={-d}_w={-w}.txt is the force on the system and consists of comma-spaced doubles.
- strain_s={-s}_r={-r}_d={-d}_w={-w}.txt is the strain on the system and consists of comma-spaced doubles.
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "cxxopts.hpp"

int area;
double *ranmarin(int ijkl, int N);
double *ranwbl(int ijkl, int N, double k, double lambda);

int main(int argc, char* argv[])
{
	cxxopts::Options options("MyProgram", "A simple program with default arguments.");

    options.add_options()
        ("t,time", "Simulation timesteps", cxxopts::value<int>()->default_value("1000"))
        ("s,size", "Simulation size (num of cells)", cxxopts::value<int>()->default_value("1000"))
        ("r,rate", "Strain rate", cxxopts::value<double>()->default_value("0.0"))
        ("d,disorder", "Disorder width", cxxopts::value<double>()->default_value("0.05"))
        ("w,weakening", "Weakening (epsilon)", cxxopts::value<double>()->default_value("0.0"));

    auto result = options.parse(argc, argv);

    int time_max = result["time"].as<int>();
    int area = result["size"].as<int>();
    double rate = result["rate"].as<double>();
    double w = result["disorder"].as<double>();
    double weakening = result["weakening"].as<double>();
    std::string stress_filename = "stress_s=" + std::to_string(area) + "_r=" + std::to_string(rate) + "_d=" + std::to_string(w) + "_w=" + std::to_string(weakening) + ".txt";
	std::string strain_filename = "strain_s=" + std::to_string(area) + "_r=" + std::to_string(rate) + "_d=" + std::to_string(w) + "_w=" + std::to_string(weakening) + ".txt";

	double consv;
	// Shear modulus
	double modulus;
	double a;
	double b;//this needs to be 1/Gamma(1+1/a), it's easier to just input it manually
	int numweak;	
	time_t t;
	srand((unsigned) time(&t));
	consv=1-1/sqrt(area);
	modulus=0.001;
	a=18.0;
	b=1.1985;
	b=1/tgamma(1+1/a);
	double *strs;
	strs=(double *) malloc(area*sizeof(double));//stress of each cell
	double *force;
	force=(double *)malloc(time_max*sizeof(double));//force on system
	double *strain;
	strain=(double *)malloc(time_max*sizeof(double));//strain on system
	int time;//current time step
	double *arr_strs;
	arr_strs=(double *)malloc(area*sizeof(double));//arrest stresses for each cell
	double *fail_strs;
	fail_strs=(double *)malloc(area*sizeof(double));//fail stresses for each cell
	int i;//iteration index
	double next_fail;//amount of stress to next failure
	int will_fail;//will a cell fail in the next time step
	double fail_amount;//force released by failures in a time step
	double *randm;
	randm=(double *)malloc(area*sizeof(double));//random number array
	FILE *fp;//file to write to

	randm=ranmarin(47,area);
		for (i=0;i<area;i++) {
			arr_strs[i]=0.1*randm[i]-0.05;//define arrest stresses randomly
		}
	free(randm);//need to free memory because it is allocated to heap
		for (i=0;i<area;i++) {	
			fail_strs[i]=1;//normalized
			strs[i]=(pow(((double)i/area),.4)*1.56585-.56585)*(fail_strs[i]-arr_strs[i])+arr_strs[i];//this distribution is an approximation of the steady state distribution. It will take some time to get to the steady state. This can be improved
		}
	time=0;
	strain[0]=0.0;
	for (i=0;i<area;i++){	
	force[time]+=strs[i];}//force is the sum of stresses	
	while (time<time_max) {
		next_fail=99.0;//this is just a large number
		for (i=0;i<area;i++) {
			if (next_fail>fail_strs[i]-strs[i])
				next_fail=fail_strs[i]-strs[i];//chose minimum stress to next failure
		}
		if (rate>0.0){
			int deltat=(int) (next_fail/(modulus*rate));
			if (deltat>0){
			for (i=0;i<deltat;i++){
				if (time==(time_max)){
					break;
				}
				force[time]=next_fail*((double) area)/((double) deltat)+force[time-1];
				strain[time]=rate;//+strain[time-1];
				time++;
			}
			}
		}
		else{
		strain[time]=next_fail*modulus+strain[time-1];
		}
		for (i=0;i<area;i++){
			strs[i]+=next_fail;//load to failure
		}
		will_fail=1;//says a cell will fail
		do {					//this loop is the failing mechanism
			
			printf("\b\b\b\b\b\b\b\b%d%% done",(time/(time_max/100)));//not necessary but tells you how far along the simulation is
						force[time]=0.0;
			for (i=0;i<area;i++){
				
				force[time]+=strs[i];}//force is the sum of stresses	
			++time;
			will_fail=0;//default to no fail
			fail_amount=0.0;
			
//			randm=ranmarin(time,area);//use this part to use a uniform distribution
//			for (i=0;i<area;i++){
//				randm[i]=1+w*randm[i]-w/2;
//			}
			
			randm=ranwbl(time,area,a,b);//use this part to use a weibull distribution
			strain[time]=0;//strain[time-1];//strain is cumulative
			for (i=0;i<area;i++){
				if (strs[i]>=fail_strs[i]){
					strain[time]+=1/(((double) area)*((double) area));//strain increases for each failure
					will_fail=1;//we now have a failure
					fail_amount+=(fail_strs[i]-arr_strs[i])*randm[i];//randomly distribute stress after failure based on arrest stress, add to stress lost
					strs[i]-=((fail_strs[i]-arr_strs[i])*randm[i])*(1+consv/(area-1));//subtract off the amount of stress lost
					if (fail_strs[i]==1)	fail_strs[i]=1-weakening*(1-arr_strs[i]);
				}

			}
			free(randm);
			for (i=0;i<area;i++){
				strs[i]+=fail_amount*consv/(area-1)+(1.0/consv-1.0)*rate;//distribute lost stress across system
			}		
		} while ((time<time_max) && will_fail);//only fail when a cell will fail and during the time period

		for (i=0;i<area;i++)
			fail_strs[i]=1;//reheal
	}
	
	const char* stress_filename_char = stress_filename.c_str();
	fp=fopen(stress_filename_char,"w");//open the stress file, change this filename to whatever you want
	for (time=0;time<time_max;time++){ fprintf(fp,"%f\n",force[time]);}//write the force array
	fclose(fp);
	
	const char* strain_filename_char = strain_filename.c_str();
	fp=fopen(strain_filename_char,"w");//open the stress file, change this filename to whatever you want
	fprintf(fp,"%f\n",strain[0]);
	for (time=1;time<time_max;time++){
		strain[time]+=strain[time-1];
	       	fprintf(fp,"%f\n",strain[time]);
	}//write the force array
	fclose(fp);

	// fp=fopen("C:\\c\\wblstrsweak0p0w4p0area1e5.txt","w");//open strain file, you can also change this
	// for (time=0;time<area;time++){ fprintf(fp,"%f,\n",strs[time]);}//write strain array
	// fclose(fp);

	free(strs);//need to free memory because heap allocation
	free(arr_strs);
	free(fail_strs);
	free(force);
	free(strain);
		}

/*
This is a random number generator. 
It has been edited. DON'T TOUCH THIS. It is the blackest of boxes. 
I don't know how it works, you won't know how it works. Seriously, just leave it be.
ijkl is a seed and N is the length of the resultant random array.
It outputs uni which is an array of random doubles betwen 0 and 1.
It is based on the code cited below.
Last edited: Alan Long 6/20/2019 and god help whoever edits it again.
*/

/*
Copyright (C) 1998 Matthew C. Kuntz, James P. Sethna, Karin A. Dahmen and John Carpenter.
This file is part of the Hysteresis program.
The Hysteresis program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.
See the file COPYING for details. 
[yeah I don't have this - Alan]
*/
double *ranmarin(int ijkl, int N)
{
double c,cd,cm,u[97];
int i97,j97,y ;
double *uni;
uni=(double *)malloc(N*sizeof(double));
double output;
int i,ii,j,jj,k,l,m ;
double s,t ;
int BIGPRIME;
BIGPRIME=899999963;
ijkl=ijkl%BIGPRIME;
	int ij=ijkl/30082;
	int kl=ijkl-30082*ij;

	i=((ij/177)%177)+2 ;
	j=(ij%177)+2 ;
	k=((kl/169)%178)+1 ;
	l=kl%169 ;
	for (ii=0;ii<97;ii++) {
		s=0.0 ;
		t=0.5 ;
		for (jj=0;jj<24;jj++) {
			m=(((i*j)%179)*k)%179 ;
			i=j;
			j=k;
			k=m;
			l=(53*l+1)%169;
			if (((l*m)%64)>=32) s+=t;
			t*=0.5;
		}
		u[ii]=s;
	}
	c=362436.0/16777216.0;
	cd=7654321.0/16777216.0;
	cm=16777213.0/16777216.0;
	i97=96;
	j97=32;
	for (y=0;y<N;y++){
		uni[y]=u[i97]-u[j97];
		if (uni[y]<0.0) uni[y]+=1.0;
		u[i97]=uni[y];
		if (--i97<0) i97=96;
		if (--j97<0) j97=96;
		c-=cd;
		if (c<0.0) c+=cm;
		uni[y]-=c;
		if (uni[y]<0.0) uni[y]+=1.0;
	}
return(uni);
}
	
/*
This is a makes a Weibull-ly distributed random number. 
It uses the uniform rng above. 
ijkl is a seed and N is the length of the resultant random array. 
k and lambda are the shape parameter and mean respectively. 
It outputs uni which is an array of random doubles.
It is based on the code cited below.
*/
double *ranwbl(int ijkl, int N, double k, double lambda)
{
	double *uni;
	double *wbl;
	wbl=(double *)malloc(N*sizeof(double));
	int i;	
	uni=ranmarin(ijkl,N);
	for(i=0;i<N;i++){
		if(uni[i]!=0.0){
		wbl[i]=lambda*pow(-1.0*log(uni[i]),1/k);
		}
		else{
			wbl[i]=10.0;
		}
	}
	free(uni);
	return(wbl);
}
