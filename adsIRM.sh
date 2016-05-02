#!/bin/bash
#PBS -N JAG227_AdsIRM
#PBS -o JAG227_Ads1.out
#PBS -e JAG227_Ads1.err
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -q shared

cd $PBS_O_WORKDIR
# change this for your environment
rm -R Out*
rm -R outputs/*
export RASPA_DIR="/ihome/cwilmer/jag227/RASPA-2.0"

	for ((i = 1; i < 79; i++))
	do
		CO2_COMP=$(awk NR==$i CO2c)
		CH4_COMP=$(awk NR==$i CH4c)
		N2_COMP=$(awk NR==$i N2c)
		C2H6_COMP=$(awk NR==$i C2H6c)



		if [ "$i" = 1 ]
		then
			grep Aver.*load Output/System_0/*00.data | awk NR==3 > co2IRMOF-1P1.txt
			grep Aver.*load Output/System_0/*00.data | awk NR==13 > ch4IRMOF-1P1.txt
      grep Aver.*load Output/System_0/*00.data | awk NR==23 > n2IRMOF-1P1.txt
      grep Aver.*load Output/System_0/*00.data | awk NR==33 > c2h6IRMOF-1P1.txt
			grep Aver.*load Output/System_0/*6.data | awk NR==3 > co2IRMOF-1P10.txt
      grep Aver.*load Output/System_0/*6.data | awk NR==13 > ch4IRMOF-1P10.txt
      grep Aver.*load Output/System_0/*6.data | awk NR==23 > n2IRMOF-1P10.txt
      grep Aver.*load Output/System_0/*6.data | awk NR==33 > c2h6IRMOF-1P10.txt

		else

			grep Aver.*load Output/System_0/*00.data | awk NR==3 >> co2IRMOF-1P1.txt
			grep Aver.*load Output/System_0/*00.data | awk NR==13 >> ch4IRMOF-1P1.txt
			grep Aver.*load Output/System_0/*00.data | awk NR==23 >> n2IRMOF-1P1.txt
			grep Aver.*load Output/System_0/*00.data | awk NR==33 >> c2h6IRMOF-1P1.txt
			grep Aver.*load Output/System_0/*6.data | awk NR==3 >> co2IRMOF-1P10.txt
      grep Aver.*load Output/System_0/*6.data | awk NR==13 >> ch4IRMOF-1P10.txt
      grep Aver.*load Output/System_0/*6.data | awk NR==23 >> n2IRMOF-1P10.txt
	    grep Aver.*load Output/System_0/*6.data | awk NR==33 >> c2h6IRMOF-1P10.txt
		fi
	done
