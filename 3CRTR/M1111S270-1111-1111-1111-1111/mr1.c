#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"
#include <math.h>

int natoms,step,natom,read_return;
float time,prec;
matrix box;
rvec *x;
XDRFILE *xtc;

main ()
{
FILE *fp=NULL;
fp=fopen("rs_1.xls", "w");   //输出文件名 

xtc=xdrfile_open ("dump_per200_2million.xtc","r");       //读入文件名 
read_xtc_natoms ("dump_per200_2million.xtc",&natoms);   //读入文件名 
x = calloc(natoms, sizeof (x[0]));


while (1)
  {
  int natom_1;
  double ax,ay,az,bx,by,bz,r,ratio,e_atom,e_vdw;
  read_return=read_xtc (xtc,natoms,&step,&time,box,x,&prec);
  if (read_return!=0) break;
  e_vdw = 0.0;
  for (natom=1;natom<=572;natom++)      //电机ID 
    {
	  ax = x[natom-1][0], ay = x[natom-1][1], az = x[natom-1][2];

	  for (natom_1=573;natom_1<=2024;natom_1++)     //转子1原子ID 
	  {
		bx = x[natom_1-1][0], by = x[natom_1-1][1], bz = x[natom_1-1][2];
		r = pow(ax-bx,2) + pow(ay-by,2) + pow(az-bz,2);
		r = pow(r,0.5);
		if (r < 10.2)     //截断半径 埃米 
			ratio = 0.3407/r,
			e_atom = 4*2.968*0.001*(pow(ratio,12) - pow(ratio,6)),
     		e_vdw = e_vdw + e_atom;
	  }
  }
    

   printf("%d  %f   %lf\n",step/200,time,e_vdw);   //每200步取一次平均 

   fprintf(fp,"%d\t %f\t %lf\t\n",step/200,time,e_vdw);
  }
fclose(fp);
xdrfile_close (xtc);
}
