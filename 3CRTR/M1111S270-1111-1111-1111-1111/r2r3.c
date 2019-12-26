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
fp=fopen("rs_3.xls", "w");   //����ļ��� 

xtc=xdrfile_open ("dump_per200_2million.xtc","r");       //�����ļ��� 
read_xtc_natoms ("dump_per200_2million.xtc",&natoms);   //�����ļ��� 
x = calloc(natoms, sizeof (x[0]));


while (1)
  {
  int natom_1;
  double ax,ay,az,bx,by,bz,r,ratio,e_atom,e_vdw;
  read_return=read_xtc (xtc,natoms,&step,&time,box,x,&prec);
  if (read_return!=0) break;
  e_vdw = 0.0;
  for (natom=2025;natom<=3476;natom++)      //ת��2ԭ��ID 
    {
	  ax = x[natom-1][0], ay = x[natom-1][1], az = x[natom-1][2];

	  for (natom_1=3477;natom_1<=4928;natom_1++)     //ת��3ԭ��ID 
	  {
		bx = x[natom_1-1][0], by = x[natom_1-1][1], bz = x[natom_1-1][2];
		r = pow(ax-bx,2) + pow(ay-by,2) + pow(az-bz,2);
		r = pow(r,0.5);
		if (r < 10.2)     //�ضϰ뾶  ���� 
			ratio = 0.3407/r,
			e_atom = 4*2.968*0.001*(pow(ratio,12) - pow(ratio,6)),
     		e_vdw = e_vdw + e_atom;
	  }
  }
    

   printf("%d  %f   %lf\n",step/200,time,e_vdw);   //ÿ200��ȡһ��ƽ�� 

   fprintf(fp,"%d\t %f\t %lf\t\n",step/200,time,e_vdw);
  }
fclose(fp);
xdrfile_close (xtc);
}
