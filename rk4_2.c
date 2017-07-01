/* Runge Kutta for a set of first order differential equations */

#include <stdio.h>
#include <math.h>

#define N 2                     
/* number of first order equations */
#define dist 0.05                
/* stepsize in t*/
#define MAX 20.0                
/* max for t */
#define PI 3.1415926535897  
FILE *output;                   
/* internal filename */

main()
{
double t, y[N];
int j;
 
void runge4(double x, double y[], double step); /* Runge-Kutta function */

double f(double x, double y[], int i);          /* function for derivatives */

output=fopen("c_program_output.dat", "w");                   /* external 
filename */

y[0]=PI/2;                                       
/* initial position */
y[1]=0.0;                                       /* initial velocity */
fprintf(output, "0\t%f\n", y[0]);
 
for (j=1; j*dist<=MAX ;j++)                     /* time loop */
{
   t=j*dist;
   runge4(t, y, dist);

   fprintf(output, "%f\t%f\n", t, y[0]);
}

fclose(output);
}

void runge4(double x, double y[], double step)
{
double h=step/2.0,                      /* the midpoint */
        t1[N], t2[N], t3[N],            /* temporary storage arrays */
        k1[N], k2[N], k3[N],k4[N];      /* for Runge-Kutta */
int i;
 
for (i=0;i<N;i++) t1[i]=y[i]+0.5*(k1[i]=step*f(x, y, i));
for (i=0;i<N;i++) t2[i]=y[i]+0.5*(k2[i]=step*f(x+h, t1, i));
for (i=0;i<N;i++) t3[i]=y[i]+    (k3[i]=step*f(x+h, t2, i));
for (i=0;i<N;i++) k4[i]=                step*f(x+step, t3, i);

for (i=0;i<N;i++) y[i]+=(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
}

double  f(double x, double y[], int i)
{
double M,D,I;
M=1.0; ///magnetism torque term
D=1.0;  ///Water drag term
I=1.0;  ///Moment of inertia term
if (i==0) return(y[1]);                 /* derivative of first equation */
if (i==1) return(-D/I*y[1]*y[1]-M/I*sin(y[0]));       
/* 
derivative of second equation */
}
