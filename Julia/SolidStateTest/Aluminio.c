/* Código C para el cálculo de la estructura de bandas del */
/* Aluminio utilizando el desarrollo de 4 ondas planas.  */  
/* (0,0,0) (0,0,2), (1,1,1) y (1,-1,1) módulo 2Pi/a */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N     5         /* Orden de la matriz más uno 4+1                */
#define M     1.32      /* Masa efectiva en el aluminio               */
#define V1 0.0295     /*Primer coeficiente de Fourier en Ry (1,1,1)    cambiar a 0.0 para modelo electrón libre*/
#define V2 0.0550    /*Segundo coeficiente de Fourier en Ry  (2,0,0) cambiar a 0.0 para modelo electrón libre    */
#define DK    0.005     /* Incremento o paso del vector recíproco            */
#define EPS   1e-6      /* Precisión para el cálculo tridiagonal                   */
#define NMAX  30        /* Número máximo de iteraciones en el cálculo tridiagonal*/
#define FILENAME  "bandas.out"  /* Archivo de datos de salida           */
#define SIGN(a) ((a)<0 ? -1 : 1)

double a[N][N];  /* Almacena la matriz simétrica del Hamiltoniano     */
double d[N];     /* Almacena los elementos de la diagonal al inicio y después devuelve los autovalores de la matriz             */
double e[N];     /* Almacena los elementos de la superdiagonal al
                    inicio y desaparece en la salida                  */

double mod2(double x,double y,double z) {
    return (x*x+y*y+z*z);/* Calcula el módulo de (k1,k2,k3)           */

}
/* Calcula la matriz hamiltoniana en el punto (k1,k2,k3)               */
void elementos(double k1,double k2,double k3) {
    double aux=1.0/M;
    a[1][1]=aux*mod2(k1,k2,k3);
    a[2][2]=aux*mod2(k1,k2,k3-2.0);
    a[3][3]=aux*mod2(k1-1.,k2-1.,k3-1.);
    a[4][4]=aux*mod2(k1-1.,k2+1.,k3-1.);
    a[1][2]=a[2][1]=a[3][4]=a[4][3]=V2;
    a[1][3]=a[3][1]=a[1][4]=a[4][1]=V1;
    a[2][3]=a[3][2]=a[2][4]=a[4][2]=V1;
}
/* Combierte la matriz hamiltoniana en tridiagonal    */

void tridiagonal(double d[N],double e[N],int n) {
   int i,j,iter=0,div;
   double c,s,t,ct2,p,se,ce,aux;
   for(j=1;j<=n;j++) {
      for(div=j;div<=n;div++) if(fabs(e[div])<EPS) break;
      while(div!=j) {
         if(div==j+1) {
            ct2=0.5*(d[j]-d[j+1])/e[j];
            t=SIGN(ct2)/(fabs(ct2)+sqrt(ct2*ct2+1.));
            c=1./sqrt(t*t+1.);
            s=t*c;
            d[j]+=t*e[j];
            d[j+1]-=t*e[j];
            e[j]=0.;
            j++;
         } else {
            if(iter++==NMAX) {
               printf("Demasiadas iteraciones en tridiagonal\n");
               exit(1);
            }
            ct2=0.5*(d[j]-d[j+1])/e[j];
            t=SIGN(ct2)/(fabs(ct2)+sqrt(ct2*ct2+1.));
            ct2=d[div]-d[j]+t*e[j];
            s=c=1.;
            p=0.;
            for(i=div-1;i>=j;i--) {
               se=s*e[i];
               ce=c*e[i];
               if(fabs(se)>=fabs(ct2)) {
                   c=ct2/se;
                   aux=sqrt(c*c+1.);
                   e[i+1]=se*aux;
                   s=1./aux;
                   c*=s;
               } else {
                   s=se/ct2;
                   aux=sqrt(s*s+1.);
                   e[i+1]=ct2*aux;
                   c=1./aux;
                   s*=c;
               }
               ct2=d[i+1]-p;
               aux=(d[i]-ct2)*s+2.*c*ce;
               p=s*aux;
               d[i+1]=ct2+p;
               ct2=c*aux-ce;
            }
            d[j]-=p;
            e[j]=ct2;
            e[div]=0.;
          }
          for(div=j;div<=n;div++) if(fabs(e[div])<EPS) break;
      }
   }
}

void householder(double a[N][N],double d[N],double e[N],int n) {
   int i,j,k;
   double s,s2,c,kappa,w[N],q[N];
   for(j=1;j<n-1;j++) {
      s2=0.;                             /* Calcula s al cuadrado     */
      for(i=j+1;i<=n;i++) s2+=a[i][j]*a[i][j];
      if(s2==0.) continue;/* Omite la transformacion ortogonal si s2=0*/
      s=-SIGN(a[j+1][j])*sqrt(s2);
      c=sqrt(2.*(s2-a[j+1][j]*s));
      w[j]=0.;                           /* Determina el vector w     */
      for(i=j+1;i<=n;i++) w[i]=a[i][j]/c;
      w[j+1]-=s/c;
      for(i=j;i<=n;i++) {                /* Determina el vector q=Aw  */
         q[i]=0.;
         for(k=j+1;k<=n;k++) q[i]+=a[i][k]*w[k];
      }
      kappa=0.;
      for(i=j;i<=n;i++) kappa+=w[i]*q[i];/* Calcula el producto wq    */
      for(i=j;i<=n;i++) q[i]-=kappa*w[i];/* Calcula q-=kappa w        */
      for(i=j;i<=n;i++)
         for(k=i;k<=n;k++) {             /* Determina A-=2(wq+qw)     */
            a[i][k]-=2.*(w[i]*q[k]+q[i]*w[k]);
            a[k][i]=a[i][k];
         }
   }
   for(i=1;i<n;i++) {
      d[i]=a[i][i];
      e[i]=a[i][i+1];
   }
   d[n]=a[n][n];
   e[n]=0.;
}

/* Calcula los autovalores en el punto (k1,k2,k3)                     */
void banda(double a[N][N],double d[N],double e[N],
           double k1,double k2,double k3)   {
   double tmp;
   int i,j;
   elementos(k1,k2,k3);
   householder(a,d,e,N-1);
   tridiagonal(d,e,N-1);
   for(i=1;i<N-1;i++)        /* Ordena los autovalores                */
      for(j=i+1;j<N;j++)
         if(d[j]<d[i]) {
            tmp=d[j];
            d[j]=d[i];
            d[i]=tmp;
         }
}
/* recorre los puntos y ejes de simetría     */
int main() {
   double k1,k2,k3,dk=DK,k=0.;
   int i;
   FILE *out;
   out=fopen(FILENAME,"w");
   for(k3=0.;k3<=1.;k3+=dk){                     /* De Gamma a X      */
      banda(a,d,e,0.,0.,k3);
      k+=1.;
      fprintf(out,"%8.6lf   ",k);
      for(i=1;i<N;i++) fprintf(out,"%8.6lf   ",d[i]);
      fprintf(out,"\n");
   }
   for(k1=0.;k1<=0.5;k1+=dk){                    /* De X a W          */
      banda(a,d,e,k1,0.,1.);
      k+=0.5;
      fprintf(out,"%8.6lf   ",k);
      for(i=1;i<N;i++) fprintf(out,"%8.6lf   ",d[i]);
      fprintf(out,"\n");
   }
   for(k2=0.,k3=1.;k2<=0.5;k2+=dk,k3-=dk){      /* De W a L           */
      banda(a,d,e,0.5,k2,k3);
      k+=0.5;
      fprintf(out,"%8.6lf   ",k);
      for(i=1;i<N;i++) fprintf(out,"%8.6lf   ",d[i]);
      fprintf(out,"\n");
   }
   for(k1=k2=k3=0.5;k1>=0.;k1-=dk,k2-=dk,k3-=dk){/* De L a Gamma      */
      banda(a,d,e,k1,k2,k3);
      k+=2.;
      fprintf(out,"%8.6lf   ",k);
      for(i=1;i<N;i++) fprintf(out,"%8.6lf   ",d[i]);
      fprintf(out,"\n");
   }
   fclose(out);
}
