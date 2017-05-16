#include <math.h>
#include <stdio.h>
#include "complex.h"

#define EPS 2.0e-6
#define MAXM 100
#define EPSS 6.e-8
#define MAXIT 100

void laguer(a,m,x,eps,polish)
fcomplex a[],*x;
int m,polish;
float eps;
{
	int j,iter;
	float err,dxold,cdx,abx;
	fcomplex sq,h,gp,gm,g2,g,b,d,dx,f,x1;
	dxold=Cabs(*x);
	for (iter=1;iter<=MAXIT;iter++) {
		b=a[m];
		err=Cabs(b);
		d=f=Complex(0.0,0.0);
		abx=Cabs(*x);
		for (j=m-1;j>=0;j--) {
			f=Cadd(Cmul(*x,f),d);
			d=Cadd(Cmul(*x,d),b);
			b=Cadd(Cmul(*x,b),a[j]);
			err=Cabs(b)+abx*err;
		}
		err *= EPSS;
		if (Cabs(b) <= err) return;
		g=Cdiv(d,b);
		g2=Cmul(g,g);
		h=Csub(g2,RCmul(2.0,Cdiv(f,b)));
		sq=Csqrt(RCmul((float) (m-1),Csub(RCmul((float) m,h),g2)));
		gp=Cadd(g,sq);
		gm=Csub(g,sq);
		if (Cabs(gp) < Cabs(gm))gp=gm;
		dx=Cdiv(Complex((float) m,0.0),gp);
		x1=Csub(*x,dx);
		if (x->r == x1.r && x->i == x1.i) return;
		*x=x1;
		cdx=Cabs(dx);
		dxold=cdx;
		if (!polish)
			if (cdx <= eps*Cabs(*x)) return;
	}
	
	//printf("Too many iterations in routine LAGUER!\n");

}

void zroots(a,m,roots,polish)
fcomplex a[],roots[];
int m,polish;
{
	int jj,j,i;
	fcomplex x,b,c,ad[MAXM];
	void laguer();

	for (j=0;j<=m;j++) ad[j]=a[j];
	for (j=m;j>=1;j--) {
		x=Complex(0.0,0.0);
		laguer(ad,j,&x,EPS,0);
		if (fabs(x.i) <= (2.0*EPS*fabs(x.r))) x.i=0.0;
		roots[j]=x;
		b=ad[j];
		for (jj=j-1;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=Cadd(Cmul(x,b),c);
		}
	}
	if (polish)
		for (j=1;j<=m;j++)
			laguer(a,m,&roots[j],EPS,1);
	for (j=2;j<=m;j++) {
		x=roots[j];
		for (i=j-1;i>=1;i--) {
			if (roots[i].r <= x.r) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
}

#undef EPS
#undef MAXM
#undef EPSS
#undef MAXIT