/**
@file roots.c
@brief Polynomial root finder from Press et al. (1991)
*/
#include <math.h>
#include <stdio.h>
#include "complex.h"
#define MAXM 100

void laguer(dcomplex a[], int m, dcomplex *x, double polyeps1, double polyeps2, int maxpolyiter, int polish) {
	/*
	
	*/
	
	int j, iter;
	double err, dxold, cdx, abx;
	dcomplex sq,h,gp,gm,g2,g,b,d,dx,f,x1;
	dxold=Cabs(*x);
	for (iter=1;iter<=maxpolyiter;iter++) {
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
		err *= polyeps2;
		if (Cabs(b) <= err) return;
		g=Cdiv(d,b);
		g2=Cmul(g,g);
		h=Csub(g2,RCmul(2.0,Cdiv(f,b)));
		sq=Csqrt(RCmul((double) (m-1),Csub(RCmul((double) m,h),g2)));
		gp=Cadd(g,sq);
		gm=Csub(g,sq);
		if (Cabs(gp) < Cabs(gm))gp=gm;
		dx=Cdiv(Complex((double) m,0.0),gp);
		x1=Csub(*x,dx);
		if (x->r == x1.r && x->i == x1.i) return;
		*x=x1;
		cdx=Cabs(dx);
		dxold=cdx;
		if (!polish)
			if (cdx <= polyeps1*Cabs(*x)) return;
	}
	
	//TODO: This happens quite often!
	//printf("Too many iterations in routine LAGUER!\n");

}

void zroots(dcomplex a[], int m, dcomplex roots[], int polish, double polyeps1, double polyeps2, int maxpolyiter) {
	/*
	
	*/
	
	int jj,j,i;
	dcomplex x,b,c,ad[MAXM];

	for (j=0;j<=m;j++) ad[j]=a[j];
	for (j=m;j>=1;j--) {
		x=Complex(0.0,0.0);
		laguer(ad,j,&x,polyeps1,polyeps2,maxpolyiter,0);
		if (fabs(x.i) <= (2.0*polyeps1*fabs(x.r))) x.i=0.0;
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
			laguer(a,m,&roots[j],polyeps1,polyeps2,maxpolyiter,1);
	for (j=2;j<=m;j++) {
		x=roots[j];
		for (i=j-1;i>=1;i--) {
			if (roots[i].r <= x.r) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
}

#undef MAXM