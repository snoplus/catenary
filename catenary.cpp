/*
 * Calculates the constants for a dual catenary with an air/D2O interface in the middle.
 *
 * For an explanation of the equations, see the Manipulator Log Book May 28/97 entry.
 *
 * Philip Harvey - 05/27/97
 *
 * Catenary equation:  y = y0 + h/w * cosh( w/h * (x-x0) )
 *
 * Where h = horizontal component of tension
 *		 w = weight per unit length of rope
 *
 * Write results to a file of the following format:
 *
 *	struct LookupRange {
 *		float	minX, maxX, incrX;
 *		float	minY, maxY, incrY;
 *		float	minLevel, maxLevel, incrLevel;
 *	};
 *
 * then the following structure repeated nx * ny times...
 *
 *	struct SourceInfo {
 *		float	x,y,level;
 *		float	length;			// length of umbilical (cm)
 *		float	driveTension;	// tension at umbilical drive (N)
 *		float	horizTension;	// horizontal tension at source (N, always +ve)
 *		float	vertTension;	// vertical tension at source (N, up is +ve)
 *		float	rubForce;		// force on rub ring (N)
 *	};
 *
 * Revisions:
 *
 * 02/25/98 - PH fixed interfaceSolve to handle x1>x3
 * 10/04/98 - PH fixed input of neck ring and feedthrough locations
 */

#ifndef LINUX
#define BORLANDC
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#ifdef BORLANDC
#include <conio.h>
#endif

#include "catenary.h"

//#define MONOTONIC		// monotonic approximations for catenary solution (test code only!)

/*
const double	kIncrX1 = 10;
const double	kIncrY1 = 10;
const double	kIncrY2 = 10;
const double	kMinX1 = -630;
const double	kMinY1 = -630;
const double	kMinY2 = -630;
const double	kMaxX1 = 0 + kIncrX1/2;
const double	kMaxY1 = 630 + kIncrY1/2;
const double	kMaxY2 = 630 + kIncrY2/2;	// water level is 1069
*/

#ifdef BORLANDC
typedef int Boolean;
const int true = 1;
const int false = 0;
#endif

#ifndef BORLANDC
static int kbhit() { return(0); }
static int getch() { return(0); }
#endif

#ifdef LINUX
typedef int Boolean;
#endif

//--------------------------------------------------------------------------------------------
/*
static double acosh(double x)
{
	if (x < 1) {
		printf("Invalid range for acosh()!\n");
		return(0);
	}
	return(log(x+sqrt(x*x-1)));
}*/
		
static char *gets_(char *buff)
{
	fgets(buff,100,stdin);
	char *pt = strchr(buff,'\0');
	while (pt>buff && (*(pt-1)=='\n' || *(pt-1)=='\r')) {
		*(--pt) = '\0';
	}
	return(buff);
}

//--------------------------------------------------------------------------------------------

enum ESolveCode {
	minTension,
	kLimitTensionSlope,
	kSetTension
};

class Catenary {
public:
			Catenary(double w=1, double h=1, double x1=0, double y1=0, double x2=1, double y2=1);
			
	void	setW(double w)					{ mW = w; mFlags = 0; }
	void	setH(double h)					{ mH = h; mFlags = 0; }
	void	setPoint1(double x1, double y1)	{ mX1 = x1; mY1 = y1; mFlags = 0; }
	void	setPoint2(double x2, double y2)	{ mX2 = x2; mY2 = y2; mFlags = 0; }
	void	setEndpoints(double x1,double y1,double x2,double y2) 
											{ setPoint1(x1,y1); setPoint2(x2,y2); }
	void	setNextSegment(Catenary *nxt=0)	{ mNext = nxt;		  }
	
	double	getX0();				// returns x0 of the catenary equation
	double	getY0();				// returns y0 of the catenary equation
	
	double	getH()		{ return mH;	}
	double	getX1()		{ return mX1;	}
	double	getY1()		{ return mY1;	}
	double	getX2()		{ return mX2;	}
	double	getY2()		{ return mY2;	}

	double	getY(double x);			// returns catenary y for a given x
	double	getX(double y, Boolean toLeft=false);
	double	getMinY();				// get the minimum point of the catenary
	double	getSlope(double x);		// reutrns slope at position x
	double	getLength();			// gets entire length of catenary system
	double	getLength(double x);	// calculate length from minimum to point x on catenary
	double	getTension(double x);	// calculate total tension at point x on catenary
	void	interfaceSolve();		// solve for a multiple catenary system
	
	static int	tensionSolve(Catenary *cat,double x1,double y1,double y2,double x3,double y3,
						double w1,double w2,ESolveCode code,double minTens=0,double minSlope=0,double minHoriz=0);
	
private:
	int			calcX0();		// solve the catenary equation for x0
	
	double		mX0, mY0;		// parameters in catenary equation
	double		mX1, mY1;		// first endpoint
	double		mX2, mY2;		// second endpoint
	double		mW;				// linear density
	double		mH;				// horizontal tension
	int			mFlags;			// flags set indicating what has been calculated
	Catenary  *	mNext;			// next segment of catenary

	enum {
		kCatFlagX0		= 0x01,
		kCatFlagY0		= 0x02,
		kCatFlagError	= 0x04,
		kCatFlagX0Y0	= kCatFlagX0 | kCatFlagY0
	};
	
};

Catenary::Catenary(double w, double h, double x1, double y1, double x2, double y2)
{
	mH = h;
	mW = w;
	mX1 = x1;
	mY1 = y1;
	mX2 = x2;
	mY2 = y2;
	mFlags = 0;
	mNext = 0;
}

// returns non-zero on error
int Catenary::calcX0()
{
	// handle the case where the rope is vertical
	if (mX1 == mX2) {
		mX0 = mX1;
		if (mY1 < mY2) mY0 = mY1;
		else		   mY0 = mY2;
		return(0);
	}
	
	double a = exp(-mW/mH*mX2) - exp(-mW/mH*mX1);
	double b = -2*mW/mH * (mY2-mY1);
	double c = exp( mW/mH*mX2) - exp( mW/mH*mX1);
	double rad = b*b - 4*a*c;

	if (rad < 0) return(-1);

	double temp = (-b+sqrt(rad)) / (2*a);
	if (temp <= 0) return(-1);

	mX0 = mH/mW * log(temp);	// solve for x0
	
	mFlags |= kCatFlagX0;
	
	return(0);
}

double Catenary::getX0()
{
	if (!(mFlags & kCatFlagX0)) calcX0();
	return(mX0);
}

double Catenary::getY0()
{
	if (!(mFlags & kCatFlagY0)) {
		if (!(mFlags & kCatFlagX0)) calcX0();
		mY0 = mY1 - mH/mW * cosh(mW/mH*(mX1-mX0));
		mFlags |= kCatFlagY0;
	}
	return(mY0);
}

double Catenary::getY(double x)
{
	return(getY0() + mH/mW * cosh(mW/mH*(x-getX0())));
}

// gets X for a given Y
// if toLeft is true, returns the lower x, otherwise returns the upper
double Catenary::getX(double y, Boolean toLeft)
{
	double tmp = mH/mW * acosh(mW/mH * (y - getY0()));
	if (toLeft) tmp = -tmp;
	return(getX0() + tmp);
}

double Catenary::getMinY()
{
	if (mX1 < mX2) {
		if (getX0() < mX1) return(mY1);
		if (getX0() > mX2) return(mY2);
	} else {
		if (getX0() < mX2) return(mY2);
		if (getX0() > mX1) return(mY1);
	}
	
	return(getY0() + mH/mW);
}

double Catenary::getSlope(double x)
{
	return(sinh(mW/mH*(x-getX0())));
}

double Catenary::getLength(double x)
{
	return(mH/mW * sinh(mW/mH * (x-getX0())));
}

double Catenary::getLength()
{
	double	length = fabs(getLength(mX2) - getLength(mX1));
	
	if (mNext) length += mNext->getLength();
	
	return(length);
}

double Catenary::getTension(double x)
{
	return(mH * cosh(mW/mH * (x-getX0())));
}

// find interface between this and the next catenary
// (mNext MUST be non-zero!)
void Catenary::interfaceSolve()
{
	if (!mNext) return;

	double			x2 = mX2;
	double			y2 = mY2;
	double			x01, x02;
	double			x2t, x2t1, x2t2;
	double			inter, slope;
	const double 	delta = 0.00001;	// change over which to evaluate slope
	const double 	tol = 0.000001;		// tolerance for convergence of interface x coordinate

	// iterate to find position where the rope hits the interface
	for (int iter=0; iter<100; ++iter) {

		// solve recursively for subsequent segments
		mNext->interfaceSolve();

		x2 -= delta / 2;	// shift back a bit

		setPoint2(x2, y2);   	   x01 = getX0();
		mNext->setPoint1(x2, y2);  x02 = mNext->getX0();

		// calculate new value for x2 at this location
		x2t1 = (mW*x01 - mNext->mW*x02) / (mW - mNext->mW);

		x2 += delta;		// shift forward

		setPoint2(x2, y2);  	   x01 = getX0();
		mNext->setPoint1(x2, y2);  x02 = mNext->getX0();

		// calculate new x2
		x2t2 = (mW*x01 - mNext->mW*x02) / (mW - mNext->mW);

		x2 -= delta / 2;	// restore original x2

		// calculate slope and intercept for newton's method
		slope = (x2t2-x2t1) / delta;
		inter = (x2t1+x2t2) / 2 - slope * x2;

		// make our best estimate of x2
		x2t = inter / (1 - slope);

		// make sure we don't go outside the valid range
		if (mX1 < mNext->mX2) {
			if (x2t >= mNext->mX2) x2t = mNext->mX2 - delta;
			if (x2t <= mX1) x2t = mX1 + delta;
		} else {
			if (x2t <= mNext->mX2) x2t = mNext->mX2 + delta;
			if (x2t >= mX1) x2t = mX1 - delta;
		}

		double diff = x2t - x2;		// calculate the difference

		x2 = x2t;					// save the new estimate of x2

		if (fabs(diff) < tol) {
			setPoint2(x2, y2);
			mNext->setPoint1(x2, y2);
			break;
		}
	}
}

//
// Find the catenaries for a specific tension situation given the endpoints,
// the location of the interface, and a strategy for choosing the tension.
//
// Inputs:	cat		- pointer to array of Catenaries for each segment of the rope
//			x1, y1	- starting point for the rope
//			y2		- y coordinate for the change in boyancy interface
//			x3, y3	- end point for the rope
//			w1		- linear weight BELOW y2
//			w2		- linear weight ABOVE y2
//			code	- one of:	minTension			- solve for minimum tension
//								kLimitTensionSlope	- limit tension and slope at x1,y1
//								kSetTension			- set rope tension to specified value, limiting slope
//			minTens	- tension used for kLimitTensionSlope and kSetTension
//			minSlope- minimum slope at x1,y1 (used for kLimitTensionSlope and kSetTension)
//			minHoriz- minimum horizontal tension component, used for kLimitTensionSlope only
//
// Returns:	number of catenaries or zero on convergence error
//
int Catenary::tensionSolve(Catenary *cat,double x1,double y1,double y2,double x3,double y3,
						   double w1,double w2,ESolveCode code,double minTens,double minSlope,double minHoriz)
{
	int				i, numCats;
	double		slopeSign;
	double		x2, x2b;
	double		h = 10;
	double		x01a[2];
	double		dh = 8;
	double		damping = 0.1;
	int				sign = 0;
#ifdef MONOTONIC
	double		oldH = h;
#else
	int				dampSign = 0;
#endif
	const double	minHHoriz = 1;	// minimum h for horizontal rope
	const double	delta_h = 0.0001;
	const double	h_convergence = 0.000005;
	const double	roundErr = 0.0001;	// minimum length of catenary

	if (x1 < x3) slopeSign = 1;
	else		 slopeSign = -1;

	cat[0].setPoint1(x1,y1);

	double dx = x3 - x1;
	double dy = y3 - y1;
	double minH = minHHoriz * (fabs(dx)/sqrt(dx*dx+dy*dy));

	// loop to find minimum tension
	for (;;) {

		cat[0].setH(h);

		// intially assume a single catenary
		cat[0].setPoint2(x3,y3);
		cat[0].setNextSegment();

		// take care of round-off errors
		if (y2+roundErr>y1 && y2-roundErr<y1) y2 = y1;

		// is more than half the height below water?
		if (y2 > (y1+y3)/2) {
			// get initial estimate for x2 from catenary in water
			cat[0].setW(w1);
			if (y2 >= y3) {
				numCats = 1;	// all under water
			} else {
				numCats = 2;
				x2 = cat[0].getX(y2, x1>x3);
			}
		} else {
			// initially assume entire catenary is in air
			cat[0].setW(w2);
			// is end of catenary submerged?
			if (y1 < y2) {
				numCats = 2;
				// get initial estimate for x2 from catenary in air
				x2 = cat[0].getX(y2, x1>x3);
			} else {
				// check to see if center of catenary is submerged
				double ymin = cat[0].getMinY(); // minimum point of catenary
				if (ymin >= y2) {	// just the one catenary
					numCats = 1;	// catenary all in air
				} else if (y1 == y2) {
					numCats = 2;	// endpoint is right on interface and catenary hangs down
					x2 = cat[0].getX(y2, x1>x3);
				} else {
					numCats = 3;	// crosses two air/water interfaces
					x2 = cat[0].getX(y2, x1>x3);
					x2b = cat[0].getX(y2, x1<=x3);	// initial guess for position of 2nd interface
				}
			}
		}

		// do calculation twice with slightly different h's to get estimate of dx0/dh
		for (int loop=0; loop<2; ++loop) {

			if (loop) h -= delta_h;
			else	  h += delta_h;
			
			// set the horizontal tensions in all catenaries
			for (i=0; i<numCats; ++i) cat[i].setH(h);
			
			if (numCats == 2) {
			
				cat[0].setPoint2(x2,y2);
				cat[1].setEndpoints(x2,y2,x3,y3);
				cat[0].setW(w1);		// source catenary is in water
				cat[1].setW(w2);
				cat[0].setNextSegment(&cat[1]);
				cat[1].setNextSegment();
				
				cat[0].interfaceSolve();			// solve the multiple catenary
				
				x2 = cat[0].getX2();	// get the interface location

			} else if (numCats == 3) {	// three catenaries
			
				cat[0].setPoint2(x2b,y2);
				cat[1].setEndpoints(x2b,y2,x2,y2);
				cat[2].setEndpoints(x2,y2,x3,y3);
				cat[0].setW(w2);
				cat[1].setW(w1);
				cat[2].setW(w2);
				cat[0].setNextSegment(&cat[1]);
				cat[1].setNextSegment(&cat[2]);
				cat[2].setNextSegment();

				cat[0].interfaceSolve();			// solve the multiple catenary
				
				x2b = cat[0].getX2();	// get the interface locations
				x2 = cat[1].getX2();
			}

			x01a[loop] = cat[0].getX0();	// save the value of x0 for this h
		}

		double	h2;
		double	tens = cat[0].getTension(x1);
		double	sl = cat[0].getSlope(x1);
		
		if (code == kLimitTensionSlope) {
			if (tens<minTens || sl*slopeSign<minSlope*slopeSign) h2 = h + dh;
			else h2 = h - dh;
			if (h2 < minHoriz) h2 = minHoriz;
			dh /= 2;
		} else if (code == kSetTension) {
			if (tens<minTens || sl*slopeSign<minSlope*slopeSign) h2 = h + dh;
			else h2 = h - dh;
			dh /= 2;
		} else {
			double dx01dh = (x01a[0]-x01a[1]) / delta_h;
			
			switch (numCats) {
				case 1:
				case 3:
					h2 = (x1 - x01a[1] * w2) / (1/tanh(w2/h*(x1-x01a[1])) - w2*dx01dh);
					break;
				case 2:
					h2 = (x1 - x01a[1] * w1) / (1/tanh(w1/h*(x1-x01a[1])) - w1*dx01dh);
					break;
			}
		}
		if (h2 < minH) h2 = minH;

		// has our tension converged?
		if (fabs((h-h2)/(h+h2)) < h_convergence) {
			for (i=0; i<numCats; ++i) cat[i].setH(h2);	// update with latest H
			break;// ALL DONE!!
		}
		if (code == minTension) {
#ifdef MONOTONIC
			int	newSign;
			if (h2 < h) newSign = -1;
			else newSign = 1;
			if (sign*newSign<0) {
				if (damping < 1e-10) {
					printf(" *%.1f*",h);
					return(0);	// give up
				}
				// increase damping factor and approach solution from the same side
				damping *= 0.1;
				h = oldH;	// restore last value of H to continue monotonic approach to solution
			} else {
				if (!sign) sign = newSign;
				oldH = h;	// save H so we can restore it if we start oscillating
				// calculate next estimate for horizontal tension
				h = h * (1-damping) + h2 * damping;
			}
#else // not MONOTONIC
			int	newSign;
			if (h2 < h) newSign = -1;
			else if (h2 > h) newSign = 1;
			else newSign = 0;
			if (sign*newSign<0) {
				if (!dampSign) dampSign = sign;
				if (sign == dampSign) {
					if (damping < 1e-10) {
						printf(" *%.1f*",h);
						return(0);	// give up
					}
					damping *= 0.1;	// reduce damping on oscillations
				}
			}
			sign = newSign;
			h = h * (1-damping) + h2 * damping;
#endif
		} else {
			h = h2;
		}
	}

	return(numCats);
}


//--------------------------------------------------------------------------------------------


int main()
{
	const double	xTop	= 0;	// always zero
	const double	minTens	= 2;		// minimum umbilical tension on source
	const double	kMinSlope = -10;

	double			slopeSign;
	double			minSlope = 0;
	Boolean			slopeSet = false;
	Boolean			hitNeck;
	Boolean			doLimit;

	double			x1, y1;				// source position
	double			y2;
	double			x3, y3;				// position of top end of rope
	double			w1, w2;				// linear densities in water and air

	Catenary		cat[10];				// catenaries for sections of umbilical (beginning at source)
	int					numCats;			// number of catenaries (1-4)
	char				filename[100];
	char				buff[256];
	char				*str1, *str2;
	char				*xStr, *yStr, *levelStr;

	double			minX1, maxX1, incrX1;
	double			minY1, maxY1, incrY1;
	double			minY2, maxY2, incrY2;
	double			yNeck, rNeck, yTop;

	SourceInfo	sourceInfo;
	LookupRange	lookup;
	int					verbose = 0;
	int					geom;
//
// This code assumes that  y3 > y1
//
	printf("\nOutput file name [catenary.dat]: ");
	gets_(filename);
	if (!filename[0]) strcpy(filename,"catenary.dat");
	FILE *fp = fopen(filename,"r+b");
	if (!fp) {

		str1 = "0.00064";
		printf("\nUmbilical linear density in water (N/cm) [%s]: ",str1);
		gets_(buff);
		if (!buff[0]) strcpy(buff,str1);
		w1 = atof(buff);

		str1 = "0.0148";
		printf("\nUmbilical linear density in air (N/cm) [%s]: ",str1);
		gets_(buff);
		if (!buff[0]) strcpy(buff,str1);
		w2 = atof(buff);

		printf("\nGeometry - 1) SNO or 2) Prototype [1]: ");
		gets_(buff);
		geom = atoi(buff);

		// get neck ring elevation
		switch (geom) {
			default:
				str1 = "580.685";
				break;
			case 2:
				str1 = "385.0";
				break;
		}
		printf("Neck ring elevation (cm) [%s]: ",str1);
		gets_(buff);
		if (!buff[0]) strcpy(buff,str1);
		yNeck = atof(buff);

		// get neck ring radius
		switch (geom) {
			default:
				str1 = "45.582";
				break;
			case 2:
				str1 = "3.85";
				break;
		}
		printf("Neck ring radius (cm) [%s]: ",str1);
		gets_(buff);
		if (!buff[0]) strcpy(buff,str1);
		rNeck = atof(buff);

		// get umbilical feedthrough elevation
		switch (geom) {
			default:
				str1 = "1450";
				break;
			case 2:
				str1 = "371.476";
				break;
		}
		printf("Umbilical feedthrough elevation (cm) [%s]: ",str1);
		gets_(buff);
		if (!buff[0]) strcpy(buff,str1);
		yTop = atof(buff);

		printf("\nSource position region - 1) Vessel or 2) Neck [1]: ");
		gets_(buff);

		switch (geom) {
			default:
				if (atoi(buff) == 2) {
					xStr = "0,80";
					yStr = "560,1450";
				} else {
					xStr = "0,620";
					yStr = "-620,620";
				}
				break;
			case 2:
				if (atoi(buff) == 2) {
					printf("I don't think so!\n");
					exit(1);
				} else {
					xStr = "0,320";
					yStr = "0,350";
				}
				break;
		}
		printf("Range of horizontal source positions (cm) [%s]: ",xStr);
		gets_(buff);
		if (!buff[0]) strcpy(buff,xStr);
		str1 = strtok(buff," ,");
		str2 = strtok(NULL," ,");
		if (!str1 || !str2) {
			printf("invalid range\n");
			exit(1);
		}
		minX1 = atof(str1);
		maxX1 = atof(str2);
		if (minX1 < 0) {
			printf("horizontal range truncated at zero for symmetry\n");
			minX1 = 0;
		}

		printf("Range of vertical source positions (cm) [%s]: ",yStr);
		gets_(buff);
		if (!buff[0]) strcpy(buff,yStr);
		str1 = strtok(buff," ,");
		str2 = strtok(NULL," ,");
		if (!str1 || !str2) {
			printf("invalid range\n");
			exit(1);
		}
		minY1 = atof(str1);
		maxY1 = atof(str2);

		printf("Grid spacing for source positions (cm) [10]: ");
		gets_(buff);
		if (!buff[0]) strcpy(buff,"10");
		incrX1 = incrY1 = atof(buff);
		if (incrX1<1 || incrX1>100) {
			printf("invalid spacing\n");
			exit(1);
		}

		printf("\nWater level region - 1) Vessel or 2) Neck [1]: ");
		gets_(buff);

		switch (geom) {
			default:
				if (atoi(buff) == 2) {
					levelStr = "580,1200";
				} else {
					levelStr = "-620,620";
				}
				break;
			case 2:
				if (atoi(buff) == 2) {
					printf("I don't think so!\n");
					exit(1);
				} else {
					levelStr = "-1000,-1000";
				}
				break;
		}
		printf("Range of water levels (cm) [%s]: ",levelStr);
		gets_(buff);
		if (!buff[0]) strcpy(buff,levelStr);
		str1 = strtok(buff," ,");
		str2 = strtok(NULL," ,");
		if (!str1 || !str2) {
			printf("invalid range\n");
			exit(1);
		}
		minY2 = atof(str1);
		maxY2 = atof(str2);

		printf("Spacing for water levels (cm) [10]: ");
		gets_(buff);
		if (!buff[0]) strcpy(buff,"10");
		incrY2 = atof(buff);
		if (incrY2<1 || incrY2>1000) {
			printf("invalid spacing\n");
			exit(1);
		}

		printf("\nCreating new output file\n");
		printf("File size will be %.2f MB\n",
				(maxX1-minX1+incrX1)/incrX1 * (maxY1-minY1+incrY1)/incrY1 *
				(maxY2-minY2+incrY2)/incrY2 * sizeof(SourceInfo) / 1048576.0);
		fp = fopen(filename,"wb");
		if (!fp) {
			printf("Error creating file %s\n",filename);
			exit(1);
		}

		printf("\nPress RTN to begin calculations");
		gets_(buff);
		printf("\n");

		lookup.minX = minX1; lookup.maxX = maxX1;  lookup.incrX = incrX1;
		lookup.minY = minY1; lookup.maxY = maxY1;  lookup.incrY = incrY1;
		lookup.minLevel = minY2; lookup.maxLevel = maxY2;  lookup.incrLevel = incrY2;
		lookup.w1 = w1; lookup.w2 = w2;
		lookup.yNeck = yNeck; lookup.rNeck = rNeck; lookup.yTop = yTop;
		if (fwrite(&lookup,sizeof(lookup),1,fp) != 1) {
			printf("error writing output file\n");
			exit(1);
		}
		fclose(fp);

		x1 = minX1;
		y1 = minY1;
		y2 = minY2;

	} else {

		printf("\nAdding to existing output file\n");
		if (fread(&lookup,sizeof(lookup),1,fp) != 1) {
			printf("error reading output file\n");
			exit(1);
		}
		minX1 = lookup.minX; maxX1 = lookup.maxX;  incrX1 = lookup.incrX;
		minY1 = lookup.minY; maxY1 = lookup.maxY;  incrY1 = lookup.incrY;
		minY2 = lookup.minLevel; maxY2 = lookup.maxLevel;  incrY2 = lookup.incrLevel;
		w1 = lookup.w1; w2 = lookup.w2;
		yNeck = lookup.yNeck; rNeck = lookup.rNeck; yTop = lookup.yTop;

		printf("X (%.0f,%.0f,%.0f)  Y (%.0f,%.0f,%.0f)  Level (%.0f,%.0f,%.0f)\n",
					minX1, maxX1, incrX1, minY1, maxY1, incrY1, minY2, maxY2, incrY2);

		fseek(fp,0L,SEEK_END);
		if (ftell(fp) < sizeof(SourceInfo)+sizeof(LookupRange)) {
			x1 = minX1;
			y1 = minY1;
			y2 = minY2;
		} else {
			fseek(fp,-(int)sizeof(SourceInfo),SEEK_END);
			if (fread(&sourceInfo,sizeof(SourceInfo),1,fp) != 1) {
				printf("error reading from file\n");
				exit(1);
			}
			x1 = sourceInfo.x;
			y1 = sourceInfo.y;
			y2 = sourceInfo.level;
		}
		fclose(fp);
		fp = fopen("minslope.dat","rb");
		if (fp && fread(&minSlope,sizeof(minSlope),1,fp)==1) {
			slopeSet = true;
		}
		if (fp) fclose(fp);

		printf("Starting from %.0f %.0f %.0f\n",x1,y1,y2);

		if ((y2+=incrY2) > maxY2+incrY2/2) {
			y2 = minY2;
			slopeSet = false;
			if ((y1+=incrY1) > maxY1+incrY1/2) {
				y1 = minY1;
				if ((x1+=incrX1) > maxX1+incrX1/2) {
					printf("All done calculations!\n");
					exit(1);
				}
			}
		}
		printf("\nPress RTN to continue, or ^C to abort. ");
		gets_(buff);
		printf("\n");
	}

	printf("Source   Water Drive Src   Horiz Vert  Seg Length  Rub    Interface\n");
	printf("X    Y    lvl  tens  tens  tens  Tens  num         force  locations\n");

	// loop through various levels for the water
	do {	// loop through x source position

		do {	// loop through y source position

			cat[0].setPoint1(x1,y1);
			doLimit = false;

			// open output file
			fp = fopen(filename,"ab");
			if (!fp) {
				printf("Error opening file!\n");
				exit(1);
			}

			do {	// loop through water levels

				if (kbhit()) {
					verbose ^= 1;
					if (getch() == 0x03) {
						printf("program halted\n");
						exit(1);
					}
				}

				if (verbose || y2 < minY2+incrY2/2) {
					printf("%4.0f %4.0f %4.0f",	x1, y1, y2);
				}

				double	length;
				double	driveTension;
				double	tens;
				double	sl;
				double	rubForce;
				double	horizTens;
				double	vertTens;

				if (fabs(x1-xTop) < 0.1) {
					// handle special case of vertical rope
					x3 = xTop; y3 = yTop;
					numCats = 1;
					length = y3 - y1;
					tens = driveTension = minTens;
					double yt = y2;
					if (yt > y3) yt = y3;
					if (yt > y1) {
						// add weight under water
						driveTension += w1 * (yt - y1);
					} else {
						yt = y1;
					}
					// add weight in air
					driveTension += w2 * (y3 - yt);
					rubForce = 0;
					horizTens = 0;
					vertTens = tens;
					if (x1>xTop) sl = -1e20;
					else	   sl = 1e20;

				} else {
					// decide if it is possible to hit neck ring
					if (y1>yNeck || (x1>-rNeck && x1<rNeck)) {
						x3 = xTop; y3 = yTop;
						hitNeck = false;
					} else {
						if (x1 < 0) x3 = -rNeck;
						else		x3 = rNeck;
						y3 = yNeck;
						hitNeck = true;		// may not actually hit neck -- we'll check later
					}

					if (x1 < x3) slopeSign = 1;
					else		 		 slopeSign = -1;

Try_Again:

					if (!doLimit) {
						numCats = Catenary::tensionSolve(cat,x1,y1,y2,x3,y3,w1,w2,minTension);

						tens = cat[0].getTension(x1);
						sl = cat[0].getSlope(x1);

						if (!slopeSet) {
							minSlope = sl;
							// limit minimum slope to prevent convergence on low-h solution
							if (minSlope * slopeSign < kMinSlope) {
								minSlope = slopeSign * kMinSlope;
							}
							slopeSet = true;
#ifndef NO_WRITES
							FILE *fp2 = fopen("minslope.dat","wb");
							if (!fp2) {
								printf("Error creating minslope file\n");
								exit(1);
							}
							fwrite(&minSlope,sizeof(minSlope),1,fp2);
							fclose(fp2);
#endif
						}
					} else {
						numCats = 0;
					}

					// if tension is too small, or the slope at the source hangs below the air-fill case,
					// solve again, limiting the tension and slope
					if (!numCats || tens<minTens || sl*slopeSign<minSlope*slopeSign) {
						doLimit = true;
						double minHoriz;
						if (numCats) {
							minHoriz = cat[0].getH();
						} else {
							minHoriz = 0;
						}
						numCats = Catenary::tensionSolve(cat,x1,y1,y2,x3,y3,w1,w2,
														 kLimitTensionSlope,minTens,minSlope,minHoriz);
						tens = cat[0].getTension(x1);
						sl = cat[0].getSlope(x1);
					}

					int		numVessel;
					double	neckTension;

					if (hitNeck) {		// determine if we really hit the neck
						numVessel = numCats;
						// must first solve for the neck catenaries
						neckTension = cat[numVessel-1].getTension(x3);
						numCats += Catenary::tensionSolve(cat+numCats,x3,y3,y2,xTop,yTop,w1,w2,kSetTension,neckTension,0);
						double theTens = cat[numVessel].getTension(x3);
						if (fabs(theTens-neckTension)/(theTens+neckTension) > 0.001) {
							printf("OOOOOPPPSSS!!\n");
							x3 = xTop; y3 = yTop;
							hitNeck = false;
							goto Try_Again;		// re-do calculations
						}

						// connect the two systems of catenaries
						cat[numVessel-1].setNextSegment(cat+numVessel);
						// if we hit the neck, then the slope will increase
						if (cat[numVessel-1].getSlope(x3)*slopeSign > cat[numVessel].getSlope(x3)*slopeSign) {
							// oops, we didn't hit the neck ring after all!
							x3 = xTop; y3 = yTop;
							hitNeck = false;
							// if this is air case, must re-calculate the minslope
							if (y2 < minY2+incrY2/2) {
								doLimit = false;
								slopeSet = false;
							}
							goto Try_Again;		// re-do calculations
						}
					}

					// calculate rubbing force on neck ring
					if (hitNeck) {
						double neckSlope1 = cat[numVessel-1].getSlope(x3);
						double neckSlope2 = cat[numVessel].getSlope(x3);
						double denom1 = sqrt(1+neckSlope1*neckSlope1);
						double denom2 = sqrt(1+neckSlope2*neckSlope2);
						double fx = neckTension * (1/denom1 - 1/denom2);
						double fy = neckTension * (neckSlope1/denom1 - neckSlope2/denom2);
						rubForce = sqrt(fx*fx + fy*fy);
					} else {
						rubForce = 0;
					}

					driveTension = cat[numCats-1].getTension(xTop);
					length = cat[0].getLength();
					horizTens = cat[0].getH();
					vertTens = sl * slopeSign * horizTens;
				}
				// print the results
				if (verbose || y2 < minY2+incrY2/2) {
					printf(" %5.2f %5.2f %5.2f %5.2f  %d %8.2f %5.2f",
							driveTension,tens,horizTens,vertTens,
							numCats,length,rubForce);
					for (int i=1; i<numCats; ++i) {
						printf("  %7.2f",cat[i].getX1());
					}
					printf("\n");
				}

				sourceInfo.x = x1;
				sourceInfo.y = y1;
				sourceInfo.level = y2;
				sourceInfo.length = length;
				sourceInfo.driveTension = driveTension;
				sourceInfo.horizTension = horizTens;
				sourceInfo.vertTension = vertTens;
				sourceInfo.rubForce = rubForce;
#ifndef NO_WRITES
				if (fwrite(&sourceInfo, sizeof(sourceInfo), 1, fp) != 1) {
					printf("Error writing to file!\n");
					exit(1);
				}
#endif

			} while ((y2+=incrY2) < maxY2+incrY2/2);
			y2 = minY2;

			fclose(fp);		// close output file to flush it

			minSlope = 0;
			slopeSet = false;

		} while ((y1+=incrY1) < maxY1+incrY1/2);
		y1 = minY1;

	} while ((x1+=incrX1) < maxX1+incrX1/2);
	return 0;
}

