/*
 *  Parses umbilical catenery data to generate lookup tables for manip
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#include "catenary.h"

const double	kNeckRingHeight = 580.7075;

struct LookupInfo {
	float	x,y;			// source position for first calculated point
	float	dx,dy;		// distance between calculated points
	float	level, dlevel;	// initial water level, and water level increment
	short	nx,ny;		// number of calculated points in x and y
	short	nlevel;		// number of water levels (one per file)
	short	inNeck;		// non-zero if source is in neck
	char	nameFormat[16];	// file name formats for these tables
};

struct SourceStats {
	float	length;			// length of umbilical (cm)
	float	driveTension;	// tension at umbilical drive (N)
	float	horizTension;	// horizontal tension at source (N)
	float	vertTension;	// vertical tension at source (N, up is +ve)
	float	rubForce;		// force on rub ring (N)
};

static char *gets_(char *buff)
{
	fgets(buff,100,stdin);
	char *pt = strchr(buff,'\0');
	while (pt>buff && (*(pt-1)=='\n' || *(pt-1)=='\r')) {
		*(--pt) = '\0';
	}
	return(buff);
}

int main(void)
{
	const short	kMaxFiles	= 50;
	FILE				**out = new FILE*[kMaxFiles];
	float				*level = new float[kMaxFiles];
	FILE				*fp;
	LookupInfo	info;
	LookupRange	range;
	SourceInfo	sourceInfo;
	SourceStats	sourceStats;

	int					index = 1;
	char				buff[512];
	char				outname[100];
	double			lastX;

	printf("\n--- Catenary lookup table parser ---\n");
	printf("\nInput file name [catenary.dat]: ");
	gets_(buff);
	if (!buff[0]) strcpy(buff,"catenary.dat");

	FILE *fpin = fopen(buff,"rb");
	if (!fpin) {
		printf("Can't open input file '%s'\n",buff);
		exit(1);
	}
	if (fread(&range,sizeof(range),1,fpin) != 1) {
		printf("Error reading from %s\n",buff);
		exit(1);
	}

	printf("\nUmbilical linear density in water: %f\n",range.w1);
	printf("Umbilical linear density in air: %f\n",range.w2);
	printf("Acrylic vessel neck ring height: %f\n",range.yNeck);
	printf("Acrylic vessel neck ring radius: %f\n",range.rNeck);
	printf("Umbilical feedthrough elevation: %f\n",range.yTop);
	printf("X (%.0f,%.0f,%.0f)  Y (%.0f,%.0f,%.0f)  Level (%.0f,%.0f,%.0f)\n",
				range.minX, range.maxX, range.incrX,
				range.minY, range.maxY, range.incrY,
				range.minLevel, range.maxLevel, range.incrLevel);

	for (;;) {
		printf("\nUmbilical code letter: ");
		gets_(buff);
		if (strlen(buff)==1) {
			buff[0] = tolower(buff[0]);
			if (buff[0]>='a' && buff[0]<='z') break;
		}
		printf("Code must be a single character a through z\n");
	}
	sprintf(outname,"%c_index.dat",buff[0]);

	fp = fopen(outname,"rb");
	if (!fp) {
		// create the output file
		fp = fopen(outname,"wb");
		if (!fp) {
			printf("Can't create index file '%s'!\n",outname);
			exit(1);
		}
		printf("Created new index file %s\n",outname);
	} else {
		printf("Modifying existing index file %s...\n\n",outname);
		while (fread(&info, sizeof(info),1,fp) == 1) {
			printf("%d) ",index++);
			printf("x=%.0f+%.0f (%d) y=%.0f+%.0f (%d) level=%.0f+%.0f (%d) %s %s\n",
					info.x, info.dx, info.nx,
					info.y, info.dy, info.ny,
					info.level, info.dlevel, info.nlevel,
					info.nameFormat,
					info.inNeck ? "in neck" : "in vessel");
		}
	}
	fclose(fp);

	if (index > 1) {
		printf("\nIndex number for new entry [%d]? ", index);
		gets_(buff);
		int tmp = atoi(buff);
		if (tmp>0 && tmp<index) index = tmp;
	}
	sprintf(info.nameFormat,"%c%d_%%.0f.dat",outname[0],index);
	printf("Format for data file names %s\n",info.nameFormat);

	double	baseLevel = -1e20;

	for (;;) {

		fseek(fpin,sizeof(range),SEEK_SET);	// rewind input file

		// parse input file

		lastX = -1e20;

		int		doOpen = 1;
		int		numOpen = 0;

		while (fread(&sourceInfo,sizeof(sourceInfo),1,fpin) == 1) {
			if (sourceInfo.x > lastX) {
				printf("%.0f  \r",sourceInfo.x);
			}
			lastX = sourceInfo.x;

			if (sourceInfo.level < baseLevel-range.incrLevel/2) continue;

			int num = -1;

			if (doOpen && numOpen < kMaxFiles) {
				if (numOpen && sourceInfo.level < level[numOpen-1]+range.incrLevel/2) {
					out[numOpen] = 0;
				} else {
					sprintf(buff,info.nameFormat,sourceInfo.level);
					out[numOpen] = fopen(buff,"wb");
				}
				if (!out[numOpen]) {
					printf("Open %d outputs (%c%d_%.0f.dat to %c%d_%.0f.dat)\n",
							numOpen,outname[0],index,level[0],outname[0],index,level[numOpen-1]);
					doOpen = 0;
				} else {
					num = numOpen;
					if (!numOpen) {
						baseLevel = sourceInfo.level;
					}
					level[numOpen++] = sourceInfo.level;
				}
			}

			if (num < 0) {
				num = (int)((sourceInfo.level - baseLevel + range.incrLevel/2) / range.incrLevel);
			}
			if (num >= numOpen) continue;
			sourceStats.length = sourceInfo.length;
			sourceStats.driveTension = sourceInfo.driveTension;
			sourceStats.horizTension = sourceInfo.horizTension;
			sourceStats.vertTension = sourceInfo.vertTension;
			sourceStats.rubForce = sourceInfo.rubForce;
			if (fwrite(&sourceStats,sizeof(sourceStats),1,out[num]) != 1) {
				printf("error writing to output file %.0f\n",level[num]);
				exit(1);
			}
		}
		printf("Last entry: %.0f %.0f %.0f",sourceInfo.x,sourceInfo.y,sourceInfo.level);
		if (sourceInfo.x+range.incrX/2 < range.maxX ||
				sourceInfo.y+range.incrY/2 < range.maxY ||
				sourceInfo.level+range.incrLevel/2 < range.maxLevel) {
			printf("  WARNING: file is incomplete!!!\n");
		} else {
			printf("\n");
		}

		// close all open files
		for (int i=0; i<numOpen;++i) fclose(out[i]);

		if (!numOpen) break;

		// set the next base level
		baseLevel = level[numOpen-1] + range.incrLevel;

		// done if we have written all levels
		if (!numOpen || baseLevel > range.maxLevel+range.incrLevel/2) break;
	}
	fclose(fpin);

	delete out;
	delete level;
//
// finally, update the output information file
//
	fp = fopen(outname,"r+b");
	if (!fp) {
		printf("Error opening %s for update\n",outname);
		exit(1);
	}

	info.x = range.minX;
	info.y = range.minY;
	info.dx = range.incrX;
	info.dy = range.incrY;
	info.level = range.minLevel;
	info.dlevel = range.incrLevel;
	info.nx = (int)((range.maxX-range.minX+range.incrX/2)/range.incrX) + 1;
	info.ny = (int)((range.maxY-range.minY+range.incrY/2)/range.incrY) + 1;
	info.nlevel = (int)((range.maxLevel-range.minLevel+range.incrLevel/2)/range.incrLevel) + 1;
	info.inNeck = ((range.minLevel+range.maxLevel)/2 > kNeckRingHeight);

	// seek to the specified index
	fseek(fp,(index-1)*sizeof(info),SEEK_SET);
	// write the entry
	if (fwrite(&info,sizeof(info),1,fp) != 1) {
		printf("Error writing info file\n");
		exit(1);
	}
	fclose(fp);

	printf("\n%d) ",index);
	printf("x=%.0f+%.0f (%d) y=%.0f+%.0f (%d) level=%.0f+%.0f (%d) %s %s\n",
			info.x, info.dx, info.nx,
			info.y, info.dy, info.ny,
			info.level, info.dlevel, info.nlevel,
			info.nameFormat,
			info.inNeck ? "in neck" : "in vessel");

	printf("\ndone.\n");
	return 0;
}

