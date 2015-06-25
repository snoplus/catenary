struct SourceInfo {
	float	x,y,level;
	float	length;			// length of umbilical (cm)
	float	driveTension;	// tension at umbilical drive (N)
	float	horizTension;	// horizontal tension at source (N, always +ve)
	float	vertTension;	// vertical tension at source (N, up is +ve)
	float	rubForce;		// force on rub ring (N)
};

struct LookupRange {
	float	minX, maxX, incrX;
	float minY, maxY, incrY;
	float	minLevel, maxLevel, incrLevel;
	float	w1, w2;
	float	yNeck, rNeck, yTop;
};

