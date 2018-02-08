#ifdef GLVX_IMPLEMENTATION

static inline int fcmp(float a, float b) {
	return fabs(a - b) < 0.00001f;
}

/* Performs the quadratic formula on the specific coefficients a, b, and c,
* i.e. in a quadratic ax^2 + bx + c.
* Returns 0 if no roots exist, and 1 if roots exist.
* If roots do exist, they are written to result0 and result1 respectively.
*/
static inline int quadraticFormula(float a, float b, float c, float *result0, float *result1) {
	/* Special cases for lines and constants */
	if (fabs(a) < 0.00001f) {
		if (fabs(b) < 0.00001f) {
			return 0;
		}
		*result0 = *result1 = -c / b;
		return 1;
	}
	float inner = b * b - 4 * a *c;
	/* We call complex roots "nonexistent" */
	if (inner < 0) return 0;
	inner = (float)sqrt(inner) / (2 * a);
	b = -b / (2 * a);
	*result0 = b + inner;
	*result1 = b - inner;
	return 1;
}

/* Implementation of Fast Inverse Square Root
* Can always be replaced with a more optimal version
*/
static inline float fastInverseSquareRoot(float number) {
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y = number;
	i = *(long *)&y;
	i = 0x5F375A86 - (i >> 1);
	y = *(float *)&i;
	y = y * (threehalfs - (x2 * y * y));

	return y;
}

static inline size_t getCurveSamples(glvxCurve c) {
	glvxExtents extents;
	glvxGetBounds(c, extents);

	return (size_t)(sqrt(extents[2] * extents[2] + extents[3] * extents[3]) * 0.25 * sampleRatio);
}

static inline void fillCurve(glvxCurve c, float t0, float t1) {
	glBegin(GL_TRIANGLE_FAN);

	float samples = getCurveSamples(c) * (t1 - t0);
	float t = t0;
	float step = (t1 - t0) / (float)samples;
	for(size_t i = 0; i <= samples; ++i) {
		glVertexVector(curveVector(c, t));
		t += step;
	}
	glVertexVector(curveVector(c, t1));

	glEnd();
}

static inline float sign(float x1, float y1, float x2, float y2, float x3, float y3) {
	return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
}

static inline int isInside(struct PolygonPoint *potentialEar, float x, float y) {
	int b1, b2, b3;

	struct PolygonPoint *p = potentialEar->previous;
	struct PolygonPoint *n = potentialEar->next;

	b1 = sign(x, y, p->x, p->y, potentialEar->x, potentialEar->y) < 0.f;
	b2 = sign(x, y, potentialEar->x, potentialEar->y, n->x, n->y) < 0.f;
	b3 = sign(x, y, n->x, n->y, p->x, p->y) < 0.f;

	return (b1 == b2) && (b2 == b3);
}


#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif