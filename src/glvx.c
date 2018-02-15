#include "GLVX.h"

#if defined(_WIN32) || defined(WIN32)
	#include <Windows.h>
#endif

#include <gl/GL.h>

#include <math.h>

struct PolygonPoint {
	struct PolygonPoint *previous;
	struct PolygonPoint *next;
	float x;
	float y;
};

static float sampleRatio = 1.f;

struct PolygonPoint *polygonPointList = NULL;
size_t polygonPointListSize = 0;

#define GLVX_IMPLEMENTATION
#include "glvxInternalVector.h"
#include "glvxInternalUtil.h"
#include "glvxInternalAdditionalStroke.h"

void glvxCleanup() {
	free(polygonPointList);
	polygonPointList = NULL;
	polygonPointListSize = 0;
}

void glvxSampleRatio(float newSampleRatio) {
	sampleRatio = newSampleRatio;
}

float glvxGetSampleRatio() {
	return sampleRatio;
}

float glvxCurveX(glvxCurve c, float t) {
	return ((c[0] * t + c[1]) * t + c[2]) * t + c[3];
}

float glvxCurveXBegin(glvxCurve c) {
	return c[3];
}

float glvxCurveXEnd(glvxCurve c) {
	return c[0] + c[1] + c[2] + c[3];
}

float glvxCurveY(glvxCurve c, float t) {
	return ((c[4] * t + c[5]) * t + c[6]) * t + c[7];
}

float glvxCurveYBegin(glvxCurve c) {
	return c[7];
}

float glvxCurveYEnd(glvxCurve c) {
	return c[4] + c[5] + c[6] + c[7];
}

void glvxGetExtents(glvxCurve c, glvxExtents extents) {
	float xmin = c[3], xmax = c[3];
	float ymin = c[7], ymax = c[7];

	/* Test curve endpoints */
	float test = glvxCurveXEnd(c), test1;
	if (test < xmin) xmin = test;
	if (test > xmax) xmax = test;
	test = glvxCurveYEnd(c);
	if (test < ymin) ymin = test;
	if (test > ymax) ymax = test;

	/* Test parametric X extrema (i.e. using the derivative) */
	if (quadraticFormula(3 * c[0], 2 * c[1], c[2], &test, &test1)) {
		if (test > 0 && test < 1) {
			test = glvxCurveX(c, test);
			if (test < xmin) xmin = test;
			if (test > xmax) xmax = test;
		}
		if (test1 > 0 && test1 < 1) {
			test1 = glvxCurveX(c, test1);
			if (test1 < xmin) xmin = test1;
			if (test1 > xmax) xmax = test1;
		}
	}

	/* Test parametric y extrema */
	if (quadraticFormula(3 * c[4], 2 * c[5], c[6], &test, &test1)) {
		if (test > 0 && test < 1) {
			test = glvxCurveY(c, test);
			if (test < ymin) ymin = test;
			if (test > ymax) ymax = test;
		}
		if (test1 > 0 && test1 < 1) {
			test1 = glvxCurveY(c, test1);
			if (test1 < ymin) ymin = test1;
			if (test1 > ymax) ymax = test1;
		}
	}

	/* Write extents */
	if (extents) {
		extents[0] = xmin;
		extents[1] = ymin;
		extents[2] = xmax;
		extents[3] = ymax;
	}
}

void glvxGetBounds(glvxCurve c, glvxExtents bounds) {
	if (!bounds) return;
	glvxGetExtents(c, bounds);
	bounds[2] -= bounds[0];
	bounds[3] -= bounds[1];
}

int glvxGetZeroCurvatureTime(glvxCurve c, float *t) {
	float commonA = -3 * c[0] * c[5] + 3 * c[1] * c[4];
	float commonB = 3 * c[2] * c[4] - 3 * c[0] * c[6];
	float commonC = c[2] * c[5] - c[1] * c[6];

	float quad0, quad1;

	if (quadraticFormula(commonA, commonB, commonC, &quad0, &quad1)) {
		if (quad0 >= 0 && quad0 <= 1) {
			*t = quad0;
			return 1;
		}
		if (quad1 >= 0 && quad1 <= 1) {
			*t = quad1;
			return 1;
		}
	}

	float xA = 3 * c[0];
	float xB = 2 * c[1];
	float xC = c[3];

	float yA = 3 * c[4];
	float yB = 2 * c[5];
	float yC = c[6];

	if (quadraticFormula(xA, xB, xC, &quad0, &quad1)) {
		if (fabs(yA * quad0 * quad0 + yB * quad0 + yC) < 0.0001) {
			if (quad0 >= 0 && quad0 <= 1) {
				*t = quad0;
				return 1;
			}
		}
		if (fabs(yA * quad1 * quad1 + yB * quad1 + yC) < 0.0001) {
			if (quad1 >= 0 && quad1 <= 1) {
				*t = quad1;
				return 1;
			}
		}
	}

	return 0;
}

int glvxGetZeroCurvaturePoint(glvxCurve c, float *x, float *y) {
	float t;
	if (glvxGetZeroCurvatureTime(c, &t)) {
		*x = glvxCurveX(c, t);
		*y = glvxCurveY(c, t);
		return 1;
	}
	return 0;
}

void glvxCalculateCurve(glvxCurve curve, float location[4], float ease[4]) {
	curve[3] = location[0];
	curve[7] = location[1];
	curve[2] = ease[0] * 3;
	curve[6] = ease[1] * 3;
	curve[1] = (location[2] - location[0]) * 3 - (ease[0] * 6) + (ease[2] * 3);
	curve[5] = (location[3] - location[1]) * 3 - (ease[1] * 6) + (ease[3] * 3);
	curve[0] = (location[2] - location[0]) - (ease[0] * 3) - curve[1];
	curve[4] = (location[3] - location[1]) - (ease[1] * 3) - curve[5];
}

void glvxStroke(glvxCurve c, float width) {
	size_t samples = getCurveSamples(c);
	if (samples < 2) return;
	float sampleSize = 1.f / samples;
	float sampleValue = sampleSize;

	/* We center the curve, and therefore use the half-width for projection */
	width /= 2;

	/* We are rendering the curve as a series of quads */
	glBegin(GL_QUAD_STRIP);

	Vector delta, orth, bi, proj;

	Vector previous = curveVector(c, 0);
	Vector current = curveVector(c, sampleValue);

	delta = subtract(current, previous);
	orth = orthogonalTo(delta);
	/* This instruction, which is used throughout this method, is simply approximately normalizing
	 * 'orth': it works well enough and is slightly faster than built in square root
	 */
	orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));

	glVertexVector(add(previous, multiply(orth, width)));
	glVertexVector(subtract(previous, multiply(orth, width)));
	
	while(samples > 1) {
		sampleValue += sampleSize;
		Vector next = curveVector(c, sampleValue);

		Vector nextDelta = subtract(current, next);
		
		/* The basic method is to get the bisector between the incoming direction
		 * and outgoing direction, and then project the curve width onto it.
		 */
		bi = bisector(delta, nextDelta);

		proj = divide(multiply(bi, width), dot(orth, bi));

		/* Add the two points (which lie on the bisector) */
		glVertexVector(add(current, proj));
		glVertexVector(subtract(current, proj));

		/* We calculate these vectors for the next loop */
		delta = multiply(nextDelta, -1);
		orth = orthogonalTo(delta);
		orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));
		current = next;
		previous = current;

		--samples;
	}

	orth = multiply(orth, width);

	/* Last point */
	current = curveVector(c, sampleValue);
	glVertexVector(add(current, orth));
	glVertexVector(subtract(current, orth));

	glEnd();
}

static inline void performEarcut(size_t count) {
	struct PolygonPoint *p = polygonPointList;

	glBegin(GL_TRIANGLES);

	size_t remainingPoints = count;

	size_t earError = 0; /* Used to prevent an infinite loop if bad geometry is specified */

	while (remainingPoints > 3) {
		struct PolygonPoint *test = p->next->next;
		int isEar = 1;

		/* Check if this is an ear */
		while (test->next != p) {
			if (isInside(p, test->x, test->y)) {
				if (test != p->next && test != p->previous)
					isEar = 0;
				break;
			}
			test = test->next;
		}

		if (isEar) {
			/* If this is an ear, cut it */
			glVertex2f(p->previous->x, p->previous->y);
			glVertex2f(p->x, p->y);
			glVertex2f(p->next->x, p->next->y);
			/* Relink the linked list */
			p->previous->next = p->next;
			p->next->previous = p->previous;
			p = p->next;
			--remainingPoints;
			earError = 0;
		}
		else {
			/* Check the next vertex */
			p = p->next;
			++earError;
			if (earError >= count) break;
		}
	}

	glVertex2f(p->previous->x, p->previous->y);
	glVertex2f(p->x, p->y);
	glVertex2f(p->next->x, p->next->y);

	glEnd();
}

void glvxEarcut(size_t count, float *polygon) {
#define list polygonPointList
#define listSize polygonPointListSize
	if (!list) {
		list = malloc(sizeof(struct PolygonPoint) * count);
		listSize = count;
	}
	if (listSize < count) {
		list = realloc(list, sizeof(struct PolygonPoint) * count);
		listSize = count;
	}

	/* Build the list */
	for (size_t i = 1; i < count - 1; ++i) {
		list[i].previous = (list + i - 1);
		list[i].next = (list + i + 1);
		list[i].x = polygon[i * 2];
		list[i].y = polygon[i * 2 + 1];
	}
	/* Special case for first point */
	list[0].x = polygon[0];
	list[0].y = polygon[1];
	list[0].previous = list + count - 1;
	list[0].next = list + 1;
	/* Special case for last point */
	list[count - 1].x = polygon[2 * (count - 1)];
	list[count - 1].y = polygon[2 * (count - 1) + 1];
	list[count - 1].previous = list + count - 2;
	list[count - 1].next = list;

	/* Perform earcut */
	performEarcut(count);

#undef list
#undef listSize
}

void glvxPaintMask(size_t count, float *curves, glvxExtents extents) {
#define list polygonPointList
#define listSize polygonPointListSize
	if (count < 1) return;

	/* Setup stencil */
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	/* We only need the first stencil bit for this algorithm */
	glStencilMask(1);
	size_t requiredPoints = count * 3;

	/* Ensure that earcut point list is good */
	if (!list) {
		list = malloc(sizeof(struct PolygonPoint) * requiredPoints);
		listSize = requiredPoints;
	}
	if (listSize < requiredPoints) {
		list = realloc(list, sizeof(struct PolygonPoint) * requiredPoints);
		listSize = requiredPoints;
	}

	/* Represents the total number of points to earcut */
	size_t totalCount = count + 1;
	/* Used to write earcut points in increasing order */
	size_t write = 0;

	/* The bounding box of all passed curves; initalize at the curves x/y begin */
	float xmin = curves[3], ymin = curves[7], xmax = curves[3], ymax = curves[7];

	float lastX = glvxCurveXBegin(curves);
	float lastY = glvxCurveYBegin(curves);

	/* Set stencil to toggle everytime a pixel is drawn */
	glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);

	for (size_t i = 0; i < count; ++i) {
		float *c = curves + i * 8;

		/* Get endpoints of curve */
		float x0 = glvxCurveXBegin(c);
		float y0 = glvxCurveYBegin(c);
		float x1 = glvxCurveXEnd(c);
		float y1 = glvxCurveYEnd(c);

		/* Determine which way this curve connects to previous curve */
		if (fcmp(x1, lastX) && fcmp(y1, lastY)) {
			list[write].x = x1;
			list[write++].y = y1;
			lastX = x0;
			lastY = y0;
		}
		else {
			/* If the curve is simply not connected, this will still work */
			list[write].x = x0;
			list[write++].y = y0;
			lastX = x1;
			lastY = y1;
		}

		/* Get curve extents */
		glvxExtents ce;
		glvxGetExtents(curves + i * 8, ce);

		/* If necessary, update bounding box */
		if (ce[0] < xmin) xmin = ce[0];
		if (ce[1] < ymin) ymin = ce[1];
		if (ce[2] > xmax) xmax = ce[2];
		if (ce[3] > ymax) ymax = ce[3];

		/* Get width and height for curve sample calculation */
		ce[2] -= ce[0];
		ce[3] -= ce[1];

		/* Simply fill along the whole curve:
		 * Any overlapping portions will be drawn twice and will cancel out
		 */
		fillCurve(c, 0, 1, CURVE_SAMPLES(ce[2], ce[3]));
	}
	list[write].x = lastX;
	list[write++].y = lastY;

	/* Build links for earcut */
	for (size_t i = 1; i < totalCount - 1; ++i) {
		list[i].previous = (list + i - 1);
		list[i].next = (list + i + 1);
	}
	list[0].previous = list + totalCount - 1;
	list[0].next = list + 1;
	list[totalCount - 1].previous = list + totalCount - 2;
	list[totalCount - 1].next = list;

	performEarcut(totalCount);

	/* Enable color mask and set stencil pass function */
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glStencilFunc(GL_EQUAL, 1, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

	/* Write extents to output */
	if (extents) {
		extents[0] = xmin;
		extents[1] = ymin;
		extents[2] = xmax;
		extents[3] = ymax;
	}
#undef list
#undef listSize
}

void glvxClearMask() {
	glDisable(GL_STENCIL_TEST);
}

void glvxFill(size_t count, float *curves) {
	glvxExtents extents;

	/* Paint stencil buffer */
	glvxPaintMask(count, curves, extents);

	/* Fill specified solid color across stencil buffer */
	glBegin(GL_QUADS);
		glVertex2f(extents[0], extents[1]);
		glVertex2f(extents[0], extents[3]);
		glVertex2f(extents[2], extents[3]);
		glVertex2f(extents[2], extents[1]);
	glEnd();

	/* Setup stencil to zero out */
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glStencilFunc(GL_ALWAYS, 0, 0xff);
	glStencilOp(GL_ZERO, GL_ZERO, GL_ZERO);

	/* Zero out the stencil over the bounding box */
	glBegin(GL_QUADS);
		glVertex2f(extents[0], extents[1]);
		glVertex2f(extents[0], extents[3]);
		glVertex2f(extents[2], extents[3]);
		glVertex2f(extents[2], extents[1]);
	glEnd();

	/* Re-enable color mask */
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	glvxClearMask();
}