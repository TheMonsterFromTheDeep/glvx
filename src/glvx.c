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
float *fillPointList = NULL;
size_t fillPointListSize = 0;

#define GLVX_IMPLEMENTATION
#include "glvxInternalVector.h"
#include "glvxInternalUtil.h"
#include "glvxInternalAdditionalStroke.h"

void glvxCleanup() {
	free(polygonPointList);
	free(fillPointList);
	polygonPointList = NULL;
	fillPointList = NULL;
	polygonPointListSize = 0;
	fillPointListSize = 0;
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

float glvxCurveY(glvxCurve c, float t) {
	return ((c[4] * t + c[5]) * t + c[6]) * t + c[7];
}

void glvxGetExtents(glvxCurve c, glvxExtents extents) {
	float xmin = c[3], xmax = c[3];
	float ymin = c[7], ymax = c[7];

	/* Test curve endpoints */
	float test = glvxCurveX(c, 1), test1;
	if (test < xmin) xmin = test;
	if (test > xmax) xmax = test;
	test = glvxCurveY(c, 1);
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

	glBegin(GL_QUAD_STRIP);

	Vector delta, orth, bi, proj;

	Vector previous = curveVector(c, 0);
	Vector current = curveVector(c, sampleValue);

	delta = subtract(current, previous);
	orth = orthogonalTo(delta);
	orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));

	glVertexVector(add(previous, multiply(orth, width)));
	glVertexVector(subtract(previous, multiply(orth, width)));
	
	while(samples > 1) {
		sampleValue += sampleSize;
		Vector next = curveVector(c, sampleValue);

		Vector nextDelta = subtract(current, next);
		
		bi = bisector(delta, nextDelta);

		proj = divide(multiply(bi, width), dot(orth, bi));

		glVertexVector(add(current, proj));
		glVertexVector(subtract(current, proj));

		delta = multiply(nextDelta, -1);
		orth = orthogonalTo(delta);
		orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));
		current = next;
		previous = current;

		--samples;
	}

	orth = multiply(orth, width);

	current = curveVector(c, sampleValue);
	glVertexVector(add(current, orth));
	glVertexVector(subtract(current, orth));

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
		list = realloc(list, listSize);
		listSize = count;
	}

	for (size_t i = 1; i < count - 1; ++i) {
		list[i].previous = (list + i - 1);
		list[i].next = (list + i + 1);
		list[i].x = polygon[i * 2];
		list[i].y = polygon[i * 2 + 1];
	}
	list[0].x = polygon[0];
	list[0].y = polygon[1];
	list[0].previous = list + count - 1;
	list[0].next = list + 1;
	list[count - 1].x = polygon[2 * (count - 1)];
	list[count - 1].y = polygon[2 * (count - 1) + 1];
	list[count - 1].previous = list + count - 2;
	list[count - 1].next = list;

	struct PolygonPoint *p = list;

	glBegin(GL_TRIANGLES);

	size_t remainingPoints = count;
	size_t earError = 0;
	while (remainingPoints > 3) {
		struct PolygonPoint *test = p->next->next;
		int isEar = 1;
		while (test->next != p) {
			if (isInside(p, test->x, test->y)) {
				if (test != p->next && test != p->previous)
					isEar = 0;
				break;
			}
			test = test->next;
		}

		if (isEar) {
			glVertex2f(p->previous->x, p->previous->y);
			glVertex2f(p->x, p->y);
			glVertex2f(p->next->x, p->next->y);
			p->previous->next = p->next;
			p->next->previous = p->previous;
			p = p->next;
			--remainingPoints;
			earError = 0;
		}
		else {
			p = p->next;
			++earError;
			if (earError >= count) break;
		}
	}

	glVertex2f(p->previous->x, p->previous->y);
	glVertex2f(p->x, p->y);
	glVertex2f(p->next->x, p->next->y);

	glEnd();

#undef list
#undef listSize
}

void glvxPaintMask(size_t count, float *curves, glvxExtents extents) {
#define list fillPointList
#define listSize fillPointListSize
	if (count < 1) return;

	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glStencilMask(1);
	glStencilFunc(GL_ALWAYS, 0, 0xff);
	glStencilOp(GL_ZERO, GL_ZERO, GL_ZERO);
	size_t requiredPoints = count * 3 * 2;
	if (!list) {
		list = malloc(sizeof(float) * requiredPoints);
		listSize = requiredPoints;
	}
	if (listSize < requiredPoints) {
		list = realloc(list, sizeof(float) * requiredPoints);
		listSize = requiredPoints;
	}
	size_t totalCount = count + 1;
	size_t write = 0;

	float xmin = curves[3], ymin = curves[7], xmax = curves[3], ymax = curves[7];

	float lastX = glvxCurveX(curves, 0);
	float lastY = glvxCurveY(curves, 0);

	for (size_t i = 0; i < count; ++i) {
		glvxExtents ce;
		glvxGetExtents(curves + i * 8, ce);
		if (ce[0] < xmin) xmin = ce[0];
		if (ce[1] < ymin) ymin = ce[1];
		if (ce[2] > xmax) xmax = ce[2];
		if (ce[3] > ymax) ymax = ce[3];
	}

	glBegin(GL_QUADS);
		glVertex2f(xmin, ymin);
		glVertex2f(xmin, ymax);
		glVertex2f(xmax, ymax);
		glVertex2f(xmax, ymin);
	glEnd();

	glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);

	for (size_t i = 0; i < count; ++i) {
		float *c = curves + i * 8;

		float x0 = glvxCurveX(c, 0);
		float y0 = glvxCurveY(c, 0);
		float x1 = glvxCurveX(c, 1);
		float y1 = glvxCurveY(c, 1);

		if (fcmp(x1, lastX) && fcmp(y1, lastY)) {
			list[write++] = x1;
			list[write++] = y1;
			lastX = x0;
			lastY = y0;
		}
		else {
			list[write++] = x0;
			list[write++] = y0;
			lastX = x1;
			lastY = y1;
		}

		float time;
		if (glvxGetZeroCurvatureTime(curves + i * 8, &time)) {
			list[write++] = glvxCurveX(c, time);
			list[write++] = glvxCurveY(c, time);
			++totalCount;

			fillCurve(c, 0, time);
			fillCurve(c, time, 1);
		}
		else {
			fillCurve(c, 0, 1);
		}
	}
	list[write++] = lastX;
	list[write++] = lastY;
	glvxEarcut(totalCount, list);

	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glStencilFunc(GL_EQUAL, 1, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

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

	glvxPaintMask(count, curves, extents);

	glBegin(GL_QUADS);
		glVertex2f(extents[0], extents[1]);
		glVertex2f(extents[0], extents[3]);
		glVertex2f(extents[2], extents[3]);
		glVertex2f(extents[2], extents[1]);
	glEnd();

	glvxClearMask();
}