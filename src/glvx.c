#include "GLVX.h"

#if defined(_WIN32) || defined(WIN32)
	#include <Windows.h>
#endif

#include <gl/GL.h>

#include <math.h>

#define GLVX_IMPLEMENTATION
#include "glvxInternalParameters.h"
#include "glvxInternalCurveProperties.h"
#include "glvxInternalVector.h"
#include "glvxInternalUtil.h"

float *curveEvenOddList = NULL;
size_t curveEvenOddListSize = 0;

void glvxCleanup() {
	free(curveEvenOddList);
	curveEvenOddList = NULL;
	curveEvenOddListSize = 0;
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

void glvxStroke(glvxCurve c) {
	size_t samples = getCurveSamples(c);
	if (samples < 2) return;
	float sampleSize = 1.f / samples;
	float sampleValue = sampleSize;

	float width = beginWidth;
	float widthStep = (endWidth - beginWidth) / samples;

	float colorStep[4];
	float color[4];

	color[0] = beginColor[0];
	color[1] = beginColor[1];
	color[2] = beginColor[2];
	color[3] = beginColor[3];

	colorStep[0] = (endColor[0] - beginColor[0]) / samples;
	colorStep[1] = (endColor[1] - beginColor[1]) / samples;
	colorStep[2] = (endColor[2] - beginColor[2]) / samples;
	colorStep[3] = (endColor[3] - beginColor[3]) / samples;

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

	glColor4f(color[0] * leftColor[0], color[1] * leftColor[1], color[2] * leftColor[2], color[3] * leftColor[3]);
	glVertexVector(add(previous, multiply(orth, width * (1 + strokeOffset))));
	glColor4f(color[0] * rightColor[0], color[1] * rightColor[1], color[2] * rightColor[2], color[3] * rightColor[3]);
	glVertexVector(subtract(previous, multiply(orth, width * (1 - strokeOffset))));
	
	while(samples > 1) {
		sampleValue += sampleSize;
		width += widthStep;
		color[0] += colorStep[0];
		color[1] += colorStep[1];
		color[2] += colorStep[2];
		color[3] += colorStep[3];
		Vector next = curveVector(c, sampleValue);

		Vector nextDelta = subtract(current, next);
		
		/* The basic method is to get the bisector between the incoming direction
		 * and outgoing direction, and then project the curve width onto it.
		 */
		bi = bisector(delta, nextDelta);

		proj = divide(multiply(bi, width), dot(orth, bi));
		if (squaredLength(proj) / (width * width) > (miterLimit * miterLimit)) {
			proj = multiply(proj, fastInverseSquareRoot(squaredLength(proj)));
			proj = multiply(proj, miterLimit * width);
		}

		/* Add the two points (which lie on the bisector) */
		glColor4f(color[0] * leftColor[0], color[1] * leftColor[1], color[2] * leftColor[2], color[3] * leftColor[3]);
		glVertexVector(add(current, multiply(proj, 1 + strokeOffset)));
		glColor4f(color[0] * rightColor[0], color[1] * rightColor[1], color[2] * rightColor[2], color[3] * rightColor[3]);
		glVertexVector(subtract(current, multiply(proj, 1 - strokeOffset)));

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
	glColor4f(color[0] * leftColor[0], color[1] * leftColor[1], color[2] * leftColor[2], color[3] * leftColor[3]);
	glVertexVector(add(current, multiply(orth, 1 + strokeOffset)));
	glColor4f(color[0] * rightColor[0], color[1] * rightColor[1], color[2] * rightColor[2], color[3] * rightColor[3]);
	glVertexVector(subtract(current, multiply(orth, 1 - strokeOffset)));

	glEnd();
}

inline void glvxEvenOdd(size_t count, float *points) {
	/* Setup stencil */
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	/* We only need the first stencil bit for this algorithm */
	glStencilMask(1);

	/* Set stencil to toggle everytime a pixel is drawn */
	glStencilFunc(GL_ALWAYS, 0, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);

	glBegin(GL_TRIANGLE_FAN);
	for(size_t i = 0; i < count; ++i) {
		glVertex2f(points[i * 2], points[i * 2 + 1]);
	}
	glEnd();
}

void glvxStencil(size_t count, float *curves, glvxExtents extents) {
#define list curveEvenOddList
#define listSize curveEvenOddListSize
	if (count < 1) return;

	/* Setup stencil */
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	/* We only need the first stencil bit for this algorithm */
	glStencilMask(1);
	size_t requiredListItems = count * 3 * 2;

	/* Ensure that earcut point list is good */
	if (!list) {
		list = malloc(sizeof(float) * requiredListItems);
		listSize = requiredListItems;
	}
	if (listSize < requiredListItems) {
		list = realloc(list, sizeof(float) * requiredListItems);
		listSize = requiredListItems;
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
	glStencilFunc(GL_ALWAYS, 0, 0xff);
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
			list[write++] = x1;
			list[write++] = y1;
			lastX = x0;
			lastY = y0;
		}
		else {
			/* If the curve is simply not connected, this will still work */
			list[write++] = x0;
			list[write++] = y0;
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
		fillCurveSpecific(c, 0, 1, CURVE_SAMPLES(ce[2], ce[3]));
	}
	list[write++] = lastX;
	list[write++] = lastY;

	glvxEvenOdd(totalCount, list);

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

void glvxDisableStencil() {
	glDisable(GL_STENCIL_TEST);
}

void glvxStencilToMask(glvxExtents extents) {
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glStencilMask(3);
	glStencilFunc(GL_NOTEQUAL, 2, 1);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

	glBegin(GL_QUADS);
		glVertex2f(extents[0], extents[1]);
		glVertex2f(extents[0], extents[3]);
		glVertex2f(extents[2], extents[3]);
		glVertex2f(extents[2], extents[1]);
	glEnd();
}

void glvxEnableMask() {
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glStencilMask(3);
	glStencilFunc(GL_EQUAL, 1, 3);
}

void glvxFill(size_t count, float *curves) {
	glvxExtents extents;

	/* Paint stencil buffer */
	glvxStencil(count, curves, extents);

	glvxPaint(extents);
}

void glvxFillRect(glvxExtents rect) {
	glvxStencilRect(rect);
	glvxPaint(rect);
}

void glvxFillMasked(size_t count, float *curves) {
	glvxExtents extents;

	glvxStencil(count, curves, extents);
	glvxEnableMask();
	glvxPaint(extents);
}

void glvxFillRectMasked(glvxExtents rect) {
	glvxStencilRect(rect);
	glvxEnableMask();
	glvxPaint(rect);
}

static inline void setGradientColor(size_t index) {
	index *= 4;
	glColor4f(gradientColors[index], gradientColors[index + 1], gradientColors[index + 2], gradientColors[index + 3]);
}

static void performSolidFill(glvxExtents extents) {
	glBegin(GL_QUADS);
		glVertex2f(extents[0], extents[1]);
		glVertex2f(extents[0], extents[3]);
		glVertex2f(extents[2], extents[3]);
		glVertex2f(extents[2], extents[1]);
	glEnd();
}

static void performGradientFill(glvxExtents extents) {
	if (gradientPoints < 1) {
		performSolidFill(extents);
		return;
	}	
	if (gradientPoints == 1) {
		setGradientColor(0);
		performSolidFill(extents);
		return;
	}

	float dx = -gradDeltaY;
	float dy = gradDeltaX;
	float deltaMag = (float)sqrt(dx * dx + dy * dy);
	dx /= deltaMag;
	dy /= deltaMag;

	Vector pointDiff, startDiff, endDiff;
	/* Comparing >= 0 to > 0 allows this to work also with horizontal and vertical gradients */
	if ((gradDeltaX >= 0) != (gradDeltaY > 0)) {
		pointDiff = make(gradBeginX - extents[0], gradBeginY - extents[1]);
		startDiff = make(-extents[0], -extents[3]);
		endDiff = make(-extents[2], -extents[1]);
	}
	else {
		pointDiff = make(gradBeginX - extents[2], gradBeginY - extents[1]);
		startDiff = make(-extents[0], -extents[1]);
		endDiff = make(-extents[2], -extents[3]);
	}

	Vector grad = make(gradDeltaX, gradDeltaY);
	Vector gradNormal = make(dx, dy);
	float dist = crossLen(grad, pointDiff) / length(grad);
	float sign = 1;
	if (gradDeltaX < 0) sign = -1;
	
	float startX = gradBeginX - dx * dist * sign;
	float startY = gradBeginY - dy * dist * sign;

	glBegin(GL_QUAD_STRIP);

	float pointSize = length(make(extents[2] - extents[0], extents[3] - extents[1]));

	{
		startDiff = add(make(startX, startY), startDiff);
		float startDist = crossLen(gradNormal, startDiff) / length(gradNormal);

		float x = startX + gradDeltaX * startDist * sign / deltaMag;
		float y = startY + gradDeltaY * startDist * sign / deltaMag;

		setGradientColor(0);
		glVertex2f(x, y);
		glVertex2f(x + dx * pointSize, y + dy * pointSize);
	}
	for (size_t i = 0; i < gradientPoints; ++i) {
		float x = startX + gradDeltaX * gradientTimes[i];
		float y = startY + gradDeltaY * gradientTimes[i];

		setGradientColor(i);
		glVertex2f(x, y);
		glVertex2f(x + dx * pointSize, y + dy * pointSize);
	}
	{
		endDiff = add(make(startX, startY), endDiff);
		float endDist = crossLen(gradNormal, endDiff) / length(gradNormal);

		float x = startX + gradDeltaX * endDist * sign / deltaMag;
		float y = startY + gradDeltaY * endDist * sign / deltaMag;

		setGradientColor(gradientPoints - 1);
		glVertex2f(x, y);
		glVertex2f(x + dx * pointSize, y + dy * pointSize);
	}

	glEnd();
}

void glvxStencilRect(glvxExtents extents) {
	/* Setup stencil */
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glStencilMask(1);
	glStencilFunc(GL_ALWAYS, 0, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);

	performSolidFill(extents);

	/* Enable color mask and set stencil pass function */
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glStencilFunc(GL_EQUAL, 1, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
}

void glvxPaint(glvxExtents extents) {
	if (fillMode == FILL_MODE_GRADIENT) {
		performGradientFill(extents);
	}
	else {
		performSolidFill(extents);
	}

	/* Setup stencil to zero out */
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glStencilFunc(GL_ALWAYS, 0, 0xff);
	glStencilOp(GL_ZERO, GL_ZERO, GL_ZERO);

	/* Zero out the stencil over the bounding box */
	performSolidFill(extents);

	/* Re-enable color mask */
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	glvxDisableStencil();
}

