#include <stdio.h>
#include <string.h>
#include <math.h>

#define GLEW_STATIC
#include <GL/glew.h>

#include "GLVX.h"

#define PI 3.14159265358979323846264338327950288

float *curveEvenOddList = NULL;
size_t curveEvenOddListSize = 0;

static GLuint hsvShader = 0;
static GLuint strokeHsvShader = 0;

#define GLVX_IMPLEMENTATION
#include "glvxInternalParameters.h"
#include "glvxInternalCurveProperties.h"
#include "glvxInternalVector.h"
#include "glvxInternalUtil.h"
#include "glvxInternalShaderLoader.h"
#include "glvxInternalColorConversion.h"

int glvxInit() {
	if (GLEW_OK != glewInit()) {
		return 1;
	}

	if (GLEW_VERSION_1_5) {
		hsvShader = compileShader(
			#include "glvxVertexShaderHSV.h"
			,
			#include "glvxFragmentShaderHSV.h"
		);
		strokeHsvShader = compileShader(
			"#version 110\n"
			"void main(void) {"
				"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;"
				"gl_FrontColor = gl_Color;"
				"gl_BackColor = gl_FrontColor;"
			"}",
			#include "glvxFragmentShaderHSV.h"
		);
	}

	return 0;
}

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

void glvxGetCircleExtents(float x, float y, float r, glvxExtents extents) {
	if (extents) {
		extents[0] = x - r;
		extents[1] = y - r;
		extents[2] = x + r;
		extents[3] = y + r;
	}
}

void glvxGetCircleBounds(float x, float y, float r, glvxExtents bounds) {
	if (!bounds) return;
	glvxGetCircleExtents(x, y, r, bounds);
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

	float leftColorStep[4];
	float leftColor[4];
	float rightColorStep[4];
	float rightColor[4];

	glvxColor myLeftBegin;
	glvxColor myLeftEnd;
	glvxColor myRightBegin;
	glvxColor myRightEnd;

	if (interpolationMode == INTERPOLATION_HSV) {
		conversionHSV(leftBeginColor, myLeftBegin);
		conversionHSV(leftEndColor, myLeftEnd);
		conversionHSV(rightBeginColor, myRightBegin);
		conversionHSV(rightEndColor, myRightEnd);

		if (GLEW_VERSION_1_5) {
			glUseProgram(strokeHsvShader);
		}
	}
	else {
		conversionPassthrough(rightBeginColor, myLeftBegin);
		conversionPassthrough(rightEndColor, myLeftEnd);
		conversionPassthrough(rightBeginColor, myRightBegin);
		conversionPassthrough(rightEndColor, myRightEnd);
	}

	leftColor[0] = myLeftBegin[0];
	leftColor[1] = myLeftBegin[1];
	leftColor[2] = myLeftBegin[2];
	leftColor[3] = myLeftBegin[3];

	leftColorStep[0] = (myLeftEnd[0] - myLeftBegin[0]) / samples;
	leftColorStep[1] = (myLeftEnd[1] - myLeftBegin[1]) / samples;
	leftColorStep[2] = (myLeftEnd[2] - myLeftBegin[2]) / samples;
	leftColorStep[3] = (myLeftEnd[3] - myLeftBegin[3]) / samples;

	rightColor[0] = myRightBegin[0];
	rightColor[1] = myRightBegin[1];
	rightColor[2] = myRightBegin[2];
	rightColor[3] = myRightBegin[3];

	rightColorStep[0] = (myRightEnd[0] - myRightBegin[0]) / samples;
	rightColorStep[1] = (myRightEnd[1] - myRightBegin[1]) / samples;
	rightColorStep[2] = (myRightEnd[2] - myRightBegin[2]) / samples;
	rightColorStep[3] = (myRightEnd[3] - myRightBegin[3]) / samples;

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

	glColor4f(leftColor[0], leftColor[1], leftColor[2], leftColor[3]);
	glVertexVector(add(previous, multiply(orth, width * (1 + strokeOffset))));
	glColor4f(rightColor[0], rightColor[1], rightColor[2], rightColor[3]);
	glVertexVector(subtract(previous, multiply(orth, width * (1 - strokeOffset))));
	
	while(samples > 1) {
		sampleValue += sampleSize;
		width += widthStep;
		leftColor[0] += leftColorStep[0];
		leftColor[1] += leftColorStep[1];
		leftColor[2] += leftColorStep[2];
		leftColor[3] += leftColorStep[3];
		rightColor[0] += rightColorStep[0];
		rightColor[1] += rightColorStep[1];
		rightColor[2] += rightColorStep[2];
		rightColor[3] += rightColorStep[3];
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
		glColor4f(leftColor[0], leftColor[1], leftColor[2], leftColor[3]);
		glVertexVector(add(current, multiply(proj, 1 + strokeOffset)));
		glColor4f(rightColor[0], rightColor[1], rightColor[2], rightColor[3]);
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
	glColor4f(leftColor[0], leftColor[1], leftColor[2], leftColor[3]);
	glVertexVector(add(current, multiply(orth, 1 + strokeOffset)));
	glColor4f(rightColor[0], rightColor[1], rightColor[2], rightColor[3]);
	glVertexVector(subtract(current, multiply(orth, 1 - strokeOffset)));

	glEnd();

	if(interpolationMode == INTERPOLATION_HSV && GLEW_VERSION_1_5) {
		/* Revert to normal hsv shader */
		glUseProgram(hsvShader);
	}
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

static inline void stencilWriteSetup() {
	glEnable(GL_STENCIL_TEST);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glStencilMask(1);
	glStencilFunc(GL_ALWAYS, 0, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);
}

static inline void stencilUseSetup() {
	/* Enable color mask and set stencil pass function */
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glStencilFunc(GL_EQUAL, 1, 0xff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
}

void glvxStencilRect(glvxExtents extents) {
	stencilWriteSetup();
	performSolidFill(extents);
	stencilUseSetup();
}

void glvxStencilRoundRect(glvxExtents extents, float r) {
	if (r < 0) r = -r;

	stencilWriteSetup();

	size_t steps = (size_t)ceil(PI * 2 * r * 0.25f * sampleRatio);
	float stepSize = (float)(PI / (steps * 2)); /* 2pi / (steps * 4) */

	glBegin(GL_TRIANGLE_FAN);

	float left = extents[0] + r;
	float right = extents[2] - r;
	float bottom = extents[1] + r;
	float top = extents[3] - r;

	for (size_t i = 0; i <= steps; ++i) {
		float t = stepSize * i;
		glVertex2d(right + r * cos(t), top + r * sin(t));
	}

	for (size_t i = steps; i <= steps * 2; ++i) {
		float t = stepSize * i;
		glVertex2d(left + r * cos(t), top + r * sin(t));
	}

	for (size_t i = steps * 2; i <= steps * 3; ++i) {
		float t = stepSize * i;
		glVertex2d(left + r * cos(t), bottom + r * sin(t));
	}

	for (size_t i = steps * 3; i <= steps * 4; ++i) {
		float t = stepSize * i;
		glVertex2d(right + r * cos(t), bottom + r * sin(t));
	}

	glEnd();

	stencilUseSetup();
}

void glvxStencilCircle(float x, float y, float r, glvxExtents extents) {
	if (r < 0) r = -r;

	stencilWriteSetup();

	size_t steps = (size_t)ceil(PI * 2 * r * sampleRatio); /* Circumference */
	float stepSize = (float)((PI * 2) / steps);

	glBegin(GL_TRIANGLE_FAN);

	for (size_t i = 0; i <= steps; ++i) {
		float t = stepSize * i;
		glVertex2d(x + r * cos(t), y + r * sin(t));
	}

	glEnd();

	stencilUseSetup();

	glvxGetCircleExtents(x, y, r, extents);
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

void glvxFillRoundRect(glvxExtents rect, float r) {
	glvxStencilRoundRect(rect, r);
	glvxPaint(rect);
}

void glvxFillCircle(float x, float y, float r) {
	glvxExtents extents;
	glvxStencilCircle(x, y, r, extents);
	glvxPaint(extents);
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

void glvxFillRoundRectMasked(glvxExtents rect, float r) {
	glvxStencilRoundRect(rect, r);
	glvxEnableMask();
	glvxPaint(rect);
}

void glvxFillCircleMasked(float x, float y, float r) {
	glvxExtents extents;
	glvxStencilCircle(x, y, r, extents);
	glvxEnableMask();
	glvxPaint(extents);
}