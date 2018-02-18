#ifndef GLVX_H_
#define GLVX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/* A fundamental unit in glvx.
 * Represents a cubic bezier curve.
 * The first four elements represent the x parameter coefficients, while the second four
 * represent the y parameter coefficients.
 * So the total parametric equation is:
 * x(t) = curve[0] * t^3 + curve[1] * t^2 + curve[2] * t + curve[3] 
 * y(t) = curve[4] * t^3 + curve[5] * t^2 + curve[6] * t + curve[7]
 */
typedef float glvxCurve[8];

/* Represents a color; used in functions that need to dynamically update color.
 * Layout is R, G, B, A.
 */
typedef float glvxColor[4];

/* Represents bounds of some sort; generally either minimum point and maximum point 
 * or point and size.
 */
typedef float glvxExtents[4];

/* Sets the global glvx sample ratio.
 * Any kind of curve that glvx draws calculates some number of samples, usually based
 * on an estimation of how many are necessary before more would become redundant.
 * These calculations are based on pixel measurements, and so the ratio needs to 
 * be changed, if, for example, glScalef() is called.
 * The sample ratio is simply multiplied onto the calculated number of curve samples.
 * It therefore begins at 1 -- total samples is equal to 100% of pixel-based sample
 * estimates.
 * Generally, the sample ratio should correspond to the scale at which objects are
 * being drawn. So if drawing at twice the size, the sample ratio should be set
 * to 2, and if drawing at half the size, the sample ratio should be set to 0.5.
 */
extern void glvxSampleRatio(float newSampleRatio);

/* Sets the miter limit. */
extern void glvxMiterLimit(float newMiterLimit);
/* Requests the currently set global sample ratio.
 * Initially returns 1.
 */
extern float glvxGetSampleRatio();

/* Strokes the specified curve with the specified width.
 * Does not affect current color or stencil masking, and so can be used with
 * either.
 */
extern void glvxStroke(glvxCurve curve);

extern void glvxLeftBeginColor(glvxColor);
extern void glvxLeftBeginColor3(float r, float g, float b);
extern void glvxLeftBeginColor4(float r, float g, float b, float a);

extern void glvxRightBeginColor(glvxColor);
extern void glvxRightBeginColor3(float r, float g, float b);
extern void glvxRightBeginColor4(float r, float g, float b, float a);

extern void glvxLeftEndColor(glvxColor);
extern void glvxLeftEndColor3(float r, float g, float b);
extern void glvxLeftEndColor4(float r, float g, float b, float a);

extern void glvxRightEndColor(glvxColor);
extern void glvxRightEndColor3(float r, float g, float b);
extern void glvxRightEndColor4(float r, float g, float b, float a);

extern void glvxLeftColor(glvxColor leftColor);
extern void glvxLeftColor3(float r, float g, float b);
extern void glvxLeftColor4(float r, float g, float b, float a);
extern void glvxRightColor(glvxColor rightColor);
extern void glvxRightColor3(float r, float g, float b);
extern void glvxRightColor4(float r, float g, float b, float a);
extern void glvxBeginColor(glvxColor beginColor);
extern void glvxBeginColor3(float r, float g, float b);
extern void glvxBeginColor4(float r, float g, float b, float a);
extern void glvxEndColor(glvxColor endColor);
extern void glvxEndColor3(float r, float g, float b);
extern void glvxEndColor4(float r, float g, float b, float a);

extern void glvxBeginWidth(float width);
extern void glvxEndWidth(float width);
extern void glvxStrokeOffset(float offset);

extern void glvxEvenOdd(size_t count, float *points);

/* Sets up the stencil buffer for filling the specified curve.
 * 'curves' is an array of floating-point values, with at least 8 * 'count' members,
 * where every 8 values are one glvxCurve. glvxPaintMask takes those curves and, wherever
 * the shape bounded by the curves should be filled, it:
 * - Writes a 1 to the last bit if the current stencil value there has a last bit of 0
 * - Writes a 0 to the last bit if the current stencil value there has a last bit of 1
 * If the curve can be drawn, 'extents' is overwritten with the bounding box of the curve,
 * so that it can be filled in, for example by drawing a rectangle over the whole bounding
 * box.
 * After glvxStencil() returns, any drawing operations used will be affected by the mask.
 * Generally glvxClearMask should be called sometime after glvxPaintMask to disable the
 * stencil buffer.
 */
extern void glvxStencil(size_t count, float *curves, glvxExtents extents);
extern void glvxStencilRect(glvxExtents);

/* Prevents subsequent drawing operations from being affected by the mask from glStencil(). 
 * Does not do any overwrite of the stencil buffer.
 */
extern void glvxDisableStencil();


/* Fills the specified shape with whatever the current color is. 
 * 'curves' is an array of floating-point values, with at least 8 * 'count' members,
 * where every 8 values are one glvxCurve. glvxFill fills the shape bounded by these
 * curves without affecting the current GL color, but with overwriting the stencil.
 * If some pixel had a stencil value ending with 1 before this function was called,
 * it will:
 * - Not be filled if it was inside the curve
 * - Not be filled if it is outside the curve's bounding box
 * - Be filled if it is inside the curve's bounding box but outside the curve
 * Will also overwrite all stencil values inside the bounding box of the curve
 * with 0.
 */
extern void glvxFill(size_t count, float *curves);
extern void glvxFillRect(glvxExtents);

extern void glvxFillMasked(size_t count, float *curves);
extern void glvxFillRectMasked(glvxExtents);

extern void glvxStencilToMask(glvxExtents);
extern void glvxEnableMask();

extern void glvxPaint(glvxExtents);

/* Cleans up internal memory used by glvx.
 * Should be called on application exit.
 */
extern void glvxCleanup();
extern int glvxInit();

extern void glvxUseHSV();
extern void glvxUseRGB();

/* Overwrites the specified glvxCurve with the curve between endpoints
 * (location[0], location[1]) and (location[2], location[3]), with ease-in and ease-out
 * handles of (ease[0], ease[1]) and (ease[2], ease[3]).
 */
extern void glvxCalculateCurve(glvxCurve curve, float location[4], float ease[4]);

/* Overwrites the specified glvxExtents object with the bounding box of
 * the specified glvxCurve (in terms of minima and maxima).
 * extents[0] = Minimum curve x
 * extents[1] = Minimum curve y
 * extents[2] = Maximum curve x
 * extents[3] = Maximum curve y
 */
extern void glvxGetExtents(glvxCurve curve, glvxExtents extents);

/* Overwrites the specified glvxExtents object with the bounding box of
 * the specified glvxCurve (in terms of width and height).
 * extents[0] = Curve x
 * extents[1] = Curve y
 * extents[2] = Curve width
 * extents[3] = Curve height
 */
extern void glvxGetBounds(glvxCurve curve, glvxExtents bounds);

/* Gets the time at which the curvature of the specified curve goes to zero.
 * Not all curves have this time; generally the ones that do are the ones with
 * two "bumps." The point inbetween those bumps is the zero curvature point,
 * and the time when the curve reaches that point is the zero curvature time.
 * If the time exists, it is written into 'time' and the function returns 1.
 * If the time does not exist, 'time' is not modified and the funciton
 * returns 0.
 */
extern int glvxGetZeroCurvatureTime(glvxCurve curve, float *time);

/* Gets the point at which the curvature of the specified curve goes to zero.
* Not all curves have this point; generally the ones that do are the ones with
* two "bumps." The point inbetween those bumps is the zero curvature point,
* and the time when the curve reaches that point is the zero curvature time.
* If the point exists, it is written into 'x' and 'y' and the functino returns 1.
* If the time does not exist, 'x' and 'y' are not modified and the funciton
* returns 0.
*/
extern int glvxGetZeroCurvaturePoint(glvxCurve c, float *x, float *y);

/* Gets the x coordinate along the curve at the specified time. */
extern float glvxCurveX(glvxCurve, float time);

/* Gets the y coordinate along the curve at the specified time. */
extern float glvxCurveY(glvxCurve, float time);

/* Gets the starting x coordinate of the curve. */
extern float glvxCurveXBegin(glvxCurve);

/* Gets the starting y coordinate of the curve. */
extern float glvxCurveYBegin(glvxCurve);

/* Gets the ending x coordinate of the curve. */
extern float glvxCurveXEnd(glvxCurve);

/* Gets the ending y coordinate of the curve. */
extern float glvxCurveYEnd(glvxCurve);

extern void glvxFillModeGradient();
extern void glvxFillModeSolid();

extern void glvxGradientBegin(float x, float y);
extern void glvxGradientEnd(float x, float y);
extern void glvxGradientDirection(float x, float y);
extern void glvxGradient(size_t count, float *colors, float *times);

#ifdef __cplusplus
}
#endif

#endif