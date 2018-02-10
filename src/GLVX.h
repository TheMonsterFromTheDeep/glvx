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
/* Requests the currently set global sample ratio.
 * Initially returns 1.
 */
extern float glvxGetSampleRatio();

/* Strokes the specified curve with the specified width.
 * Does not affect current color or stencil masking, and so can be used with
 * either.
 */
extern void glvxStroke(glvxCurve curve, float width);

/* Strokes the specified curve, varying width linearly between the specified endpoint widths.
 * The curve will have widthStart at its starting endpoint, and widthEnd at its closing endpoint.
 * Starting and closing endpoints can be polled with glvxCurveXBegin(), glvxCurveYBegin(),
 * glvxCurveXEnd(), and glvxCurveYEnd().
 * Does not affect current color or stencil masking, and so can be used with either.
 */
extern void glvxStrokew(glvxCurve curve, float widthStart, float widthEnd);

/* Strokes the specified curve with the specified width, varying color linearly between the specified endpoint colors.
 * The curve will have colorStart at its starting endpoint, and colorEnd at its closing endpoint.
 * Starting and closing endpoints can be polled with glvxCurveXBegin(), glvxCurveYBegin(),
 * glvxCurveXEnd(), and glvxCurveYEnd().
 * Sets OpenGL color at all sample points along the curve, to color linearly interpolated between the specified
 * endpoint colors.
 */
extern void glvxStrokec(glvxCurve curve, float width, glvxColor colorStart, glvxColor colorEnd);

/* Strokes the specified curve, varying with and color linearly between the specified endpoint widths and colors.
 * The curve will have widthStart and colorStart at its starting endpoint, and widthEnd and colorEnd at its closing endpoint.
 * Starting and closing endpoints can be polled with glvxCurveXBegin(), glvxCurveYBegin(),
 * glvxCurveXEnd(), and glvxCurveYEnd().
 * Sets OpenGL color at all sample points along the curve, to color linearly interpolated between the specified
 * endpoint colors.
 */
extern void glvxStrokewc(glvxCurve curve, float widthStart, float widthEnd, glvxColor colorStart, glvxColor colorEnd);

/* Fills the specified simple polygon.
 * The implementation of this method is not guaranteed to be earcut, as long as it behaves equivalently or better.
 * This means it is incorrect to call this method with a non-simple polygon and expect a particular result.
 * 'points' specifies the number of points in the polygon, and 'polygon' is a list of x, y pairs.
 * 'polygon' should have at least twice 'points' number of elements.
 */
extern void glvxEarcut(size_t points, float *polygon);

/* Sets up the stencil buffer for filling the specified curve.
 * 'curves' is an array of floating-point values, with at least 8 * 'count' members,
 * where every 8 values are one glvxCurve. glvxPaintMask takes those curves and fills
 * the stencil value with 1 wherever the shape bounded by those curves should be filled.
 * If the curve can be drawn, 'extents' is overwritten with the bounding box of the curve,
 * so that it can be filled in, for example by drawing a rectangle over the whole bounding
 * box.
 * After glvxPaintMask() returns, any drawing operations used will be affected by the mask.
 * Generally glvxClearMask should be called after glvxPaintMask() and before any other operations,
 * and anything else that uses the stencil buffer should not assume it is in any particular
 * state.
 */
extern void glvxPaintMask(size_t count, float *curves, glvxExtents extents);

/* Prevents subsequent drawing operations from being affected by the mask from glPaintMask(). */
extern void glvxClearMask();

/* Fills the specified shape with whatever the current color is. 
 * 'curves' is an array of floating-point values, with at least 8 * 'count' members,
 * where every 8 values are one glvxCurve. glvxFill fills the shape bounded by these
 * curves without affecting the current GL color, but with overwriting the stencil.
 * Users should make no assumptions about the state of the stencil after this
 * method returns, except that it will have no effect until set up again.
 */
extern void glvxFill(size_t count, float *curves);

/* Cleans up internal memory used by glvx.
 * Should be called on application exit.
 */
extern void glvxCleanup();

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

#ifdef __cplusplus
}
#endif

#endif