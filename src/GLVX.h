#ifndef GLVX_H_
#define GLVX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

typedef float glvxCurve[8];
typedef float glvxColor[4];
typedef float glvxExtents[4];

extern void glvxSampleRatio(float newSampleRatio);
extern float glvxGetSampleRatio();

extern void glvxStroke(glvxCurve curve, float width);
extern void glvxStrokew(glvxCurve curve, float widthStart, float widthEnd);
extern void glvxStrokec(glvxCurve curve, float width, glvxColor color0, glvxColor color1);
extern void glvxStrokewc(glvxCurve curve, float widthStart, float widthEnd, glvxColor color0, glvxColor color1);

extern void glvxEarcut(size_t points, float *polygon);
extern void glvxPaintMask(size_t count, float *curves, glvxExtents extents);
extern void glvxClearMask();
extern void glvxFill(size_t count, float *curves);
extern void glvxCleanup();

extern void glvxCalculateCurve(glvxCurve curve, float location[4], float ease[4]);

void glvxGetExtents(glvxCurve curve, glvxExtents extents);
void glvxGetBounds(glvxCurve curve, glvxExtents bounds);
extern int glvxGetZeroCurvatureTime(glvxCurve curve, float *time);
extern int glvxGetZeroCurvaturePoint(glvxCurve c, float *x, float *y);

#ifdef __cplusplus
}
#endif

#endif