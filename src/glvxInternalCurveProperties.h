#ifdef GLVX_IMPLEMENTATION

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
#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif