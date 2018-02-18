#ifdef GLVX_IMPLEMENTATION

static inline void conversionPassthrough(glvxColor in, glvxColor out) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
	out[3] = in[3];
}

static inline float min_(float a, float b) {
	return a > b ? b : a;
}

static inline float max_(float a, float b) {
	return a > b ? a : b;
}

static inline void conversionHSV(glvxColor in, glvxColor out) {
	float minv, maxv, delta;

	out[3] = in[3];

	minv = min_(min_(in[0], in[1]), in[2]);
	maxv = max_(max_(in[0], in[1]), in[2]);
	out[2] = maxv;

	delta = maxv - minv;

	if (maxv != 0)
		out[1] = delta / maxv;
	else {
		out[1] = 0;
		out[0] = 0;
		return;
	}

	if (in[0] == maxv)
		out[0] = (in[1]- in[2]) / delta;
	else if (in[1]== maxv)
		out[0] = 2 + (in[2]- in[0]) / delta;
	else
		out[0] = 4 + (in[0] - in[1]) / delta;

	out[0] *= 60;
	if (out[0] < 0)
		out[0] += 360;

	out[0] /= 360;
}

#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif