#ifdef GLVX_IMPLEMENTATION

typedef struct {
	float x, y;
} Vector;

static inline Vector make(float x, float y) { return (Vector) { x, y }; }

static inline Vector add(Vector a, Vector b) {
	return make(a.x + b.x, a.y + b.y);
}

static inline Vector subtract(Vector a, Vector b) {
	return make(a.x - b.x, a.y - b.y);
}

static inline Vector multiply(Vector a, float s) {
	return make(a.x * s, a.y * s);
}

static inline Vector divide(Vector a, float s) {
	return make(a.x / s, a.y / s);
}

static inline Vector orthogonalTo(Vector a) {
	return make(-a.y, a.x);
}

static inline float squaredLength(Vector v) {
	return v.x * v.x + v.y * v.y;
}

static inline float length(Vector v) {
	return (float)sqrt(v.x * v.x + v.y * v.y);
}

static inline Vector bisector(Vector a, Vector b) {
	Vector bisec = add(multiply(a, length(b)), multiply(b, length(a)));
	if (squaredLength(bisec) < 0.0000001f) return orthogonalTo(a);
	return bisec;
}

static inline float dot(Vector a, Vector b) {
	return a.x * b.x + a.y * b.y;
}

static inline void glVertexVector(Vector v) {
	glVertex2f(v.x, v.y);
}

static inline Vector curveVector(glvxCurve c, float t) {
	return make(
		((c[0] * t + c[1]) * t + c[2]) * t + c[3],
		((c[4] * t + c[5]) * t + c[6]) * t + c[7]
	);
}

#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif