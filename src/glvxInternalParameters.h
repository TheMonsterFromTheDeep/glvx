#ifdef GLVX_IMPLEMENTATION

#define FILL_MODE_SOLID 0
#define FILL_MODE_GRADIENT 1

static float sampleRatio = 1.f;
static float miterLimit = 5.f;

static glvxColor leftColor = { 1.f, 1.f, 1.f, 1.f };
static glvxColor rightColor = { 1.f, 1.f, 1.f, 1.f };
static glvxColor beginColor = { 0.f, 0.f, 0.f, 1.f };
static glvxColor endColor = { 0.f, 0.f, 0.f, 1.f };

static float beginWidth = 1.f;
static float endWidth = 1.f;
static float strokeOffset = 0.f;

static float gradBeginX = 0;
static float gradBeginY = 0;
static float gradDeltaX = 1;
static float gradDeltaY = 0;

static int fillMode = FILL_MODE_SOLID;

static float *gradientTimes;
static float *gradientColors;
static size_t gradientPoints = 0;

#define SET_COLOR(dest, src)\
	dest[0] = src[0];\
	dest[1] = src[1];\
	dest[2] = src[2];\
	dest[3] = src[3]

void glvxSampleRatio(float newSampleRatio) {
	sampleRatio = newSampleRatio;
}

float glvxGetSampleRatio() {
	return sampleRatio;
}

void glvxMiterLimit(float newMiterLimit) {
	miterLimit = newMiterLimit;
}

void glvxLeftColor(glvxColor newColor) {
	SET_COLOR(leftColor, newColor);
}

void glvxRightColor(glvxColor newColor) {
	SET_COLOR(rightColor, newColor);
}

void glvxBeginColor(glvxColor newColor) {
	SET_COLOR(beginColor, newColor);
}

void glvxEndColor(glvxColor newColor) {
	SET_COLOR(endColor, newColor);
}

void glvxBeginWidth(float newWidth) {
	/* Because we use half width anyways, calculate that here */
	beginWidth = newWidth * 0.5f;
}

void glvxEndWidth(float newWidth) {
	endWidth = newWidth * 0.5f;
}

void glvxStrokeOffset(float newOffset) {
	strokeOffset = newOffset;
}

void glvxFillModeGradient() { fillMode = FILL_MODE_GRADIENT; }

void glvxFillModeSolid() { fillMode = FILL_MODE_SOLID; }

void glvxGradient(size_t count, float *colors, float *times) {
	gradientTimes = times;
	gradientColors = colors;
	gradientPoints = count;
}

void glvxGradientBegin(float x, float y) {
	gradBeginX = x;
	gradBeginY = y;
}

void glvxGradientEnd(float x, float y) {
	gradDeltaX = x - gradBeginX;
	gradDeltaY = y - gradBeginY;
}

void glvxGradientDirection(float x, float y) {
	gradDeltaX = x;
	gradDeltaY = y;
}

#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif