#ifdef GLVX_IMPLEMENTATION

#define FILL_MODE_SOLID 0
#define FILL_MODE_GRADIENT 1

#define INTERPOLATION_RGB 0
#define INTERPOLATION_HSV 1

static float sampleRatio = 1.f;
static float miterLimit = 5.f;

static glvxColor leftBeginColor = { 0.f, 0.f, 0.f, 1.f };
static glvxColor rightBeginColor = { 0.f, 0.f, 0.f, 1.f };
static glvxColor leftEndColor = { 0.f, 0.f, 0.f, 1.f };
static glvxColor rightEndColor = { 0.f, 0.f, 0.f, 1.f };

static float beginWidth = 1.f;
static float endWidth = 1.f;
static float strokeOffset = 0.f;

static float gradBeginX = 0;
static float gradBeginY = 0;
static float gradDeltaX = 1;
static float gradDeltaY = 0;

static int fillMode = FILL_MODE_SOLID;

static int interpolationMode = INTERPOLATION_RGB;

static float *gradientTimes;
static float *gradientColors;
static size_t gradientPoints = 0;

#define SET_COLOR(dest, src)\
	dest[0] = src[0];\
	dest[1] = src[1];\
	dest[2] = src[2];\
	dest[3] = src[3]

#define SET_COLOR_RGB(dest)\
	dest[0] = r;\
	dest[1] = g;\
	dest[2] = b;\
	dest[3] = 1.f

#define SET_COLOR_RGBA(dest)\
	dest[0] = r;\
	dest[1] = g;\
	dest[2] = b;\
	dest[3] = a

void glvxSampleRatio(float newSampleRatio) {
	sampleRatio = newSampleRatio;
}

float glvxGetSampleRatio() {
	return sampleRatio;
}

void glvxMiterLimit(float newMiterLimit) {
	miterLimit = newMiterLimit;
}

void glvxLeftBeginColor(glvxColor newColor) { SET_COLOR(leftBeginColor, newColor); }
void glvxLeftBeginColor3(float r, float g, float b) { SET_COLOR_RGB(leftBeginColor); }
void glvxLeftBeginColor4(float r, float g, float b, float a) { SET_COLOR_RGBA(leftBeginColor); }

void glvxRightBeginColor(glvxColor newColor) { SET_COLOR(rightBeginColor, newColor); }
void glvxRightBeginColor3(float r, float g, float b) { SET_COLOR_RGB(rightBeginColor); }
void glvxRightBeginColor4(float r, float g, float b, float a) { SET_COLOR_RGBA(rightBeginColor); }

void glvxLeftEndColor(glvxColor newColor) { SET_COLOR(leftEndColor, newColor); }
void glvxLeftEndColor3(float r, float g, float b) { SET_COLOR_RGB(leftEndColor); }
void glvxLeftEndColor4(float r, float g, float b, float a) { SET_COLOR_RGBA(leftEndColor); }

void glvxRightEndColor(glvxColor newColor) { SET_COLOR(rightEndColor, newColor);  }
void glvxRightEndColor3(float r, float g, float b) { SET_COLOR_RGB(rightEndColor); }
void glvxRightEndColor4(float r, float g, float b, float a) { SET_COLOR_RGBA(rightEndColor); }

void glvxLeftColor(glvxColor newColor) {
	glvxLeftBeginColor(newColor);
	glvxLeftEndColor(newColor);
}

void glvxRightColor(glvxColor newColor) {
	glvxRightBeginColor(newColor);
	glvxRightEndColor(newColor);
}

void glvxBeginColor(glvxColor newColor) {
	glvxLeftBeginColor(newColor);
	glvxRightBeginColor(newColor);
}

void glvxEndColor(glvxColor newColor) {
	glvxLeftEndColor(newColor);
	glvxRightEndColor(newColor);
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

void glvxUseHSV() {
	if (GLEW_VERSION_1_5) {
		interpolationMode = INTERPOLATION_HSV;
		glUseProgram(hsvShader);
	}
}

void glvxUseRGB() {
	if (GLEW_VERSION_1_5) {
		interpolationMode = INTERPOLATION_RGB;
		glUseProgram(0);
	}
}

#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif