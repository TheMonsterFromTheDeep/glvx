#ifdef GLVX_IMPLEMENTATION

void glvxStrokew(glvxCurve c, float widthStart, float widthEnd) {
	size_t samples = getCurveSamples(c);
	if (samples < 2) return;
	float sampleSize = 1.f / samples;
	float sampleValue = sampleSize;

	widthEnd = (widthEnd - widthStart) / (samples * 2);
	widthStart /= 2;

	glBegin(GL_QUAD_STRIP);

	Vector delta, orth, bi, proj;

	Vector previous = curveVector(c, 0);
	Vector current = curveVector(c, sampleValue);

	delta = subtract(current, previous);
	orth = orthogonalTo(delta);
	orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));

	glVertexVector(add(previous, multiply(orth, widthStart)));
	glVertexVector(subtract(previous, multiply(orth, widthStart)));

	while (samples > 1) {
		sampleValue += sampleSize;
		widthStart += widthEnd;
		Vector next = curveVector(c, sampleValue);

		Vector nextDelta = subtract(current, next);

		bi = bisector(delta, nextDelta);

		proj = divide(multiply(bi, widthStart), dot(orth, bi));

		glVertexVector(add(current, proj));
		glVertexVector(subtract(current, proj));

		delta = multiply(nextDelta, -1);
		orth = orthogonalTo(delta);
		orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));
		current = next;
		previous = current;

		--samples;
	}

	orth = multiply(orth, widthStart);

	current = curveVector(c, sampleValue);
	glVertexVector(add(current, orth));
	glVertexVector(subtract(current, orth));

	glEnd();
}

void glvxStrokec(glvxCurve c, float width, glvxColor color0, glvxColor color1) {
	size_t samples = getCurveSamples(c);
	if (samples < 2) return;
	float sampleSize = 1.f / samples;
	float sampleValue = sampleSize;

	/* We center the curve, and therefore use the half-width for projection */
	width /= 2;

	float colorStep[4];
	float color[4];

	color[0] = color0[0];
	color[1] = color0[1];
	color[2] = color0[2];
	color[3] = color0[3];

	colorStep[0] = (color1[0] - color0[0]) / samples;
	colorStep[1] = (color1[1] - color0[1]) / samples;
	colorStep[2] = (color1[2] - color0[2]) / samples;
	colorStep[3] = (color1[3] - color0[3]) / samples;

	glBegin(GL_QUAD_STRIP);

	Vector delta, orth, bi, proj;

	Vector previous = curveVector(c, 0);
	Vector current = curveVector(c, sampleValue);

	delta = subtract(current, previous);
	orth = orthogonalTo(delta);
	orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));

	glColor4f(color[0], color[1], color[2], color[3]);

	glVertexVector(add(previous, multiply(orth, width)));
	glVertexVector(subtract(previous, multiply(orth, width)));

	while (samples > 1) {
		sampleValue += sampleSize;
		color[0] += colorStep[0];
		color[1] += colorStep[1];
		color[2] += colorStep[2];
		color[3] += colorStep[3];
		glColor4f(color[0], color[1], color[2], color[3]);

		Vector next = curveVector(c, sampleValue);

		Vector nextDelta = subtract(current, next);

		bi = bisector(delta, nextDelta);

		proj = divide(multiply(bi, width), dot(orth, bi));

		glVertexVector(add(current, proj));
		glVertexVector(subtract(current, proj));

		delta = multiply(nextDelta, -1);
		orth = orthogonalTo(delta);
		orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));
		current = next;
		previous = current;

		--samples;
	}

	orth = multiply(orth, width);

	current = curveVector(c, sampleValue);
	glVertexVector(add(current, orth));
	glVertexVector(subtract(current, orth));

	glEnd();
}

void glvxStrokewc(glvxCurve c, float widthStart, float widthEnd, glvxColor color0, glvxColor color1) {
	size_t samples = getCurveSamples(c);
	if (samples < 2) return;
	float sampleSize = 1.f / samples;
	float sampleValue = sampleSize;

	widthEnd = (widthEnd - widthStart) / (samples * 2);
	widthStart /= 2;

	float colorStep[4];
	float color[4];

	color[0] = color0[0];
	color[1] = color0[1];
	color[2] = color0[2];
	color[3] = color0[3];

	colorStep[0] = (color1[0] - color0[0]) / samples;
	colorStep[1] = (color1[1] - color0[1]) / samples;
	colorStep[2] = (color1[2] - color0[2]) / samples;
	colorStep[3] = (color1[3] - color0[3]) / samples;

	glBegin(GL_QUAD_STRIP);

	Vector delta, orth, bi, proj;

	Vector previous = curveVector(c, 0);
	Vector current = curveVector(c, sampleValue);

	delta = subtract(current, previous);
	orth = orthogonalTo(delta);
	orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));

	glVertexVector(add(previous, multiply(orth, widthStart)));
	glVertexVector(subtract(previous, multiply(orth, widthStart)));

	while (samples > 1) {
		sampleValue += sampleSize;
		widthStart += widthEnd;
		color[0] += colorStep[0];
		color[1] += colorStep[1];
		color[2] += colorStep[2];
		color[3] += colorStep[3];
		glColor4f(color[0], color[1], color[2], color[3]);
		Vector next = curveVector(c, sampleValue);

		Vector nextDelta = subtract(current, next);

		bi = bisector(delta, nextDelta);

		proj = divide(multiply(bi, widthStart), dot(orth, bi));

		glVertexVector(add(current, proj));
		glVertexVector(subtract(current, proj));

		delta = multiply(nextDelta, -1);
		orth = orthogonalTo(delta);
		orth = multiply(orth, fastInverseSquareRoot(squaredLength(orth)));
		current = next;
		previous = current;

		--samples;
	}

	orth = multiply(orth, widthStart);

	current = curveVector(c, sampleValue);
	glVertexVector(add(current, orth));
	glVertexVector(subtract(current, orth));

	glEnd();
}

#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif