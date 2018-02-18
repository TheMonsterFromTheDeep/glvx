#ifdef GLVX_IMPLEMENTATION

static inline GLuint compileShader(const GLchar *vertexSource, const GLchar *fragmentSource) {
	GLint vertexLength = (GLint)strlen(vertexSource);
	GLint fragmentLength = (GLint)strlen(fragmentSource);

	GLint vertexShader, fragmentShader;
	vertexShader = glCreateShader(GL_VERTEX_SHADER);
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(vertexShader, 1, &vertexSource, &vertexLength);
	glShaderSource(fragmentShader, 1, &fragmentSource, &fragmentLength);

	glCompileShaderARB(vertexShader);
	glCompileShaderARB(fragmentShader);

	GLint compiled;

	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled) {
		return 0;
	}

	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled) {
		return 0;
	}

	GLuint shaderProgram;
	shaderProgram = glCreateProgram();

	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);

	glLinkProgram(shaderProgram);
	
	return shaderProgram;
}

#else
	#error "This file is part of the GLVX implementation. It is not an include file."
#endif