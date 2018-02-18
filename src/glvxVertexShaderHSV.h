"#version 110\n"

"vec4 rgbToHsv(in vec4 rgba) {"
	"vec4 res;"
	"res.w = rgba.w;"
	"float minv, maxv, delta;"

	"minv = min(min(rgba.x, rgba.y), rgba.z);"
	"maxv = max(max(rgba.x, rgba.y), rgba.z);"
	"res.z = maxv;"

	"delta = maxv - minv;"

	"if (maxv != 0.f) res.y = delta / maxv;"
	"else {"
		"res.y = 0.f;"
		"res.x = 0.f;"
		"return res;"
	"}"

	"if (rgba.x == maxv) { res.x = (rgba.y - rgba.z) / delta; }"
	"else if (rgba.y == maxv) { res.x = 2.f + (rgba.z - rgba.x) / delta; }"
	"else { res.x = 4.f + (rgba.x - rgba.y) / delta; }"

	"res.x *= 60.f;"
	"if (res.x < 0.f) { res.x += 360.f; }"

	"res.x /= 360.f;"

	"return res;"
"}"

"void main(void) {"
	"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;"
	"gl_FrontColor = rgbToHsv(gl_Color);"
	"gl_BackColor = gl_FrontColor;"
"}"