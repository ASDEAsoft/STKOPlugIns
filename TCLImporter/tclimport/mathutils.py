import math

def v3norm(x):
	return math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

def v3normalize(x):
	n = v3norm(x)
	if n > 0.0:
		return (x[0]/n, x[1]/n, x[2]/n)
	else:
		return (0.0, 0.0, 0.0)

def v3cross(x, y):
	return (
		x[1]*y[2] - x[2]*y[1],
		x[2]*y[0] - x[0]*y[2],
		x[0]*y[1] - x[1]*y[0])

def v3sub(x, y):
	return (x[0]-y[0], x[1]-y[1], x[2]-y[2])

def v3dot(x, y):
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]