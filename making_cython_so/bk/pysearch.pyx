cdef extern from "search.h":
	int main_gut(int st, int se)

def py_main_gut(int st, int se):
	main_gut(st, se)
