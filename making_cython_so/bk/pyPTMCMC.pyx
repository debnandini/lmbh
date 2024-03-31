cdef extern from "PTMCMC.h":
	int main_gut(int seg, int rep, int delf)

def py_main_gut(int seg, int rep, int delf):
	main_gut(seg, rep, delf)
