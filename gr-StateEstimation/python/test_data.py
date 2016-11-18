#!usr/bin/env python

def main():
	#input data
	size = 512
	signal = []
	f = open("input_test_data.txt", 'r')
	for i in range(1, size+1):
		line = f.readline()
		signal.append(float(line))
	f.close()
	input_signal = tuple(signal)
	#src = gr.vector_source_f(input_signal, False, size)
	
	print input_signal

if __name__ == "__main__":
	main()