all: matrix-search

matrix-search: clean
	g++ MatrixSearch.cpp -o matrix-search
	
clean:
	rm matrix-search -f
