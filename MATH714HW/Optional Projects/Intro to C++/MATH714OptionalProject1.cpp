// MATH741OptionalProject1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
using namespace std;


int PascalTriangle_Coefficients(int n, int k)
{
   /*
   Instructions:
    (a) Write a program that can evaluate the numbers $c_{n,k}$ in Pascal's triangle using a recurrence relation.
   Define $c_{n,0} = 1$ for all $n \geq 0$ and $c_{0,k}=0$ for $k \geq 0$. Then the entries for $n + 1$ and be computed from the entries for $n$ using 
   \[
		c_{n+1,k} = c_{n,k} + c_{n,k-1}.
	\]
	Print a table of entries in Pascal's triangle for $0 \leq n \leq 7$ and $0 \leq k \leq 7$. You should find that $c_{n,k} = 0$ for $k > n$, and those terms can be left blank in your table.


	(b) Extend your program so that it prints ASCII art, displaying a "." character for each even number and a "#" character for each odd number. Hence the first few lines would be
	#
	##
	#.#
	####
	Using your program, extend this output to $n = 31$.
   */ 
	if (k == 1) return 1;
	if (n == 0) return 0;
	return PascalTriangle_Coefficients(n - 1, k) + PascalTriangle_Coefficients(n - 1, k - 1);
}

vector<vector<int>> Generate_Pascal_Triangle(int numRows) {
	vector<vector<int>> triangle;
	for (int i = 0; i < numRows; i++) {
		vector<int> row(i + 1, 1); // Initialize a row with 1s
		for (int j = 1; j < i; j++) {
			row[j] = triangle[i - 1][j - 1] + triangle[i - 1][j]; // Fill in the values based on the previous row
		}
		triangle.push_back(row);
	}
	return triangle;
}

void Q1a() {
	vector<vector<int>> pascalTriangle7 = Generate_Pascal_Triangle(7);
	// Print the pascalTriangle up to n = 7
	for (int n = 0; n < pascalTriangle7.size(); n++)
	{
		cout << n + 1 << ": ";
		for (int k = 0; k < pascalTriangle7[n].size(); k++)
		{
			cout << pascalTriangle7[n][k] << " ";
		}
		cout << endl;
	}
}

void Q1b() {
	vector<vector<int>> pascalTriangle = Generate_Pascal_Triangle(32);
	// Print the pascalTriangle generated using Generate_Pascal_Triangle function
	for (int n = 0; n < pascalTriangle.size(); n++)
	{
		cout << n + 1 << ": ";
		for (int k = 0; k < pascalTriangle[n].size(); k++)
		{
			int temp = pascalTriangle[n][k];
			if (temp % 2 == 0) {
				cout << ".";
			}
			else {
				cout << "#";
			}
		}
		cout << endl;
	}
}

void mat_mul(const double* A, const double* B, double* C, int m, int n, int p) {
	// Initialize C to zero
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			C[i * p + j] = 0;
		}
	}
	// Perform matrix multiplication
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			for (int k = 0; k < n; k++) {
				C[i * p + j] += A[i * n + k] * B[k * p + j];
			}
		}
	}
}

void Q2a() {
	const double A[6] = { 0, 1, 2, 0, 1, 2 }; // 2x3 matrix
	const double B[12] = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 }; // 3x4 matrix
	double C[8]; // Resulting 2x4 matrix
	mat_mul(A, B, C, 2, 3, 4);

	// Print the resulting matrix C
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 4; j++) {
			cout << C[i * 4 + j] << " ";
		}
		cout << endl;
	}

}

void ConwaysGameOfLife(int rows, int cols, vector<vector<int>>& grid) {
	// Create a copy of the grid to store the next state
	vector<vector<int>> nextGrid = grid;
	// Define the directions for the 8 neighbors
	int directions[8][2] = { {-1, -1}, {-1, 0}, {-1, 1},
							 {0, -1},          {0, 1},
							 {1, -1}, {1, 0}, {1, 1} };
	// Iterate through each cell in the grid
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			int liveNeighbors = 0;
			// Count live neighbors
			for (auto& dir : directions) {
				int newRow = r + dir[0];
				int newCol = c + dir[1];
				if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols) {
					liveNeighbors += grid[newRow][newCol];
				}
			}
			// Apply the rules of Conway's Game of Life
			if (grid[r][c] == 1) { // Cell is currently alive
				if (liveNeighbors < 2 || liveNeighbors > 3) {
					nextGrid[r][c] = 0; // Cell dies
				}
			}
			else { // Cell is currently dead
				if (liveNeighbors == 3) {
					nextGrid[r][c] = 1; // Cell becomes alive
				}
			}
		}
	}
	// Update the original grid to the next state
	grid = nextGrid;
}

void Q3a() {
	/* 
	Initialize a grid of size 80x40 with the following pattern:
	0,0,0,1,0,0,0
	0,0,1,0,1,0,0
	0,1,1,0,1,1,0
	Centered at the bottom of the grid.
	I will be printing the grid to the console after 1, 2, 4, 50, 100, and 200 iterations.
	*/
	int rows = 40;
	int cols = 80;
	vector<vector<int>> grid(rows, vector<int>(cols, 0));
	// Initialize the glider pattern at the bottom center
	grid[rows - 3][cols / 2] = 1;
	grid[rows - 2][cols / 2 - 1] = 1;
	grid[rows - 2][cols / 2 + 1] = 1;
	grid[rows - 1][cols / 2 - 1] = 1;
	grid[rows - 1][cols / 2 - 2] = 1;
	grid[rows - 1][cols / 2 + 1] = 1;
	grid[rows - 1][cols / 2 + 2] = 1;
	// Print the initial state
	cout << "Initial State:" << endl;
	for (const auto& row : grid) {
		for (const auto& cell : row) {
			cout << (cell ? '#' : '.');
		}
		cout << endl;
	}
	// Define the number of iterations to print
	vector<int> iterationsToPrint = { 1, 2, 4, 50, 100, 200 };
	int currentIteration = 0;
	for (int i = 1; i <= 200; i++) {
		ConwaysGameOfLife(rows, cols, grid);
		if (i == iterationsToPrint[currentIteration]) {
			cout << "After " << i << " iterations:" << endl;
			for (const auto& row : grid) {
				for (const auto& cell : row) {
					cout << (cell ? '#' : '.');
				}
				cout << endl;
			}
			currentIteration++;
			if (currentIteration >= iterationsToPrint.size()) break; // Stop if all specified iterations are printed
		}
	}
}


int main()
{
	vector<vector<int>> pascalTriangle = Generate_Pascal_Triangle(32);
	vector<vector<int>> pascalTriangle7 = Generate_Pascal_Triangle(7);

	// Below is the inefficient method using the PascalTriangle_Coefficients function

	// Fill the 2d vector using PascalTriangle_Coefficients function
	//for (int n = 0; n <= 31; n++)
	//{
	//	cout << n + 1 << ": ";
	//	for (int k = 0; k <= n; k++)
	//	{
	//		C[n][k] = PascalTriangle_Coefficients(n, k);
	//		// Print the coefficients
	//		if (C[n][k] != 0) {
	//			int temp = C[n][k];

	//			if (temp % 2 == 0) {
	//				cout << ".";
	//			}
	//			else {
	//				cout << "#";
	//			}
	//		}

	//		else
	//			cout << "  "; // Print space for 0 entries
	//	}
	//	cout << endl;
	//}


	// Uncomment for answer to Q1a
	// Q1a();

	// Uncomment for answer to Q1b
	// Q1b();


	/*
	\begin{ enumerate } [a)]
		\item Write a function that takes in two arrays and preforms matrix multiplication. \\
		Your function's signature should be:
		\begin{ verbatim }
	void mat_mul(const double* A, const double* B, double* C, int m, int n, int p);
	\end{ verbatim }
		The function should compute $C = AB$, where $A \in \rr^ { m\times n }, B \in \rr^ { n \times p }$, and $C \in \rr^ { m \times p }$.Test your program on the matrices
		\[
			A = \begin{ pmatrix }
				0 & 1 & 2 \\
				0 & 1 & 2
				\end{ pmatrix },
				\quad
				B = \begin{ pmatrix }
				0 & 1 & 2 & 3 \\
				0 & 1 & 2 & 3 \\
				0 & 1 & 2 & 3
				\end{ pmatrix },
				\]
		which can be generated using the provided \verb | gradient_matrix | function.
		\item \textbf{ Optional. } Measure the time $T(m)$ to generate two random matrices $A, B \in \rr^ { m\times m }$ and multiply them together using your routine.Calculate $T(m)$ for $m = 100, 200, 300, \dots, 1600$.Use linear regression on fit the data to
		\[
			T(m) = Cm ^ \alpha
				\]
		for constants $C$ and $\alpha$, comment on whether the value $\alpha$ is consistent with your matrix multiplication algorithm.
			\end{ enumerate }

	*/

	// Uncomment below for answer to Q2a
	// Q2a();

	// Uncomment below for answer to Q2b (WARNING: This has not been implemented yet)
	// Q2b();	// Optional. Come back to this if time permits.


	// Uncomment below for answer to Q3a
	Q3a();





	return 0;
}

