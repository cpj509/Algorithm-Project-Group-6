#include <iostream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

#define UNKNOWN_VALUE -1
#define NOT_A_CELL 0

// total number of cells
int totalNumberOfCells;
// List of the numbers that appears on the finalBoard (for the generation it's 1 and final number)
vector<int>given;
// Dimensions of the finalBoard
int width, height;

// Board for the solving(filling)
vector<vector<int> >blankBoard;
// Board for the generation
vector<vector<int> >finalBoard;

// This funciton initializes blankBoard with the user input
// and returns number of cells for the calculation of the final number
int setup() {
	ifstream infile("input.txt");
	totalNumberOfCells = 0;
	infile >> width >> height;
	for (int i = 0; i < width + 2; ++i) {
		blankBoard.push_back(vector<int>());
		for (int j = 0; j < height + 2; ++j) {
			blankBoard[i].push_back(NOT_A_CELL);
		}
	}
	for (int i = 0; i < width; ++i) {
		finalBoard.push_back(vector<int>());
		for (int j = 0; j < height; ++j) {
			finalBoard[i].push_back(UNKNOWN_VALUE);
		}
	}

	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			int temp;
			infile >> temp;
			if (temp != NOT_A_CELL) {
				totalNumberOfCells++;
				// Change 1 to UNKNOWN_VALUE(-1) in order to escape error
				if (temp == 1) {
					temp = UNKNOWN_VALUE;
				}
			}
			else {
				finalBoard[i][j] = NOT_A_CELL;
			}
			blankBoard[i + 1][j + 1] = temp;
		}
	}
	return totalNumberOfCells;
}

// This funcion displays final finalBoard
void displayFinalBoard() {
	cout << endl;
	ofstream outfile("output.txt");
	cout << "GENERATED BOARD:" << endl;
	outfile << width << " " << height << endl;
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			cout << finalBoard[i][j] << " ";
			outfile << finalBoard[i][j] << " ";
		}
		cout << endl;
		outfile << endl;
	}

}

bool solve(int currentI, int currentJ, int n, int next) {
	if (n > given[given.size() - 1]) {
		return true;
	}
	if (blankBoard[currentI][currentJ] != UNKNOWN_VALUE && blankBoard[currentI][currentJ] != n) {
		return false;
	}
	if (blankBoard[currentI][currentJ] == UNKNOWN_VALUE && given[next] == n) {
		return false;
	}

	int back = UNKNOWN_VALUE;
	if (blankBoard[currentI][currentJ] == n) {
		next++;
		back = n;
	}

	blankBoard[currentI][currentJ] = n;
	for (int i = -1; i < 2; ++i) {
		for (int j = -1; j < 2; ++j) {
			if (solve(currentI + i, currentJ + j, n + 1, next)) {
				return true;
			}
		}
	}

	blankBoard[currentI][currentJ] = back;
	return false;
}

// This method inserts 1 and final number, 
// also chooses coordinates for them.
// Then calls solving algo
void prepareBoard() {
	// Insert 1 and final number to 'given'
	given.push_back(1);
	given.push_back(totalNumberOfCells);

	srand(time(NULL));
	int startI, startJ, finalI, finalJ;

	// Generate valid coordinates for the first number
	while (true) {
		startI = rand() % width + 1;
		startJ = rand() % height + 1;
		if (blankBoard[startI][startJ] != NOT_A_CELL) {
			break;
		}
	}

	// Generate valid coordinates for the final number
	while (true) {
		finalI = rand() % width + 1;
		finalJ = rand() % height + 1;
		if (!(finalI == startI && finalJ == startJ) && blankBoard[finalI][finalJ] != NOT_A_CELL) {
			break;
		}
	}

	cout << 1 << " at (" << startI << ":" << startJ << ")" << endl;
	cout << totalNumberOfCells << " at (" << finalI << ":" << finalJ << ")" << endl;

	// Insert first and final numbers to the blankBoard
	blankBoard[startI][startJ] = 1;
	blankBoard[finalI][finalJ] = totalNumberOfCells;

	// Call solving algo
	solve(startI, startJ, 1, 0);
}

// This function builds the final board
// by copying only 1/2 of all numbers to the finalBoard
void buildFinalBoard() {
	bool copyFlag = true;
	for (int t = 1; t <= totalNumberOfCells; ++t) {
		if (copyFlag) {
			for (int i = 1; i < width + 2; ++i) {
				bool done = false;
				for (int j = 1; j < height + 2; ++j) {
					if (blankBoard[i][j] == t) {
						finalBoard[i - 1][j - 1] = t;
						done = true;
						break;
					}
				}
				if (done) {
					break;
				}
			}
		}
		copyFlag = !copyFlag;
	}
}

int main() {
	setup();
	prepareBoard();
	buildFinalBoard();
	displayFinalBoard();
	return 0;
}