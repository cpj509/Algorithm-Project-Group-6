#include <iostream>
#include <vector>
#include <algorithm> 
#include <fstream>
#include <iomanip>

using namespace std;

#define UNKNOWN_VALUE -1
#define NOT_A_CELL 0

vector<int>given;
int width, height;
int startI, startJ;
vector<vector<int> >board;

void displayFinalBoard() {
	ofstream outfile("solve.txt");
	cout << "Board(solve):" << endl;
	outfile << width << " " << height << endl;
	cout << "  ";
	for (int i = 1; i < width + 1; ++i) {
		for (int j = 1; j < height + 1; ++j) {
			cout << board[i][j] << setw(4);
			outfile << board[i][j] << " ";
		}
		cout << endl;
		outfile << endl;
	}
}

void setup() {
	ifstream infile("output.txt");
	infile >> width >> height;
	for (int i = 0; i < width + 2; ++i) {
		board.push_back(vector<int>());
		for (int j = 0; j < height + 2; ++j) {
			board[i].push_back(NOT_A_CELL);
		}
	}

	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			int temp;
			infile >> temp;
			if (temp != 0 && temp != -1) {
				given.push_back(temp);
				if (temp == 1) {
					startI = i + 1;
					startJ = j + 1;
				}
			}
			board[i + 1][j + 1] = temp;
		}
	}

	sort(given.begin(), given.end());
}

bool solve(int currentI, int currentJ, int n, int next) {
	if (n > given[given.size() - 1]) {
		return true;
	}
	if (board[currentI][currentJ] != UNKNOWN_VALUE && board[currentI][currentJ] != n) {
		return false;
	}
	if (board[currentI][currentJ] == UNKNOWN_VALUE && given[next] == n) {
		return false;
	}

	int back = UNKNOWN_VALUE;
	if (board[currentI][currentJ] == n) {
		next++;
		back = n;
	}

	board[currentI][currentJ] = n;
	for (int i = -1; i < 2; ++i) {
		for (int j = -1; j < 2; ++j) {
			if (solve(currentI + i, currentJ + j, n + 1, next)) {
				return true;
			}
		}
	}

	board[currentI][currentJ] = back;
	return false;
}

int main() {
	setup();
	solve(startI, startJ, 1, 0);
	displayFinalBoard();

	return 0;
}
