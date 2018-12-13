#include <iostream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <math.h>
#include <stack>
#include <float.h>
#include <set>
#include <cstring>
#include <iomanip>
using namespace std;

#define INPUT_UNKNOWN_VALUE 1
#define UNKNOWN_VALUE -1
#define NOT_A_CELL 0

// total number of cells
int totalNumberOfCells = 0;
// Dimensions of the finalBoard
int width, height;
// Main board
vector<vector<int> >board;
// Creating a shortcut for int, int pair type 
typedef pair<int, int> Pair;
// Final Path
vector<Pair>finalPath;
// Start and End points
int startI, startJ, finalI, finalJ;

void readInput() {
	//Input the user input file 'input.txt'
	ifstream infile("input.txt");

	infile >> height >> width;
	for (int i = 0; i < height; ++i) {
		board.push_back(vector<int>());
		for (int j = 0; j < width; ++j) {
			int temp;
			infile >> temp;
			if (temp == INPUT_UNKNOWN_VALUE) {
				totalNumberOfCells++;
				board[i].push_back(UNKNOWN_VALUE);
			}
			else {
				board[i].push_back(temp);
			}
		}
	}
}

void initStartAndEndCells() {
	srand(time(NULL));

	// Generate valid coordinates for the first number
	while (true) {
		startI = rand() % height;
		startJ = rand() % width;
		if (board[startI][startJ] == UNKNOWN_VALUE) {
			break;
		}
	}

	// Generate valid coordinates for the final number
	while (true) {
		finalI = rand() % height;
		finalJ = rand() % width;
		if (!(finalI == startI && finalJ == startJ) && board[finalI][finalJ] == UNKNOWN_VALUE) {
			break;
		}
	}
	board[startI][startJ] = 1;
	board[finalI][finalJ] = totalNumberOfCells;
	cout << "START: " << "(" << startI + 1 << ", " << startJ + 1 << ")" << endl;
	cout << "END: " << "(" << finalI + 1 << ", " << finalJ + 1 << ")" << endl;
}

void outputBoard() {
	ofstream outfile("output.txt");
	outfile << height << " " << width << endl;
	cout << "Board(generate):" << endl;
	cout << "  ";
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			if (board[i][j] / 10 > 0 || board[i][j] < 0) {
				cout << board[i][j] << setw(4);
				outfile << board[i][j] << " ";
			}
			else {
				cout << board[i][j] << setw(4);
				outfile << board[i][j] << " ";
			}
		}
		cout << endl;
		outfile << endl;
	}
}

bool boardFilled() {
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			if (board[i][j] == UNKNOWN_VALUE) {
				return false;
			}
		}
	}
	return true;
}

// A Utility Function to check whether given cell (row, col) 
// is a valid cell or not. 
bool isValid(int row, int col) {
	// Returns true if row number and column number 
	// is in range 
	return (row >= 0) && (row < height) &&
		(col >= 0) && (col < width);
}

/**
 * Returns coordinates of neighbour cell(for 2 included cells) that is empty
 * if there's no empty neighbour cells, returns (-1, -1)
 **/
Pair twoCellsEmptyNeighbourIndex(Pair firstCell, Pair secondCell) {
	if (abs(firstCell.first - secondCell.first) > 1 || abs(firstCell.second - secondCell.second) > 1) {
		// cout<<"Two cells: ("<<firstCell.second<<":"<<firstCell.first<<") and ("<<secondCell.second<<":"<<secondCell.first<<") are not connected!"<<endl;
		return make_pair(-1, -1);
	}
	for (int i = -1; i < 2; i++) {
		for (int j = -1; j < 2; j++) {
			if (i == 0 && j == 0) {
				continue;
			}
			if (!isValid(firstCell.first + i, firstCell.second + j)) {
				continue;
			}
			for (int ii = -1; ii < 2; ii++) {
				for (int jj = -1; jj < 2; jj++) {
					if (ii == 0 && jj == 0) {
						continue;
					}
					if (!isValid(secondCell.first + ii, secondCell.second + jj)) {
						continue;
					}
					if (firstCell.first + i == secondCell.first + ii && firstCell.second + j == secondCell.second + jj) {
						if (board[firstCell.first + i][firstCell.second + j] == UNKNOWN_VALUE) {
							// cout<<"neighbour for ("<<firstCell.second<<":"<<firstCell.first<<") and ("<<secondCell.second<<":"<<secondCell.first<<") is ("<<firstCell.second+j<<":"<<firstCell.first+i<<")"<<endl;
							return make_pair(firstCell.first + i, firstCell.second + j);
						}
					}
				}
			}
		}
	}
	return make_pair(-1, -1);
}

Pair emptyNeighbourIndex(Pair current) {
	for (int i = -1; i < 2; i++) {
		for (int j = -1; j < 2; j++) {
			if (i == 0 && j == 0) {
				continue;
			}
			if (!isValid(current.first + i, current.second + j)) {
				continue;
			}
			if (board[current.first + i][current.second + j] == UNKNOWN_VALUE) {
				return make_pair(current.first + i, current.second + j);
			}
		}
	}
	return make_pair(-1, -1);
}

bool findFinalPath() {
	// Fill traced cells by 1 temporary
	for (int i = 0; i < finalPath.size(); ++i) {
		board[finalPath[i].first][finalPath[i].second] = 1;
	}

	bool success;
	while (!boardFilled()) {
		success = false;
		for (int i = 0; i < finalPath.size() - 1; ++i) {
			Pair emptyNeighbour = twoCellsEmptyNeighbourIndex(finalPath[i], finalPath[i + 1]);
			if (emptyNeighbour.first != -1 && emptyNeighbour.second != -1) {
				finalPath.insert(finalPath.begin() + i + 1, emptyNeighbour);
				board[emptyNeighbour.first][emptyNeighbour.second] = 1;
				success = true;
				break;
			}
		}
		if (!success) {
			// cout<<"FAILED INTERNAL"<<endl;
			break;
		}
	}

	if (success) {
		// cout<<"Total Success"<<endl;
		return true;
	}
	// If there's incompleted cells on the edge
	while (!boardFilled()) {
		success = true;
		Pair startNeighbourIndex = emptyNeighbourIndex(finalPath[0]);
		if (startNeighbourIndex.first != -1 && startNeighbourIndex.second != -1) {
			board[startNeighbourIndex.first][startNeighbourIndex.second] = 1;
			finalPath.insert(finalPath.begin(), startNeighbourIndex);
			continue;
		}
		Pair endNeighbourIndex = emptyNeighbourIndex(finalPath[finalPath.size() - 1]);
		if (endNeighbourIndex.first != -1 && endNeighbourIndex.second != -1) {
			board[endNeighbourIndex.first][endNeighbourIndex.second] = 1;
			finalPath.insert(finalPath.end(), endNeighbourIndex);
		}
		else {
			success = false;
			break;
		}
	}
	if (success) {
		// cout<<"Total Success"<<endl;
		return true;
	}
	else {
		// cout<<"Total Fail"<<endl;
		return false;
	}
}

void cutHints() {
	if (totalNumberOfCells % 2 != 0) {
		for (int i = 1; i < finalPath.size() - 1; i += 2) {
			board[finalPath[i].first][finalPath[i].second] = -1;
		}
	}
	else {
		for (int i = 1; i < finalPath.size() - 2; i += 2) {
			board[finalPath[i].first][finalPath[i].second] = -1;
		}
		board[finalPath[finalPath.size() - 2].first][finalPath[finalPath.size() - 2].second] = -1;
	}
}

void reset() {
	finalPath.clear();
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			if (board[i][j] != NOT_A_CELL) {
				board[i][j] = UNKNOWN_VALUE;
			}
		}
	}
}

void fillBoard() {
	int counter = 1;
	for (int i = 0; i < finalPath.size(); ++i) {
		board[finalPath[i].first][finalPath[i].second] = counter++;
	}
}

/**
 * Code block for implementation A* alghorithm
 **/

 // Creating a shortcut for pair<int, pair<int, int>> type 
typedef pair<double, pair<int, int> > pPair;
// A structure to hold the neccesary parameters 
struct cell {
	// Row and Column index of its parent 
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1 
	int parent_i, parent_j;
	// f = g + h 
	double f, g, h;
};

// A Utility Function to check whether the given cell is 
// blocked or not 
bool isUnBlocked(int row, int col) {
	// Returns true if the cell is not blocked else false 
	return board[row][col] != NOT_A_CELL;
}

// A Utility Function to check whether destination cell has 
// been reached or not 
bool isDestination(int row, int col, Pair dest) {
	return row == dest.first && col == dest.second;
}

// A Utility Function to calculate the 'h' heuristics. 
double calculateHValue(int row, int col, Pair dest) {
	// Return using the distance formula 
	return ((double)sqrt((row - dest.first)*(row - dest.first)
		+ (col - dest.second)*(col - dest.second)));
}

// A Utility Function to trace the path from the source 
// to destination 
void tracePath(vector<vector<cell> >cellDetails, Pair dest) {
	// printf ("\nThe Init Path is "); 
	int row = dest.first;
	int col = dest.second;

	stack<Pair> Path;

	while (!(cellDetails[row][col].parent_i == row
		&& cellDetails[row][col].parent_j == col)) {
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row][col].parent_i;
		int temp_col = cellDetails[row][col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	while (!Path.empty()) {
		pair<int, int> p = Path.top();
		finalPath.push_back(p);
		Path.pop();
		// printf("-> (%d,%d) ",p.second,p.first); 
	}
	// cout<<endl;

	return;
}

// A Function to find the shortest path between 
// a given source cell to a destination cell according 
// to A* Search Algorithm 
void connectStartToEnd() {
	Pair src = make_pair(startI, startJ);
	Pair dest = make_pair(finalI, finalJ);
	// If the source is out of range 
	if (!isValid(src.first, src.second)) {
		// printf ("Source is invalid\n"); 
		return;
	}
	// If the destination is out of range 
	if (!isValid(dest.first, dest.second)) {
		// printf ("Destination is invalid\n"); 
		return;
	}
	// Either the source or the destination is blocked 
	if (!isUnBlocked(src.first, src.second) ||
		!isUnBlocked(dest.first, dest.second)) {
		// printf ("Source or the destination is blocked\n"); 
		return;
	}
	// If the destination cell is the same as source cell 
	if (isDestination(src.first, src.second, dest)) {
		// printf ("We are already at the destination\n"); 
		return;
	}

	// Create a closed list and initialise it to false which means 
	// that no cell has been included yet 
	// This closed list is implemented as a boolean 2D array 
	bool closedList[height][width];
	memset(closedList, false, sizeof(closedList));
	// Declare a 2D array of structure to hold the details 
	//of that cell 
	vector<vector<cell> >cellDetails;
	for (int i = 0; i < height; ++i) {
		cellDetails.push_back(vector<cell>());
		for (int j = 0; j < width; ++j) {
			cell c;
			cellDetails[i].push_back(c);
		}
	}

	int i, j;

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			cellDetails[i][j].f = FLT_MAX;
			cellDetails[i][j].g = FLT_MAX;
			cellDetails[i][j].h = FLT_MAX;
			cellDetails[i][j].parent_i = -1;
			cellDetails[i][j].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node 
	i = src.first, j = src.second;
	cellDetails[i][j].f = 0.0;
	cellDetails[i][j].g = 0.0;
	cellDetails[i][j].h = 0.0;
	cellDetails[i][j].parent_i = i;
	cellDetails[i][j].parent_j = j;

	/*
	 Create an open list having information as-
	 <f, <i, j>>
	 where f = g + h,
	 and i, j are the row and column index of that cell
	 Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	 This open list is implenented as a set of pair of pair.*/
	set<pPair> openList;

	// Put the starting cell on the open list and set its 
	// 'f' as 0 
	openList.insert(make_pair(0.0, make_pair(i, j)));

	// We set this boolean value as false as initially 
	// the destination is not reached. 
	bool foundDest = false;

	while (!openList.empty())
	{
		pPair p = *openList.begin();

		// Remove this vertex from the open list 
		openList.erase(openList.begin());

		// Add this vertex to the closed list 
		i = p.second.first;
		j = p.second.second;
		closedList[i][j] = true;

		/*
		 Generating all the 8 successor of this cell

			 N.W   N   N.E
			   \   |   /
				\  |  /
			 W----Cell----E
				  / | \
				/   |  \
			 S.W    S   S.E

		 Cell-->Popped Cell (i, j)
		 N -->  North       (i-1, j)
		 S -->  South       (i+1, j)
		 E -->  East        (i, j+1)
		 W -->  West           (i, j-1)
		 N.E--> North-East  (i-1, j+1)
		 N.W--> North-West  (i-1, j-1)
		 S.E--> South-East  (i+1, j+1)
		 S.W--> South-West  (i+1, j-1)*/

		 // To store the 'g', 'h' and 'f' of the 8 successors 
		double gNew, hNew, fNew;

		//----------- 1st Successor (North) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j].parent_i = i;
				cellDetails[i - 1][j].parent_j = j;
				// printf ("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}
			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j] == false &&
				isUnBlocked(i - 1, j) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i - 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j].f == FLT_MAX ||
					cellDetails[i - 1][j].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i - 1, j)));

					// Update the details of this cell 
					cellDetails[i - 1][j].f = fNew;
					cellDetails[i - 1][j].g = gNew;
					cellDetails[i - 1][j].h = hNew;
					cellDetails[i - 1][j].parent_i = i;
					cellDetails[i - 1][j].parent_j = j;
				}
			}
		}

		//----------- 2nd Successor (South) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j].parent_i = i;
				cellDetails[i + 1][j].parent_j = j;
				// printf("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}
			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j] == false &&
				isUnBlocked(i + 1, j) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i + 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j].f == FLT_MAX ||
					cellDetails[i + 1][j].f > fNew)
				{
					openList.insert(make_pair(fNew, make_pair(i + 1, j)));
					// Update the details of this cell 
					cellDetails[i + 1][j].f = fNew;
					cellDetails[i + 1][j].g = gNew;
					cellDetails[i + 1][j].h = hNew;
					cellDetails[i + 1][j].parent_i = i;
					cellDetails[i + 1][j].parent_j = j;
				}
			}
		}

		//----------- 3rd Successor (East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i, j + 1) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i][j + 1].parent_i = i;
				cellDetails[i][j + 1].parent_j = j;
				// printf("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i][j + 1] == false &&
				isUnBlocked(i, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i][j + 1].f == FLT_MAX ||
					cellDetails[i][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i, j + 1)));

					// Update the details of this cell 
					cellDetails[i][j + 1].f = fNew;
					cellDetails[i][j + 1].g = gNew;
					cellDetails[i][j + 1].h = hNew;
					cellDetails[i][j + 1].parent_i = i;
					cellDetails[i][j + 1].parent_j = j;
				}
			}
		}

		//----------- 4th Successor (West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i, j - 1) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i][j - 1].parent_i = i;
				cellDetails[i][j - 1].parent_j = j;
				// printf("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i][j - 1] == false &&
				isUnBlocked(i, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i][j - 1].f == FLT_MAX ||
					cellDetails[i][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i, j - 1)));

					// Update the details of this cell 
					cellDetails[i][j - 1].f = fNew;
					cellDetails[i][j - 1].g = gNew;
					cellDetails[i][j - 1].h = hNew;
					cellDetails[i][j - 1].parent_i = i;
					cellDetails[i][j - 1].parent_j = j;
				}
			}
		}

		//----------- 5th Successor (North-East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j + 1) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j + 1].parent_i = i;
				cellDetails[i - 1][j + 1].parent_j = j;
				// printf ("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j + 1] == false &&
				isUnBlocked(i - 1, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i - 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j + 1].f == FLT_MAX ||
					cellDetails[i - 1][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i - 1, j + 1)));

					// Update the details of this cell 
					cellDetails[i - 1][j + 1].f = fNew;
					cellDetails[i - 1][j + 1].g = gNew;
					cellDetails[i - 1][j + 1].h = hNew;
					cellDetails[i - 1][j + 1].parent_i = i;
					cellDetails[i - 1][j + 1].parent_j = j;
				}
			}
		}

		//----------- 6th Successor (North-West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i - 1, j - 1) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i - 1, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i - 1][j - 1].parent_i = i;
				cellDetails[i - 1][j - 1].parent_j = j;
				// printf ("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i - 1][j - 1] == false &&
				isUnBlocked(i - 1, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i - 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i - 1][j - 1].f == FLT_MAX ||
					cellDetails[i - 1][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew, make_pair(i - 1, j - 1)));
					// Update the details of this cell 
					cellDetails[i - 1][j - 1].f = fNew;
					cellDetails[i - 1][j - 1].g = gNew;
					cellDetails[i - 1][j - 1].h = hNew;
					cellDetails[i - 1][j - 1].parent_i = i;
					cellDetails[i - 1][j - 1].parent_j = j;
				}
			}
		}

		//----------- 7th Successor (South-East) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j + 1) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j + 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j + 1].parent_i = i;
				cellDetails[i + 1][j + 1].parent_j = j;
				// printf ("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j + 1] == false &&
				isUnBlocked(i + 1, j + 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i + 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j + 1].f == FLT_MAX ||
					cellDetails[i + 1][j + 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i + 1, j + 1)));

					// Update the details of this cell 
					cellDetails[i + 1][j + 1].f = fNew;
					cellDetails[i + 1][j + 1].g = gNew;
					cellDetails[i + 1][j + 1].h = hNew;
					cellDetails[i + 1][j + 1].parent_i = i;
					cellDetails[i + 1][j + 1].parent_j = j;
				}
			}
		}

		//----------- 8th Successor (South-West) ------------ 

		// Only process this cell if this is a valid one 
		if (isValid(i + 1, j - 1) == true)
		{
			// If the destination cell is the same as the 
			// current successor 
			if (isDestination(i + 1, j - 1, dest) == true)
			{
				// Set the Parent of the destination cell 
				cellDetails[i + 1][j - 1].parent_i = i;
				cellDetails[i + 1][j - 1].parent_j = j;
				// printf("The destination cell is found\n"); 
				tracePath(cellDetails, dest);
				foundDest = true;
				return;
			}

			// If the successor is already on the closed 
			// list or if it is blocked, then ignore it. 
			// Else do the following 
			else if (closedList[i + 1][j - 1] == false &&
				isUnBlocked(i + 1, j - 1) == true)
			{
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i + 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to 
				// the open list. Make the current square 
				// the parent of this square. Record the 
				// f, g, and h costs of the square cell 
				//                OR 
				// If it is on the open list already, check 
				// to see if this path to that square is better, 
				// using 'f' cost as the measure. 
				if (cellDetails[i + 1][j - 1].f == FLT_MAX ||
					cellDetails[i + 1][j - 1].f > fNew)
				{
					openList.insert(make_pair(fNew,
						make_pair(i + 1, j - 1)));

					// Update the details of this cell 
					cellDetails[i + 1][j - 1].f = fNew;
					cellDetails[i + 1][j - 1].g = gNew;
					cellDetails[i + 1][j - 1].h = hNew;
					cellDetails[i + 1][j - 1].parent_i = i;
					cellDetails[i + 1][j - 1].parent_j = j;
				}
			}
		}
	}

	// When the destination cell is not found and the open 
	// list is empty, then we conclude that we failed to 
	// reach the destiantion cell. This may happen when the 
	// there is no way to destination cell (due to blockages) 
	if (foundDest == false)
		// printf("Failed to find the Destination Cell\n"); 

		return;
}

int main() {
	int attempt = 1;
	readInput();
	while (true) {
		initStartAndEndCells();
		connectStartToEnd();
		if (findFinalPath()) {
			fillBoard();
			//cout<<"Attempt #"<<attempt<<" successful"<<endl;
			cutHints();
			outputBoard();
			break;
		}
		else {
			cout << "Attempt #" << attempt++ << " failed" << endl;
			reset();
		}
	}
}
