#ifndef GUARD_List
#define GUARD_List

// Note: This is simply an implementation of a doubly-linked list.
// Each element of the doubly-linked list is called a Cell.
// The List itself contains pointers to the initial and final Cells.

class Cell {
public:
	int Beg; //flag for is it the first cell(1=yes 0=no)
	int End; //flag for is it the last cell (1=yes 0=no)
	Cell *before;
	Cell *after;
	void *contents;
};

class List {
public:
	Cell *first;
	Cell *last;
	int length;
	void deall();
	List();
	~List();
};

typedef void *voidP;

void CombineLists(List *list1,List *list2, List **list);
void PrintList(List *list);
void DelCell(List *list,Cell *cell);
void AddCellToEnd(List *list, void *p);
void AddCellToBeg(List *list, void *p);
void AddCellAfter(List *list,Cell *oldcell,void *p);
void AddCellBefore(List *list,Cell *oldcell,void *p);
void ListToVector(List *list,voidP **p,int *n);

#endif



