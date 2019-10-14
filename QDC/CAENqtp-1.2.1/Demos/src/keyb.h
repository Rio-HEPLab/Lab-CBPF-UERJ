#ifdef LINUX

/*****************************************************************************/
/*  GETCH  */
/*****************************************************************************/
int getch(void);


/*****************************************************************************/
/*  KBHIT  */
/*****************************************************************************/
int kbhit();

#else
#include <conio.h>
#include <ctype.h>

// getch and kbhit was declared deprecated in Microsoft Visual C++. The ISO C++
// conformant name with '_' preceding the name is suggested.
#define getch _getch
#define kbhit _kbhit

#endif

/*****************************************************************************/
/*  GETSTRING  */
/*****************************************************************************/
int getstring( char* string, int max_len);

