#ifndef BGENERAL_H_
#define BGENERAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <list>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <assert.h>

#ifndef TARGET_TYPE
#define TARGET_TYPE double
#endif

using namespace std;

#define sq(x)		((x)*(x))

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

// Check versus Velocity Verlet in cart coords
enum{
  SERV_QX, SERV_QY, SERV_QZ,
  SERV_VX, SERV_VY, SERV_VZ,
  SERV_AX, SERV_AY, SERV_AZ,
  VASS_QX, VASS_QY, VASS_QZ,
  VASS_VX, VASS_VY, VASS_VZ,
  VASS_AX, VASS_AY, VASS_AZ,
  SAYS_QX, SAYS_QY, SAYS_QZ,
  SAYS_VX, SAYS_VY, SAYS_VZ,
  SAYS_AX, SAYS_AY, SAYS_AZ
};


/**************************************
 * 		General Functions             *
 **************************************/

inline int RandomIntRange(int min, int max)
{
  assert(max > min);
  return rand()%(max-min+1) + min;
}

inline float RandomRealRange(float min, float max)
{
  assert(max > min);
  return min + (((float)rand()) / (float)RAND_MAX) * (max-min);
}

inline double RandomRealRange(double min, double max)
{
  assert(max > min);
  return min + (((double)rand()) / (double)RAND_MAX) * (max-min);
}

inline long double RandomRealRange(long double min, long double max)
{
  assert(max > min);
  return min + (((long double)rand()) / (long double)RAND_MAX) * (max-min);
}


double round(double r);

/*
 *  Aminoacids 3 letter notation to 1 letter notation
 */
char aa321 (const char *aa);

/*
 * aa321 inverse
 * ?????????????????????????
 * ?????????????????????????
 */
char aa123 (char *dest, char aa);
char aa123 (std::string& dest, char aa);

/*
 *  If not, puts '/' at the end of RESULT
 */
void bCheck_path (char *result, const char *path);

/*
 * Puts START-STOP SRC substring into DEST
 * Character is counted from 0
 * Pointers have to be allocated
 */
int bExtractStr (char *dest, const char *src, int start, int stop);

/*
 * Tolower bExtractStr
 */
int bExtractTolowerStr ( char *dest, const char *src, int start, int stop );

/*
 * Modifies all characters in a string to NULL
 * Returns how many characters has modified
 */
int bZeroStr(char *dest);

/*
 * Modifies NO_ELEM elements in an array to NULL
 * Returns how many elements has modified
 */int bZeroCharArray(char *array, int no_elem);


/*
 * Awk Substr
 */
int bSubstr (char *dest, const char *src, int start, int no_chars);

/*
 * Left trim
 */
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

/*
 * Right trim
 */
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

/*
 * Trim from both ends
 */
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

/*
 * Decimal prefix of zeros to limit
 */
string decimal_prefix(double inp_no, long int limit);


#endif /*BGENERAL_H_*/
