/* A Bison parser, made by GNU Bison 3.3.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2019 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.3.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         gmsh_yyparse
#define yylex           gmsh_yylex
#define yyerror         gmsh_yyerror
#define yydebug         gmsh_yydebug
#define yynerrs         gmsh_yynerrs

#define yylval          gmsh_yylval
#define yychar          gmsh_yychar

/* First part of user prologue.  */
#line 1 "Gmsh.y" /* yacc.c:337  */

// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <sstream>
#include <map>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "GmshConfig.h"
#include "GmshMessage.h"
#include "fullMatrix.h"
#include "MallocUtils.h"
#include "ListUtils.h"
#include "TreeUtils.h"
#include "StringUtils.h"
#include "Numeric.h"
#include "Context.h"
#include "GModel.h"
#include "GModelIO_GEO.h"
#include "GModelIO_OCC.h"
#include "GeoDefines.h"
#include "ExtrudeParams.h"
#include "Options.h"
#include "Parser.h"
#include "OpenFile.h"
#include "CommandLine.h"
#include "FunctionManager.h"
#include "ColorTable.h"
#include "OS.h"
#include "CreateFile.h"
#include "gmshSurface.h"
#include "gmshLevelset.h"
#include "fullMatrix.h"

#if defined(HAVE_MESH)
#include "Generator.h"
#include "Field.h"
#include "BackgroundMesh.h"
#include "HighOrder.h"
#endif

#if defined(HAVE_POST)
#include "PView.h"
#include "PViewDataList.h"
#endif

#if defined(HAVE_PLUGINS)
#include "PluginManager.h"
#endif

#if defined(HAVE_OPENGL)
#include "drawContext.h"
#endif

#if defined(HAVE_POPPLER)
#include "gmshPopplerWrapper.h"
#endif

#define MAX_RECUR_TESTS 100
#define MAX_RECUR_LOOPS 100

// global parser variables
std::string gmsh_yyname;
int gmsh_yyerrorstate = 0;
int gmsh_yyviewindex = 0;
std::map<std::string, gmsh_yysymbol> gmsh_yysymbols;
std::map<std::string, std::vector<std::string> > gmsh_yystringsymbols;
std::string gmsh_yyfactory;
NameSpaces gmsh_yynamespaces;

// static parser variables (accessible only in this file)
#if defined(HAVE_POST)
static PViewDataList *ViewData = 0;
#endif
static std::vector<double> ViewCoord;
static std::vector<double> *ViewValueList = 0;
static int *ViewNumList = 0;
static ExtrudeParams extr;
static gmshSurface *myGmshSurface = 0;
static int statusImbricatedTests[MAX_RECUR_TESTS];
static int ImbricatedLoop = 0, ImbricatedTest = 0;
static fpos_t yyposImbricatedLoopsTab[MAX_RECUR_LOOPS];
static int yylinenoImbricatedLoopsTab[MAX_RECUR_LOOPS];
static double LoopControlVariablesTab[MAX_RECUR_LOOPS][3];
static std::string LoopControlVariablesNameTab[MAX_RECUR_LOOPS];
static std::string struct_name, struct_namespace;
static int flag_tSTRING_alloc = 0;
static int dim_entity;

static std::map<std::string, std::vector<double> > floatOptions;
static std::map<std::string, std::vector<std::string> > charOptions;
static int flag_Enum, member_ValMax;

void init_options(int member_ValMax_ = 0)
{
  floatOptions.clear(); charOptions.clear();
  flag_Enum = 0; member_ValMax = member_ValMax_;
}

// parser functions defined at the end of this file
void yyerror(const char *s);
void yymsg(int level, const char *fmt, ...);
char *strsave(char *ptr);
void skip(const char *skip, const char *until);
void skipTest(const char *skip, const char *until,
              const char *until2, int l_until2_sub, int *type_until2);
void assignVariable(const std::string &name, int index, int assignType,
                    double value);
void assignVariables(const std::string &name, List_T *indices, int assignType,
                     List_T *values);
void incrementVariable(const std::string &name, int index, double value);
int printListOfDouble(char *format, List_T *list, char *buffer);
fullMatrix<double> ListOfListOfDouble2Matrix(List_T *list);
void ListOfDouble2Vector(List_T *list, std::vector<int> &v);
void ListOfDouble2Vector(List_T *list, std::vector<double> &v);
void ListOfShapes2VectorOfPairs(List_T *list, std::vector<std::pair<int, int> > &v);
void VectorOfPairs2ListOfShapes(const std::vector<std::pair<int, int> > &v, List_T *list);
void addPeriodicEdge(int, int, const std::vector<double>&);
void addPeriodicFace(int, int, const std::map<int, int>&);
void addPeriodicFace(int, int, const std::vector<double>&);
void computeAffineTransformation(SPoint3&, SPoint3&, double, SPoint3&,
                                 std::vector<double>&);
void addEmbedded(int dim, std::vector<int> tags, int dim2, int tag2);
void getAllElementaryTags(int dim, List_T *in);
void getAllPhysicalTags(int dim, List_T *in);
void getElementaryTagsForPhysicalGroups(int dim, List_T *in, List_T *out);
void getElementaryTagsInBoundingBox(int dim, double x1, double y1, double z1,
                                    double x2, double y2, double z2, List_T *out);
void getParentTags(int dim, List_T *in, List_T *out);
void getBoundingBox(int dim, List_T *in, List_T *out);
void setVisibility(int dim, int visible, bool recursive);
void setVisibility(const std::vector<std::pair<int, int> > &dimTags, int visible,
                   bool recursive);
void setColor(const std::vector<std::pair<int, int> > &dimTags, unsigned int val,
              bool recursive);

double treat_Struct_FullName_Float
  (char* c1, char* c2, int type_var = 1, int index = 0,
   double val_default = 0., int type_treat = 0);
double treat_Struct_FullName_dot_tSTRING_Float
  (char* c1, char* c2, char* c3, int index = 0,
   double val_default = 0., int type_treat = 0);
List_T * treat_Struct_FullName_dot_tSTRING_ListOfFloat
  (char* c1, char* c2, char* c3);
int treat_Struct_FullName_dot_tSTRING_Float_getDim
  (char* c1, char* c2, char* c3);
char* treat_Struct_FullName_String
  (char* c1, char* c2, int type_var = 1, int index = 0,
   char* val_default = nullptr, int type_treat = 0);
char* treat_Struct_FullName_dot_tSTRING_String
  (char* c1, char* c2, char* c3, int index = 0,
   char* val_default = nullptr, int type_treat = 0);
List_T * treat_Struct_FullName_dot_tSTRING_ListOfString
  (char* c1, char* c2, char* c3);

struct doubleXstring{
  double d;
  char *s;
};


#line 243 "Gmsh.tab.cpp" /* yacc.c:337  */
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "Gmsh.tab.hpp".  */
#ifndef YY_GMSH_YY_GMSH_TAB_HPP_INCLUDED
# define YY_GMSH_YY_GMSH_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int gmsh_yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    tDOUBLE = 258,
    tSTRING = 259,
    tBIGSTR = 260,
    tEND = 261,
    tAFFECT = 262,
    tDOTS = 263,
    tSCOPE = 264,
    tPi = 265,
    tMPI_Rank = 266,
    tMPI_Size = 267,
    tEuclidian = 268,
    tCoordinates = 269,
    tTestLevel = 270,
    tExp = 271,
    tLog = 272,
    tLog10 = 273,
    tSqrt = 274,
    tSin = 275,
    tAsin = 276,
    tCos = 277,
    tAcos = 278,
    tTan = 279,
    tRand = 280,
    tAtan = 281,
    tAtan2 = 282,
    tSinh = 283,
    tCosh = 284,
    tTanh = 285,
    tFabs = 286,
    tAbs = 287,
    tFloor = 288,
    tCeil = 289,
    tRound = 290,
    tFmod = 291,
    tModulo = 292,
    tHypot = 293,
    tList = 294,
    tLinSpace = 295,
    tLogSpace = 296,
    tListFromFile = 297,
    tCatenary = 298,
    tPrintf = 299,
    tError = 300,
    tStr = 301,
    tSprintf = 302,
    tStrCat = 303,
    tStrPrefix = 304,
    tStrRelative = 305,
    tStrReplace = 306,
    tAbsolutePath = 307,
    tDirName = 308,
    tStrSub = 309,
    tStrLen = 310,
    tFind = 311,
    tStrFind = 312,
    tStrCmp = 313,
    tStrChoice = 314,
    tUpperCase = 315,
    tLowerCase = 316,
    tLowerCaseIn = 317,
    tTextAttributes = 318,
    tBoundingBox = 319,
    tDraw = 320,
    tSetChanged = 321,
    tToday = 322,
    tFixRelativePath = 323,
    tCurrentDirectory = 324,
    tSyncModel = 325,
    tNewModel = 326,
    tOnelabAction = 327,
    tOnelabRun = 328,
    tCodeName = 329,
    tCpu = 330,
    tMemory = 331,
    tTotalMemory = 332,
    tCreateTopology = 333,
    tCreateGeometry = 334,
    tRenumberMeshNodes = 335,
    tRenumberMeshElements = 336,
    tDistanceFunction = 337,
    tDefineConstant = 338,
    tUndefineConstant = 339,
    tDefineNumber = 340,
    tDefineStruct = 341,
    tNameStruct = 342,
    tDimNameSpace = 343,
    tAppend = 344,
    tDefineString = 345,
    tSetNumber = 346,
    tSetTag = 347,
    tSetString = 348,
    tPoint = 349,
    tCircle = 350,
    tEllipse = 351,
    tCurve = 352,
    tSphere = 353,
    tPolarSphere = 354,
    tSurface = 355,
    tSpline = 356,
    tVolume = 357,
    tBox = 358,
    tCylinder = 359,
    tCone = 360,
    tTorus = 361,
    tEllipsoid = 362,
    tQuadric = 363,
    tShapeFromFile = 364,
    tRectangle = 365,
    tDisk = 366,
    tWire = 367,
    tGeoEntity = 368,
    tCharacteristic = 369,
    tLength = 370,
    tParametric = 371,
    tElliptic = 372,
    tRefineMesh = 373,
    tAdaptMesh = 374,
    tRelocateMesh = 375,
    tReorientMesh = 376,
    tSetFactory = 377,
    tThruSections = 378,
    tWedge = 379,
    tFillet = 380,
    tChamfer = 381,
    tPlane = 382,
    tRuled = 383,
    tTransfinite = 384,
    tPhysical = 385,
    tCompound = 386,
    tPeriodic = 387,
    tParent = 388,
    tUsing = 389,
    tPlugin = 390,
    tDegenerated = 391,
    tRecursive = 392,
    tRotate = 393,
    tTranslate = 394,
    tSymmetry = 395,
    tDilate = 396,
    tExtrude = 397,
    tLevelset = 398,
    tAffine = 399,
    tBooleanUnion = 400,
    tBooleanIntersection = 401,
    tBooleanDifference = 402,
    tBooleanSection = 403,
    tBooleanFragments = 404,
    tThickSolid = 405,
    tRecombine = 406,
    tSmoother = 407,
    tSplit = 408,
    tDelete = 409,
    tCoherence = 410,
    tIntersect = 411,
    tMeshAlgorithm = 412,
    tReverseMesh = 413,
    tLayers = 414,
    tScaleLast = 415,
    tHole = 416,
    tAlias = 417,
    tAliasWithOptions = 418,
    tCopyOptions = 419,
    tQuadTriAddVerts = 420,
    tQuadTriNoNewVerts = 421,
    tRecombLaterals = 422,
    tTransfQuadTri = 423,
    tText2D = 424,
    tText3D = 425,
    tInterpolationScheme = 426,
    tTime = 427,
    tCombine = 428,
    tBSpline = 429,
    tBezier = 430,
    tNurbs = 431,
    tNurbsOrder = 432,
    tNurbsKnots = 433,
    tColor = 434,
    tColorTable = 435,
    tFor = 436,
    tIn = 437,
    tEndFor = 438,
    tIf = 439,
    tElseIf = 440,
    tElse = 441,
    tEndIf = 442,
    tExit = 443,
    tAbort = 444,
    tField = 445,
    tReturn = 446,
    tCall = 447,
    tSlide = 448,
    tMacro = 449,
    tShow = 450,
    tHide = 451,
    tGetValue = 452,
    tGetStringValue = 453,
    tGetEnv = 454,
    tGetString = 455,
    tGetNumber = 456,
    tUnique = 457,
    tHomology = 458,
    tCohomology = 459,
    tBetti = 460,
    tExists = 461,
    tFileExists = 462,
    tGetForced = 463,
    tGetForcedStr = 464,
    tGMSH_MAJOR_VERSION = 465,
    tGMSH_MINOR_VERSION = 466,
    tGMSH_PATCH_VERSION = 467,
    tGmshExecutableName = 468,
    tSetPartition = 469,
    tNameToString = 470,
    tStringToName = 471,
    tAFFECTPLUS = 472,
    tAFFECTMINUS = 473,
    tAFFECTTIMES = 474,
    tAFFECTDIVIDE = 475,
    tOR = 476,
    tAND = 477,
    tEQUAL = 478,
    tNOTEQUAL = 479,
    tLESSOREQUAL = 480,
    tGREATEROREQUAL = 481,
    tLESSLESS = 482,
    tGREATERGREATER = 483,
    tPLUSPLUS = 484,
    tMINUSMINUS = 485,
    UNARYPREC = 486
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 166 "Gmsh.y" /* yacc.c:352  */

  char *c;
  int i;
  unsigned int u;
  double d;
  double v[5];
  Shape s;
  List_T *l;
  struct TwoChar c2;

#line 529 "Gmsh.tab.cpp" /* yacc.c:352  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE gmsh_yylval;

int gmsh_yyparse (void);

#endif /* !YY_GMSH_YY_GMSH_TAB_HPP_INCLUDED  */



#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  5
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   16636

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  254
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  113
/* YYNRULES -- Number of rules.  */
#define YYNRULES  610
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  2177

#define YYUNDEFTOK  2
#define YYMAXUTOK   486

/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  ((unsigned) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   239,     2,   251,     2,   236,   238,     2,
     244,   245,   234,   232,   253,   233,   250,   235,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     226,     2,   228,   221,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   246,     2,   247,   243,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   248,   237,   249,   252,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   222,   223,   224,   225,
     227,   229,   230,   231,   240,   241,   242
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   268,   268,   269,   274,   276,   280,   281,   282,   283,
     302,   303,   304,   305,   306,   307,   308,   309,   310,   311,
     312,   313,   314,   315,   316,   317,   321,   325,   332,   337,
     342,   356,   369,   382,   410,   424,   437,   450,   469,   474,
     475,   476,   477,   478,   482,   484,   489,   491,   497,   601,
     496,   619,   626,   637,   636,   654,   661,   672,   671,   688,
     705,   728,   727,   741,   742,   743,   744,   745,   749,   750,
     756,   756,   757,   757,   763,   764,   765,   766,   771,   777,
     839,   854,   883,   893,   898,   906,   911,   919,   928,   933,
     945,   962,   968,   977,   995,  1013,  1022,  1034,  1039,  1047,
    1067,  1090,  1110,  1118,  1125,  1132,  1154,  1177,  1215,  1236,
    1248,  1262,  1262,  1264,  1266,  1275,  1285,  1284,  1305,  1304,
    1322,  1332,  1331,  1345,  1347,  1355,  1361,  1366,  1393,  1395,
    1398,  1400,  1404,  1405,  1409,  1421,  1434,  1449,  1458,  1471,
    1473,  1477,  1478,  1483,  1491,  1500,  1508,  1522,  1540,  1544,
    1551,  1560,  1563,  1569,  1573,  1585,  1588,  1595,  1618,  1634,
    1650,  1687,  1727,  1743,  1759,  1791,  1807,  1824,  1840,  1890,
    1908,  1914,  1920,  1927,  1958,  1973,  1995,  2018,  2041,  2064,
    2088,  2112,  2136,  2162,  2179,  2195,  2213,  2231,  2238,  2247,
    2246,  2276,  2278,  2280,  2282,  2284,  2292,  2294,  2296,  2298,
    2306,  2308,  2310,  2318,  2320,  2322,  2324,  2334,  2350,  2366,
    2382,  2398,  2414,  2431,  2468,  2490,  2514,  2515,  2520,  2523,
    2527,  2544,  2564,  2584,  2603,  2630,  2649,  2670,  2685,  2701,
    2719,  2770,  2791,  2813,  2836,  2941,  2957,  2992,  3014,  3036,
    3042,  3057,  3085,  3097,  3106,  3113,  3125,  3145,  3149,  3154,
    3158,  3163,  3170,  3177,  3184,  3196,  3269,  3287,  3312,  3327,
    3360,  3372,  3396,  3400,  3405,  3412,  3417,  3427,  3432,  3438,
    3446,  3450,  3454,  3458,  3462,  3471,  3535,  3551,  3568,  3585,
    3607,  3629,  3664,  3672,  3680,  3686,  3693,  3700,  3720,  3746,
    3758,  3769,  3787,  3805,  3824,  3823,  3848,  3847,  3874,  3873,
    3898,  3897,  3920,  3936,  3953,  3970,  3993,  4021,  4024,  4030,
    4042,  4062,  4066,  4070,  4074,  4078,  4082,  4086,  4090,  4099,
    4112,  4113,  4114,  4115,  4116,  4120,  4121,  4122,  4123,  4124,
    4127,  4151,  4170,  4193,  4196,  4212,  4215,  4232,  4235,  4241,
    4244,  4251,  4254,  4261,  4283,  4324,  4368,  4407,  4432,  4444,
    4456,  4468,  4480,  4489,  4519,  4545,  4571,  4603,  4630,  4656,
    4682,  4708,  4734,  4756,  4767,  4815,  4869,  4884,  4896,  4907,
    4914,  4929,  4943,  4944,  4945,  4949,  4955,  4967,  4985,  5013,
    5014,  5015,  5016,  5017,  5018,  5019,  5020,  5021,  5028,  5029,
    5030,  5031,  5032,  5033,  5034,  5035,  5036,  5037,  5038,  5039,
    5040,  5041,  5042,  5043,  5044,  5045,  5046,  5047,  5048,  5049,
    5050,  5051,  5052,  5053,  5054,  5055,  5056,  5057,  5058,  5059,
    5060,  5061,  5062,  5063,  5064,  5065,  5074,  5075,  5076,  5077,
    5078,  5079,  5080,  5081,  5082,  5083,  5084,  5089,  5088,  5096,
    5098,  5103,  5108,  5112,  5117,  5122,  5126,  5130,  5134,  5138,
    5142,  5146,  5152,  5168,  5173,  5179,  5185,  5204,  5225,  5258,
    5262,  5267,  5271,  5275,  5279,  5284,  5289,  5299,  5309,  5314,
    5325,  5334,  5339,  5344,  5372,  5373,  5379,  5380,  5386,  5385,
    5408,  5410,  5415,  5424,  5426,  5432,  5433,  5438,  5442,  5446,
    5450,  5454,  5461,  5465,  5469,  5473,  5480,  5485,  5492,  5497,
    5501,  5506,  5510,  5518,  5529,  5533,  5537,  5551,  5559,  5567,
    5574,  5584,  5607,  5612,  5618,  5623,  5629,  5640,  5646,  5652,
    5658,  5668,  5678,  5688,  5700,  5704,  5709,  5721,  5725,  5729,
    5733,  5751,  5759,  5767,  5796,  5806,  5822,  5833,  5838,  5842,
    5846,  5858,  5862,  5874,  5891,  5901,  5905,  5920,  5925,  5932,
    5936,  5941,  5955,  5971,  5975,  5979,  5983,  5987,  5995,  6001,
    6010,  6014,  6018,  6026,  6032,  6038,  6042,  6050,  6058,  6065,
    6074,  6078,  6082,  6097,  6111,  6125,  6137,  6153,  6162,  6171,
    6181,  6192,  6200,  6208,  6212,  6231,  6238,  6244,  6251,  6259,
    6258,  6268,  6292,  6294,  6300,  6305,  6307,  6312,  6317,  6322,
    6324,  6328,  6340,  6354,  6358,  6365,  6373,  6381,  6392,  6394,
    6397
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "tDOUBLE", "tSTRING", "tBIGSTR", "tEND",
  "tAFFECT", "tDOTS", "tSCOPE", "tPi", "tMPI_Rank", "tMPI_Size",
  "tEuclidian", "tCoordinates", "tTestLevel", "tExp", "tLog", "tLog10",
  "tSqrt", "tSin", "tAsin", "tCos", "tAcos", "tTan", "tRand", "tAtan",
  "tAtan2", "tSinh", "tCosh", "tTanh", "tFabs", "tAbs", "tFloor", "tCeil",
  "tRound", "tFmod", "tModulo", "tHypot", "tList", "tLinSpace",
  "tLogSpace", "tListFromFile", "tCatenary", "tPrintf", "tError", "tStr",
  "tSprintf", "tStrCat", "tStrPrefix", "tStrRelative", "tStrReplace",
  "tAbsolutePath", "tDirName", "tStrSub", "tStrLen", "tFind", "tStrFind",
  "tStrCmp", "tStrChoice", "tUpperCase", "tLowerCase", "tLowerCaseIn",
  "tTextAttributes", "tBoundingBox", "tDraw", "tSetChanged", "tToday",
  "tFixRelativePath", "tCurrentDirectory", "tSyncModel", "tNewModel",
  "tOnelabAction", "tOnelabRun", "tCodeName", "tCpu", "tMemory",
  "tTotalMemory", "tCreateTopology", "tCreateGeometry",
  "tRenumberMeshNodes", "tRenumberMeshElements", "tDistanceFunction",
  "tDefineConstant", "tUndefineConstant", "tDefineNumber", "tDefineStruct",
  "tNameStruct", "tDimNameSpace", "tAppend", "tDefineString", "tSetNumber",
  "tSetTag", "tSetString", "tPoint", "tCircle", "tEllipse", "tCurve",
  "tSphere", "tPolarSphere", "tSurface", "tSpline", "tVolume", "tBox",
  "tCylinder", "tCone", "tTorus", "tEllipsoid", "tQuadric",
  "tShapeFromFile", "tRectangle", "tDisk", "tWire", "tGeoEntity",
  "tCharacteristic", "tLength", "tParametric", "tElliptic", "tRefineMesh",
  "tAdaptMesh", "tRelocateMesh", "tReorientMesh", "tSetFactory",
  "tThruSections", "tWedge", "tFillet", "tChamfer", "tPlane", "tRuled",
  "tTransfinite", "tPhysical", "tCompound", "tPeriodic", "tParent",
  "tUsing", "tPlugin", "tDegenerated", "tRecursive", "tRotate",
  "tTranslate", "tSymmetry", "tDilate", "tExtrude", "tLevelset", "tAffine",
  "tBooleanUnion", "tBooleanIntersection", "tBooleanDifference",
  "tBooleanSection", "tBooleanFragments", "tThickSolid", "tRecombine",
  "tSmoother", "tSplit", "tDelete", "tCoherence", "tIntersect",
  "tMeshAlgorithm", "tReverseMesh", "tLayers", "tScaleLast", "tHole",
  "tAlias", "tAliasWithOptions", "tCopyOptions", "tQuadTriAddVerts",
  "tQuadTriNoNewVerts", "tRecombLaterals", "tTransfQuadTri", "tText2D",
  "tText3D", "tInterpolationScheme", "tTime", "tCombine", "tBSpline",
  "tBezier", "tNurbs", "tNurbsOrder", "tNurbsKnots", "tColor",
  "tColorTable", "tFor", "tIn", "tEndFor", "tIf", "tElseIf", "tElse",
  "tEndIf", "tExit", "tAbort", "tField", "tReturn", "tCall", "tSlide",
  "tMacro", "tShow", "tHide", "tGetValue", "tGetStringValue", "tGetEnv",
  "tGetString", "tGetNumber", "tUnique", "tHomology", "tCohomology",
  "tBetti", "tExists", "tFileExists", "tGetForced", "tGetForcedStr",
  "tGMSH_MAJOR_VERSION", "tGMSH_MINOR_VERSION", "tGMSH_PATCH_VERSION",
  "tGmshExecutableName", "tSetPartition", "tNameToString", "tStringToName",
  "tAFFECTPLUS", "tAFFECTMINUS", "tAFFECTTIMES", "tAFFECTDIVIDE", "'?'",
  "tOR", "tAND", "tEQUAL", "tNOTEQUAL", "'<'", "tLESSOREQUAL", "'>'",
  "tGREATEROREQUAL", "tLESSLESS", "tGREATERGREATER", "'+'", "'-'", "'*'",
  "'/'", "'%'", "'|'", "'&'", "'!'", "tPLUSPLUS", "tMINUSMINUS",
  "UNARYPREC", "'^'", "'('", "')'", "'['", "']'", "'{'", "'}'", "'.'",
  "'#'", "'~'", "','", "$accept", "All", "GeoFormatItems", "GeoFormatItem",
  "SendToFile", "Printf", "View", "Views", "ElementCoords",
  "ElementValues", "Element", "$@1", "$@2", "Text2DValues", "Text2D",
  "$@3", "Text3DValues", "Text3D", "$@4", "InterpolationMatrix", "Time",
  "$@5", "NumericAffectation", "NumericIncrement", "LP", "RP",
  "Affectation", "Comma", "DefineConstants", "$@6", "$@7", "$@8",
  "UndefineConstants", "Enumeration", "FloatParameterOptionsOrNone",
  "FloatParameterOptionsOrNone_NoComma", "FloatParameterOptions",
  "FloatParameterOption", "CharParameterOptionsOrNone",
  "CharParameterOptions", "CharParameterOption",
  "PhysicalId_per_dim_entity", "SurfaceConstraints", "CircleOptions",
  "Shape", "$@9", "GeoEntity", "GeoEntity123", "GeoEntity12",
  "GeoEntity02", "Transform", "MultipleShape", "ListOfShapes", "LevelSet",
  "Delete", "Colorify", "SetPartition", "Visibility", "Command", "Slide",
  "Loop", "Extrude", "$@10", "$@11", "$@12", "$@13", "ExtrudeParameters",
  "ExtrudeParameter", "BooleanOperator", "BooleanOption", "Boolean",
  "BooleanShape", "TransfiniteType", "TransfiniteArrangement",
  "TransfiniteCorners", "RecombineAngle", "PeriodicTransform",
  "Constraints", "Coherence", "HomologyCommand", "Homology", "FExpr",
  "FExpr_Single", "$@14", "GetForced_Default", "GetForcedStr_Default",
  "DefineStruct", "$@15", "Struct_FullName", "tSTRING_Member", "Append",
  "AppendOrNot", "VExpr", "VExpr_Single", "RecursiveListOfListOfDouble",
  "ListOfDouble", "ListOfDoubleOrAll", "FExpr_Multi",
  "RecursiveListOfDouble", "ColorExpr", "ListOfColor",
  "RecursiveListOfColor", "StringExprVar", "StringExpr", "$@16",
  "NameStruct_Arg", "Str_BracedRecursiveListOfStringExprVar",
  "BracedOrNotRecursiveListOfStringExprVar",
  "BracedRecursiveListOfStringExprVar", "RecursiveListOfStringExprVar",
  "MultiStringExprVar", "StringIndex", "String__Index", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,    63,   476,   477,   478,   479,    60,   480,    62,   481,
     482,   483,    43,    45,    42,    47,    37,   124,    38,    33,
     484,   485,   486,    94,    40,    41,    91,    93,   123,   125,
      46,    35,   126,    44
};
# endif

#define YYPACT_NINF -1836

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-1836)))

#define YYTABLE_NINF -558

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
   12583,    55,    44, 12739, -1836, -1836,  -117,    99,    63,  -122,
     -93,     4,   184,   209,   219,   223,     5,   255,   293,   311,
     321,   -45,   148,    23,  -112,   488,  -112,   136,   160,   215,
      51,   228,   235,    62,   270,   299,   310,   316,   340,   347,
     369,   375,   417,   433,   197,   379,   474,   456,   435,   476,
     590,   460,  6658,   479,   501,   552,   600,   -46,   384,   467,
      36,   471,   485,   721,   -69,   596,  -124,  -124,   635,   169,
     408,   665, -1836, -1836, -1836, -1836, -1836,   584,   250,   738,
     817,    19,    68,   835,   860,   377,   957,   970,   976,  5690,
     981,   745,   748,   749,    20,    67, -1836,   751,   752, -1836,
   -1836,   992,   993,   754, -1836,  2817,   758,  3249,    21,    25,
   -1836, -1836, -1836, 11480,   757, -1836, -1836, -1836, -1836, -1836,
     759, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836,
   -1836,   -44, -1836, -1836, -1836, -1836,    16, -1836,   999,   760,
    5446,   329,   761,  1004, 11480, 12912, 12912, -1836, 11480, -1836,
   -1836, -1836, -1836, 12912, -1836, -1836, -1836, -1836, -1836, -1836,
     762,   769,  1006, -1836, -1836,  3758,   772,   773,   774,   775,
      23, 11480, 11480, 11480,   776, 11480, 11480, 11480,   777, 11480,
   11480, 11480, 11480, 11480, 11480, 11480, 12912, 11480, 11480, 11480,
   11480,  5690,   779, -1836,  9310, -1836, -1836, -1836,   778,  5690,
    6900, 12912, -1836, -1836, -1836, -1836, -1836,  -112,  -112,  -112,
    -112,  -112,  -112,  -112,  -112,  -112,  -112,  -112,  -112,  -112,
    -112,  -112,  -112,  -112,  -112,  -112,  -112,  -112,  -112,  -112,
     192,  -112,  -112,  -112,  -112,  -112,   780,  -112,  -112,   781,
     467, -1836, -1836, -1836,  -112,  -112,    31,   846,   847,   848,
     783,  6900,   910,   467,   467,   793,  -112,  -112,   795,   797,
     799, -1836, -1836, -1836, 11480,  7142, 11480, 11480,  7384,    23,
     864,    35, -1836, -1836,   800, -1836,  4626, -1836, -1836, -1836,
   -1836, -1836,    91, 11480,  9310,  9310,   803,   805,  7626,  5690,
    5690,  5690, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836,
     802,  7868,   804,  4443,  1047,  6900,   806,    20,   807,   808,
    -124,  -124,  -124, 11480, 11480,   -62, -1836,     9,  -124, 10202,
     362,   -24,   815,   816,   818,   819,   820,   821,   822,  9310,
   11480,  5690,  5690,  5690,   823,     1,  1062,   824, -1836,  1063,
    1065, -1836,   825,   826,   829, -1836, -1836,   833,  5690,   868,
     869,   870, -1836, 11480,  5932, -1836,  1066,  1086, 11480, 11480,
   11480,   281, 11480,   845, -1836,   936, 11480, 11480, 11480, -1836,
   -1836, 11480, -1836,  -112,  -112,  -112,   875,   876,   877,  -112,
    -112,  -112,  -112,  -112,  -112,  -112, -1836,  -112, -1836, -1836,
   -1836,  -112,  -112,   878,   880,  -112,   881, -1836,   882,  1120,
    1121,   883, -1836, -1836,  1123,  1122,  1147,  1148,  -112, 11480,
   14252,    98, 12912,  9310, 11480, -1836, -1836,  6900,  6900, -1836,
     912,  3758,   467,  1151, -1836, -1836, -1836, -1836, -1836, -1836,
   11480, 11480,    39,  6900,  1156,   513,   915,   544,   916,  1159,
      30,   918, -1836,   919, 12933, 11480, -1836,   952,  -177, -1836,
      50,  -163,  6352, -1836,   -58, -1836,    90,  -189,   121,  1078,
   -1836,    23,   938, 11480, 11480, 11480, 11480,   939, 14768, 14793,
   14818, 11480, 14843, 14868, 14893, 11480, 14918, 14943, 14968, 14993,
   15018, 15043, 15068,   924, 15093, 15118, 15143,  4564,  1186, 11480,
    9310,  4664, -1836,   114, 11480,  1188,  1190,   953, 11480, 11480,
   11480, 11480, 11480, 11480, 11480, 11480, 11480, 11480, 11480, 11480,
   11480, 11480, 11480, 11480,  9310, 11480, 11480, 11480, 11480, 11480,
   11480,  9310,  9310,   949, 11480, 11480, 12912, 11480, 12912,  6900,
   12912, 12912, 12912,   951, 11480,    56, -1836, 10281, 11480,  6900,
    5690,  6900, 12912, 12912,  9310,    23,  3758,    23,   958,  9310,
     958, -1836,   958, 15168, -1836,   182,   955,   103,  1138, -1836,
    1195, 11480, 11480, 11480, 11480, 11480, 11480, 11480, 11480, 11480,
   11480, 11480, 11480, 11480, 11480,  8110, 11480, 11480, 11480, 11480,
   11480,    23, 11480, 11480,  1202, -1836,   739, 15193,   268,   296,
   11480, 11480, 11480, -1836,  1200,  1201,  1201,   966, 11480, 11480,
    1208,  9310,  9310, 14280,   971,  1209, -1836,   972, -1836, -1836,
    -114, -1836, -1836,  6594,  6836,  -124,  -124,   329,   329,   -64,
   10202, 10202, 11480,  1837,   -48, -1836, 11480, 11480, 11480, 11480,
   11480, 11480, 11480, 11480, 11480,   414, 15218,  1212,  1216,  1221,
   11480,  1217, 11480, -1836, 11480, 10643, -1836, -1836,  9310,  9310,
    9310, 11480,  1224, 11480, 11480, 11480, 15243,   985, -1836, -1836,
   15268, 15293, 15318,  1054,  7078, -1836,   987,  4963, 15343, 15368,
   14363, 12912, 12912, 12912, 12912, 12912, 12912, 12912, 12912, 12912,
   11480, 12912, 12912, 12912, 12912,    12,  3758, 12912, 12912, 12912,
      23,    23, -1836, -1836,  9310, -1836,   988, 10966, -1836,   996,
   11896, 11480,   958, 11480, -1836,    23, 11480, 11480,  1202,   989,
     452, 15393, 11289,   997,   481, 11480,  1235,   998,  6900, 15418,
   14390,   141,  1000,  1240,  1242, -1836, -1836, -1836,  9310,   166,
   11480, -1836, -1836, -1836,    23, 11480, 11480,  1202,  1005, -1836,
    1008,   -31,   467,    36,   467, -1836,  1009, 13469, -1836,   115,
    9310,    23, 11480, 11480,  1252,  1247,  9310, 11480,  1253, 12912,
      23, 10518,  1252,  1254, -1836,    23,  1256, 12912, 11480,  1018,
    1019, -1836, 11480,  7320,  7562,  7804,  8046,  3758,  1257,  1259,
    1261, 15443,  1263,  1264,  1266, 15468,  1267,  1268,  1269,  1270,
    1271,  1274,  1275, -1836,  1276,  1277,  1278, -1836, 11480, 15493,
    9310,  1032,  9310, 13498, -1836, -1836,  1280, 14336, 14336, 14336,
   14336, 14336, 14336, 14336, 14336, 14336, 14336, 14336,  8288, 14336,
   14336, 14336, 14336,   914,   472, 14336, 14336, 14336,  8525,  8762,
    9004,  4664,  1043,  1042,    75,  9310,  9246,  9578,   472,  9908,
     472,  1037,  1038,  1039,   -15,  9310, 16393, -1836,   472,  1044,
   13527, 13556, -1836, -1836,  1041,  -150,   472,  -149,  1050,   257,
     483,  1289, -1836,  1252,   472,  1051,  1049,  5140,  5371,   625,
    1084,  1461,  1461,   475,   475,   475,   475,   475,   475,   419,
     419,  9310,    13, -1836,    13,    13,   958,   958,   958,  1052,
   15518, 14417,   410,   420,  9310, -1836,  1293,  1056,  1058, 15543,
   15568, 15593, 11480,  6900,  1319,  1320,  9965, 13585, 15618, -1836,
     484,   490,  9310,  1075, -1836, 11953, -1836, 12020, 12077,  -124,
   11480, 11480, -1836, -1836,  1079,  1081, 10202,  4877,  1198,   295,
    -124, 12144, 15643, 13614, 15668, 15693, 15718, 15743, 15768, 15793,
   15818,  1085,  1327, 11480,  1329, -1836, 11480, 15843, -1836, 14444,
   14471, -1836,   491,   492,   493, 13643, -1836, 14498, 14525, 10149,
   -1836, -1836,  1330,  1331,  1332,  1090, 11480, 12201, 11480, 11480,
   -1836, -1836,    34,   -81,   277,   -81,  1091,  1096,  1092,   472,
     472,  1094, 10233,   472,   472,   472,   472, 11480,   472,  1334,
   -1836,  1095,  1099,   358,   -75,  1102,   502, -1836, -1836, -1836,
   -1836, 14336,    13, 12268,  1100,   506,  1101,  1169,  1345,  1203,
   10604,  1105,  1107,  1350,  6900, 13672, -1836, 11480,  1352,   158,
      73,  3758, 11480,  1354,  1357,    29, -1836,   510,  1316,  1318,
    6900, 13701,    10,  1115, 15868, 14552,   451, 11480, 11480,  1124,
    1118,  1125,  1119,  8352, -1836, -1836, -1836, -1836, 12912,   142,
    1131, 15893, 14579, -1836,  1126, -1836,   163, 10474, -1836, -1836,
   -1836,  1137, -1836,  1127, -1836,    74, -1836, -1836, 16393, -1836,
    1368, 14336, 11480, 11480, 11480, 11480,   472,  -124,  6900,  6900,
    1367,  6900,  6900,  6900,  1384,  6900,  6900,  6900,  6900,  6900,
    6900,  6900,  6900,  6900,  6900,  1249,  1385,  9310,  4664, -1836,
   -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836,
   -1836, -1836, -1836, 11480, -1836, -1836, -1836, -1836, -1836, -1836,
   -1836, -1836, -1836, 11480, 11480, 11480, -1836, -1836, -1836,   535,
   11480, 11480, -1836, 11480, -1836,  6900, 12912, 12912, -1836,   536,
    1140, -1836, -1836, -1836,  1213, 11480, 11480, -1836, -1836, -1836,
    1252, -1836,  1252, 11480, 11480,  1149, -1836,  6900,  -112, -1836,
   11480, -1836, 11480, 11480,   537,  1252,   410,    -7, -1836, 11480,
   11480,   472,   543,  6900,  9310,  9310,  1389,  1390,  1391,  3978,
   -1836, -1836,  1393, -1836,  1155, 16393,  1152, -1836,  1394,  1400,
    1401,   545,  1405, -1836, 12325, -1836, -1836,   -39, 10797, 11111,
   -1836, -1836, 13730,  -176,  1298,  1407, 10841,  1165,  1408,  1174,
      32,    46,   -72, -1836,   -30, -1836,   295,  1412,  1414,  1415,
    1417,  1418,  1419,  1420,  1421,  1422,   329,  6900, 16393, -1836,
    1273,  1173,  1426,  1427,  1428,  1335,  1429, -1836,  1431,  1433,
   11480,  6900,  6900,  6900,  1436, 11436, -1836,  5626,   719,    49,
    1437, -1836,  9310, -1836, -1836, -1836, -1836, 12912, -1836, -1836,
   11480, 12912, -1836, -1836, -1836, -1836, 16393, -1836,  1191,  1192,
   12912, -1836, 12912, -1836,  1252, 12912,  1199, -1836,  1193, -1836,
    1252, 11480, 11480,  1204,   467,  1205, 10927, -1836,  1380,  1206,
    6900, -1836,  1194, -1836, 13759, -1836, -1836, 11480,  1441,    42,
   11480,  1442,  1444,  1847, -1836,  1449,    20,  1445,  1218,   472,
    -112,  -112,  1452, -1836, -1836,  1215,  1219,  1214, -1836,  1459,
   -1836, -1836, -1836, -1836, -1836,  1252,   356, 11523, 11480, 14606,
   15918, 11480,  8589, 11480,  9310,  1220,   553,  1460,   134,  1252,
   -1836,  1239, 11480,  1484, 11480,  1252, 11164,  9547,   472,  4927,
    1244,  1260, -1836,  1487, 15943, 15968, 15993, 16018,  1512,    22,
    1340,  1340,  6900,  1513,  1514,  1517,  6900,   -99,  1518,  1519,
    1520,  1523,  1524,  1525,  1526,  1528,  1529, -1836,  1531,   554,
   14336, 14336, 14336, 14336,   472, 11694, 12379, 12863,  1292,   472,
     472, -1836,  1368,   472, 16043, 14336,  1294,   147, 16393, 14336,
   -1836,  1532,   472, 12908, 16393, 16393, -1836,   549, -1836,  1536,
   -1836, 16068, 14633, -1836,   472,  1535,   556,   561,  6900,  6900,
    6900,  1538,  1537, -1836,   205, 11480,  6900,  1296,  1299,  1539,
     633, -1836, 11480, 11480, 11480,  1300,  1301,  1302,  1306, -1836,
    2029,  6900, -1836, 11480, -1836,  1548, -1836,  1550, -1836, -1836,
   10202,    -5,  6174, -1836,  1309,  1311,  1312,  1313,  1314,  1315,
    8826,  1321,  1559, -1836,  9310, -1836, -1836, -1836,  1323, 11480,
   -1836, -1836, 14660,  1561,  1562,  1395, -1836, 11480, 11480, 11480,
   -1836,  1566,  1568,  1569,   889,   396,  1324,  5395,  1326, 11480,
      38,   472,  1338,   472,  1347, -1836, -1836,  3758,   599, 11480,
    1328, -1836, -1836,  2212, -1836, -1836,  1369,  1574, -1836,  2487,
   -1836,   201,  1337,  1613,  2875, -1836, -1836, -1836,    20, -1836,
     562, -1836, 11480,   205, 11698, 11770, -1836,  1381, 11480, 11480,
    6900,  1378, -1836,   478,  1623,  1622, 16093,  1625,  1289, 16118,
    1383,   563, 16143,   568,  1629,  1630, -1836, -1836, 12912,  1396,
    1633, 16168, -1836, 12941,  1398, -1836,  5204, 16393, -1836,  1632,
    -112,  7384, -1836, -1836, -1836, -1836,  1368, -1836,  1640,  1641,
    1644,  1645, -1836, -1836,  -124,  1647,  1648,  1649, -1836, -1836,
   -1836,  1651,   -43,  1560,  1653, -1836, -1836, -1836, -1836, -1836,
   -1836, -1836, -1836, -1836,  1624,  1413, -1836, -1836, -1836, -1836,
   -1836, 11480, 11480, 11480, -1836, -1836, -1836,  1260, -1836, -1836,
   -1836, -1836, 11480,  1423,  1410, -1836, -1836, 11480, 11480, 11480,
     472,   410, -1836, -1836, -1836, -1836,  1416,  1425,  1659,   -99,
    1660, 11480, -1836,  6900, 16393,   984,  9310,  9310, 11480, -1836,
    9965, 13788, 16193,  4362,   329,   329, 11480, 11480, -1836,   399,
    1424, 16218, -1836, -1836, 13817,   -26, -1836,  1663,  1664,  6900,
    -124,  -124,  -124,  -124,  -124,  6416,  1665, -1836, -1836,   570,
   11480,  3201,  1668, -1836, -1836,  6900,  5868,   874, 16243, -1836,
   -1836, -1836, -1836,  9635, -1836, 12912, 11480, -1836, 12912, 16393,
    9877,  3758,  1430, -1836, -1836, -1836, -1836,  1434,  1447, 11480,
   11480, 13846, 11480, 11289, -1836, 11289,  6900, -1836, -1836,  3758,
   11480,  1669,  1674,    29, -1836,  1675, -1836,    20, 14687,  6900,
   12912,  1678,   472, -1836,  1448,   472, 11480, 12974, 13007,   578,
   -1836, 11480, 11480,   480, -1836,  1432, -1836,  1391,  1680,  1696,
    1394,  1699, -1836, -1836,  1700, 11480, -1836, -1836, 11480, 11243,
   -1836, -1836,  1463, 11770,   580,  4515,  1701, -1836, -1836, -1836,
   -1836, -1836,   673, -1836, -1836, -1836, -1836,  1465,  1467,  1470,
   -1836,  1702,  6900, 14336, 14336, 13040, 14336, -1836,  1464, 13073,
   16268, 14714, -1836, -1836,  9310,  9310, -1836,  1713, -1836, 16393,
    1714,  1473, -1836,   581,   586, 14308,  3763,  1716,  1475, -1836,
   -1836, 11480,  1476,  1477, 13875, 14741,  1718,  6900,  1720,  1480,
   11480, -1836, -1836,   616,   -22,   -16,   137,   139,   173,  9068,
     195, -1836,  1723, 13904, -1836, -1836,  1554, -1836, 11480, 11480,
   -1836, -1836,  9310,  3797,  1726,  1488, 14336,   472, 16393, -1836,
   -1836, -1836, -1836,    38, -1836,  3758, -1836, 13933,  1485,  1486,
    1498,  1730,  4094, -1836,  1741,  1731, -1836, -1836,  1501,  1747,
     617, -1836,  1748,  1749,   204, 16393, 11480, 11480,  1508,  6900,
     618, 16393, 16293, -1836, -1836, -1836, -1836, 16318, 13106, -1836,
    1140,  1192,  6900,   472, -1836, 11480,  3758,    23,  9310,  9310,
   11480,  1751,   623, -1836, -1836, 11480,  1410, -1836, 11480, -1836,
   -1836,   624,   626, -1836, -1836,  6900,   470,   521,  9310, -1836,
   -1836,   329,  6110, -1836, -1836, -1836,  1752, -1836,  1510,  6900,
   -1836, 13962,  1754,  9310,  -124,  -124,  -124,  -124,  -124, -1836,
   -1836, 11480, 13991, 14020,   631, -1836, -1836, -1836, -1836, -1836,
   -1836,  1516,  1756,  1515, -1836,  1758, -1836, -1836,    20, -1836,
    1585, -1836, -1836, -1836, -1836, -1836, 11480, 13139, 13172,  6900,
   -1836,  1760, 11480,  1522, -1836, 11480,  1527,  1530, -1836, -1836,
    4594, -1836,  1533,   632,   637, 14049, -1836,  1540, 13205,  1541,
   13238, -1836,  1542,   639,  1543,  -124,  6900,  1761,  1544,  -124,
    1762,   640,  1534, -1836, 11480, -1836,  1765,  1638, 12392,  1546,
   -1836,   645,   200,   211,   218,   243,   256,  4144, -1836, -1836,
    1767,  1768, -1836, -1836, -1836,  1769, -1836,  1547, 16393, 11480,
   11480,   648, -1836, 16393, 13271, -1836, -1836,  1140,  3758,  1551,
   -1836, -1836, -1836, 11480, 11480, -1836, 11480,  9310,  1772,  -124,
     122, -1836, -1836,  -124,   124, -1836,  1774, -1836, 14078, -1836,
   11480, -1836,   295, -1836,  1778,  9310,  9310,  9310,  9310,  9068,
   -1836, -1836, -1836, 11289, -1836, 11480, 16343, 13304,    41, 11480,
    1549, -1836, -1836, 13337, 13370, 13403,   654, -1836,   267, -1836,
     324, -1836, -1836, -1836,  4337,   332, 12449, -1836,   655,   663,
     666,   674,   326,   675,  1552,   676, -1836, 11480, -1836,  6900,
   14107, -1836, 11480, 11480, 11480, -1836,  -124,  -124, -1836, -1836,
   -1836,   295,  1779,  1781,  1783,  1791,  9310,  1793,  1794,  1796,
    1557, 16368,   677,  1800, 14136, 14336, 13436,   330,   333,   441,
   -1836, -1836, -1836, -1836,   682, -1836, -1836, -1836, 12912, -1836,
    1563, -1836,  1801, -1836, 11480, 11480, 11480, -1836,  1802,   684,
   -1836,  1564,  6900, -1836, 14165, 14194, 14223, -1836,  1803, 12912,
   12912,   685, -1836,  1804,  1807, -1836, -1836,   710, -1836,  1808,
   -1836, -1836,  1809, 12912, -1836, -1836, -1836
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     2,     3,     1,   608,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   191,     0,     0,
     192,     0,     0,   193,     0,   194,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   320,   321,   322,   323,   324,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   281,     0,     0,   289,
     290,     0,     0,     0,   284,     0,     0,     0,     0,     0,
     372,   373,   374,     0,     0,     5,     6,     7,     8,    10,
       0,    11,    24,    12,    13,    14,    15,    23,    22,    21,
      16,     0,    17,    18,    19,    20,     0,    25,     0,   609,
       0,   218,     0,     0,     0,     0,     0,   266,     0,   268,
     269,   264,   265,     0,   270,   271,   272,   273,   113,   123,
     608,   485,   480,    70,    71,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   274,     0,   203,   204,   205,     0,     0,
       0,     0,   426,   427,   429,   430,   428,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   434,   435,   436,     0,     0,   191,   196,   197,   198,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   431,   432,   433,     0,     0,     0,     0,     0,     0,
       0,     0,   520,   521,     0,   522,   498,   379,   439,   442,
     303,   499,   480,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   191,   192,   193,   194,   189,   196,   197,   198,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   487,     0,     0,   218,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   608,     0,     0,   218,     0,
       0,   369,     0,     0,     0,   200,   201,     0,     0,     0,
       0,     0,   506,     0,     0,   504,     0,     0,     0,     0,
       0,   608,     0,     0,   543,     0,     0,     0,     0,   262,
     263,     0,   560,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   562,     0,   586,   564,
     565,     0,     0,     0,     0,     0,     0,   563,     0,     0,
       0,     0,   282,   283,     0,   218,     0,   218,     0,     0,
       0,   480,     0,     0,     0,   218,   375,     0,     0,    76,
       0,    63,     0,     0,    64,    65,    66,    67,    68,    69,
      70,    71,     0,     0,     0,     0,     0,     0,     0,   549,
     480,     0,   217,     0,   216,     0,   170,     0,     0,   549,
     550,     0,     0,   598,     0,   599,   550,   111,   111,     0,
     478,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   537,   538,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    70,    71,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   513,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   382,     0,
     381,   507,   383,     0,   500,     0,     0,   480,     0,   515,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    70,    71,     0,   456,     0,     0,     0,     0,
       0,     0,     0,   304,     0,   337,   337,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   218,     0,   218,   218,
       0,   489,   488,     0,     0,     0,     0,   218,   218,     0,
       0,     0,     0,   300,     0,   218,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   339,     0,     0,
       0,     0,     0,   243,     0,     0,   241,   370,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   347,   261,
       0,     0,     0,     0,     0,   218,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   286,   285,     0,   248,     0,     0,   250,     0,
       0,     0,   381,     0,   218,     0,     0,     0,     0,     0,
       0,     0,   325,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    80,    72,    73,     0,     0,
       0,   259,    38,   255,     0,     0,     0,     0,     0,   213,
       0,     0,     0,     0,     0,   219,     0,     0,   171,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   112,     0,     0,     0,   483,     0,
       0,   481,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   331,     0,     0,     0,   195,     0,     0,
       0,     0,     0,     0,   365,   366,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   480,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   437,   455,     0,     0,
       0,     0,   516,   517,     0,     0,     0,     0,     0,   474,
       0,   380,   501,     0,     0,     0,     0,   509,     0,   399,
     398,   396,   397,   392,   394,   393,   395,   401,   400,   385,
     384,     0,   386,   508,   387,   390,   388,   389,   391,   481,
       0,     0,   482,   459,     0,   523,     0,     0,     0,     0,
       0,     0,     0,     0,   335,     0,     0,     0,     0,   368,
       0,     0,     0,     0,   367,     0,   218,     0,     0,     0,
       0,     0,   491,   490,     0,     0,     0,     0,     0,     0,
       0,   294,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   354,     0,     0,   242,     0,
       0,   237,     0,     0,     0,     0,   364,     0,     0,     0,
     380,   505,     0,     0,     0,     0,     0,     0,     0,     0,
     287,   288,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     589,     0,     0,     0,   476,     0,     0,   247,   251,   249,
     253,     0,   386,     0,   481,   459,   610,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    87,     0,     0,   380,
       0,    63,     0,     0,     0,     0,    79,     0,    63,    64,
       0,     0,     0,   481,     0,     0,   459,     0,     0,     0,
     189,     0,     0,     0,   605,    28,    26,    27,     0,     0,
       0,     0,     0,   482,   553,    29,     0,     0,   256,   600,
     601,     0,   602,   553,    74,   114,    75,   124,   484,   486,
     130,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   539,   540,
     206,     9,   403,   404,   405,   406,   407,   408,   409,   410,
     411,   425,   412,     0,   414,   415,   416,   417,   418,   536,
     419,   420,   421,     0,     0,     0,   528,   527,   526,     0,
       0,     0,   533,     0,   471,     0,     0,     0,   473,     0,
     128,   454,   512,   511,   199,     0,     0,   440,   535,   445,
       0,   451,     0,     0,     0,     0,   502,     0,     0,   452,
       0,   514,     0,     0,     0,     0,   444,   443,   466,    70,
      71,     0,     0,     0,     0,     0,     0,     0,   380,   333,
     338,   336,     0,   346,     0,   148,   149,   199,   380,     0,
       0,     0,     0,   238,     0,   252,   254,     0,     0,     0,
     207,   209,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   307,     0,   291,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   218,     0,   340,   353,
       0,     0,     0,     0,     0,     0,     0,   202,     0,     0,
       0,     0,     0,     0,     0,     0,   244,     0,     0,     0,
       0,   576,     0,   583,   572,   573,   574,     0,   588,   587,
       0,     0,   577,   578,   579,   585,   592,   591,     0,   139,
       0,   566,     0,   568,     0,     0,     0,   561,     0,   246,
       0,     0,     0,     0,     0,     0,     0,   326,     0,     0,
       0,   376,     0,   606,     0,   101,    63,     0,     0,     0,
       0,     0,     0,     0,    95,     0,     0,     0,     0,     0,
       0,     0,     0,   558,    48,     0,     0,     0,    61,     0,
      39,    40,    41,    42,    43,     0,   444,   443,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     552,   551,     0,     0,     0,     0,     0,     0,     0,   135,
       0,   131,   132,     0,     0,     0,     0,     0,     0,     0,
     155,   155,     0,     0,     0,     0,     0,   151,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   343,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   519,     0,     0,     0,     0,     0,   474,   475,     0,
     447,     0,     0,     0,   510,   402,   503,   460,   458,     0,
     457,     0,     0,   524,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   345,     0,     0,     0,     0,     0,     0,
       0,   245,     0,     0,     0,     0,     0,     0,     0,   312,
       0,     0,   311,     0,   314,     0,   316,     0,   301,   308,
       0,     0,     0,   236,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   355,     0,   240,   239,   371,     0,     0,
      35,    36,     0,     0,     0,     0,   544,     0,     0,     0,
     277,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   476,   477,   570,     0,   460,     0,
       0,   218,   327,     0,   328,   218,     0,     0,   559,     0,
      86,     0,     0,     0,     0,    84,    91,    93,     0,   547,
       0,    99,     0,     0,     0,     0,    81,     0,     0,     0,
       0,     0,    34,   460,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    30,    31,     0,   554,
       0,     0,    32,     0,   554,   603,     0,   115,   120,     0,
       0,     0,   134,   137,   138,   479,     0,    77,     0,     0,
       0,     0,    78,   157,     0,     0,     0,     0,   158,   173,
     174,     0,     0,     0,     0,   159,   184,   175,   179,   180,
     176,   177,   178,   165,     0,     0,   413,   422,   423,   424,
     529,     0,     0,     0,   469,   470,   472,   129,   438,   468,
     441,   446,     0,     0,   474,   185,   453,     0,    70,    71,
       0,   465,   461,   463,   530,   181,     0,     0,     0,   151,
       0,     0,   344,     0,   150,     0,     0,     0,     0,   260,
       0,     0,     0,     0,   218,   218,     0,     0,   313,   498,
       0,     0,   315,   317,     0,     0,   295,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   212,   182,     0,
       0,     0,     0,   162,   163,     0,     0,     0,     0,   102,
     103,   104,   108,     0,   584,     0,     0,   582,     0,   593,
       0,     0,   140,   141,   590,   567,   569,     0,     0,     0,
       0,     0,     0,   325,   329,   325,     0,   377,    85,    63,
       0,     0,     0,     0,    83,     0,   545,     0,     0,     0,
       0,     0,     0,   596,   595,     0,     0,     0,     0,     0,
     496,     0,     0,   465,   257,   461,   258,     0,     0,     0,
       0,     0,   223,   220,     0,     0,   557,   555,     0,     0,
     116,   121,     0,     0,     0,   537,   538,   133,   348,   349,
     350,   351,   156,   160,   161,   166,   183,     0,     0,     0,
     168,     0,     0,     0,     0,     0,     0,   448,     0,     0,
       0,     0,   525,   467,     0,     0,   167,     0,   186,   334,
       0,     0,   187,     0,     0,     0,     0,     0,     0,   495,
     494,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   235,   234,     0,     0,     0,     0,     0,     0,     0,
       0,   225,     0,     0,   352,    37,     0,   542,     0,     0,
     279,   278,     0,     0,     0,     0,     0,     0,   143,   144,
     147,   146,   145,     0,   571,     0,   607,     0,     0,     0,
       0,     0,     0,    96,     0,     0,    97,   548,     0,     0,
       0,    88,     0,     0,     0,    44,     0,     0,     0,     0,
       0,    46,     0,   224,   221,   222,    33,     0,     0,   604,
     128,   139,     0,     0,   136,     0,     0,     0,     0,     0,
       0,     0,     0,   531,   532,     0,   474,   449,     0,   462,
     464,     0,     0,   169,   190,     0,   341,   341,     0,   109,
     110,   218,     0,   210,   211,   302,     0,   309,     0,     0,
     218,     0,     0,     0,     0,     0,     0,     0,     0,   215,
     214,     0,     0,     0,     0,   105,   106,   575,   581,   580,
     142,     0,     0,     0,   330,     0,    92,    94,     0,   100,
       0,    82,   597,    89,    90,    49,     0,     0,     0,     0,
     497,     0,     0,   462,   556,     0,     0,     0,   118,   594,
       0,   125,     0,     0,     0,     0,   172,     0,     0,     0,
       0,   305,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   493,     0,   319,     0,     0,   296,     0,
     226,     0,     0,     0,     0,     0,     0,     0,   541,   280,
       0,     0,   363,   218,   378,     0,   546,     0,    45,     0,
       0,     0,    62,    47,     0,   117,   122,   128,     0,     0,
     153,   154,   152,     0,     0,   450,     0,     0,     0,     0,
       0,   342,   356,     0,     0,   357,     0,   208,     0,   310,
       0,   292,     0,   218,     0,     0,     0,     0,     0,     0,
     164,   107,   276,   325,    98,     0,     0,     0,     0,     0,
       0,   126,   127,     0,     0,     0,     0,   188,     0,   360,
       0,   361,   362,   492,     0,     0,   298,   229,     0,     0,
       0,     0,     0,     0,     0,     0,    53,     0,    59,     0,
       0,   119,     0,     0,     0,   306,     0,     0,   318,   297,
     293,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     230,   231,   232,   233,     0,   227,   332,    50,     0,    57,
       0,   267,     0,   534,     0,     0,     0,   299,     0,     0,
      51,     0,     0,   275,     0,     0,     0,   228,     0,     0,
       0,     0,   518,     0,     0,    54,    52,     0,    55,     0,
     358,   359,     0,     0,    60,    58,    56
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1836, -1836, -1836, -1836,   482, -1836, -1836, -1836, -1836,  -248,
   -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836, -1836,
   -1836, -1836,  -710,  -109,  4025,  3321, -1836,  1363, -1836, -1836,
   -1836, -1836, -1836, -1836, -1835, -1836,   430,   258,   -56, -1836,
     -20, -1836,   198,   465,  1825, -1836,   294,   -54, -1836, -1836,
       0,  -610,  -301, -1836, -1836, -1836, -1836, -1836, -1836, -1836,
   -1836,  1827, -1836, -1836, -1836, -1836, -1212, -1211,  1828, -1690,
    1829, -1836, -1836, -1836,  1237, -1836,   -73, -1836, -1836, -1836,
   -1836,  2052, -1836, -1836, -1388,   341,  1833, -1836,     6,  -687,
   -1836, -1836,    78, -1836, -1647,   503,  -174,  2718,  2725,  -305,
     126, -1836,   701,   -42, -1836, -1836,   150,   318, -1642,  -121,
    1082, -1836,    -3
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     3,   115,  1048,   116,   117,  1032,  1864,  1870,
    1320,  1527,  2017,  2149,  1321,  2120,  2167,  1322,  2151,  1323,
    1324,  1531,   433,   585,   586,  1118,   118,   765,   457,  1880,
    2027,  1881,   458,  1754,  1393,  1350,  1351,  1352,  1491,  1692,
    1693,  1184,  1584,  1575,   745,   597,   270,   271,   348,   199,
     272,   443,   444,   122,   123,   124,   125,   126,   127,   128,
     129,   273,  1216,  2052,  2111,   929,  1212,  1213,   274,  1011,
     275,   133,  1422,  1182,   904,   944,  1987,   134,   135,   136,
     137,   276,   277,  1140,  1155,  1276,   278,   770,   279,   893,
     769,   460,   612,   316,  1729,   355,   356,   281,   555,   363,
    1307,  1520,   453,   449,  1269,   988,  1564,  1722,  1723,   973,
     455,   139,   411
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     140,  1449,   607,   121,  1451,   641,   301,   924,   925,  1613,
     147,  1022,  1813,  1848,  1314,  1849,   160,   488,   623,  1030,
     162,  1005,   416,   335,   361,   495,   404,   160,  1573,   161,
     406,   434,   454,   160,   160,  1582,   536,   645,  1444,   734,
     559,  1249,  1690,   721,     5,  1966,  1511,  2098,  1841,   282,
    1036,  1767,  1446,  1481,   287,   174,   615,   616,   763,   751,
     160,     4,  1205,   399,   764,   403,   178,  1054,   749,   287,
    1768,   365,   340,  1436,   341,  1063,   750,   288,   339,  1206,
    1296,  1347,   755,  1583,   581,   306,   282,  1207,  1208,  1209,
     756,   364,  1039,  1210,  1211,   726,  1149,   727,   439,   760,
     581,  1150,   400,  1146,   697,   336,   700,   705,   311,   312,
     307,  1883,   461,   143,   712,   594,   595,   596,   615,   616,
     313,  1045,   145,   462,   314,  1892,   308,   309,  2079,  1205,
    2081,   141,   163,   297,   164,   142,   298,   440,   299,   919,
    1547,   442,   450,   450,   315,   317,  1206,   320,  1021,   300,
     456,   146,  1482,  1483,  1207,  1208,  1209,   637,   638,   639,
    1210,  1211,   342,   144,   726,  1296,   727,   467,   615,   616,
     615,   616,   759,  1028,   652,  1274,  1158,  1448,  1275,  1315,
    1316,  1317,  1318,   450,   615,   616,   617,   758,   282,   926,
     149,   282,  2070,   615,   616,   759,   282,   282,   450,  1445,
     414,   158,   615,   616,   415,   930,   615,   616,  1709,   337,
     615,   616,  1296,  1447,  1432,   150,   615,   616,   722,   723,
     626,  1512,  1513,  1450,   627,   151,  1778,  1810,   114,   152,
    1138,  1923,  1691,   428,   429,   114,   114,  1924,   759,   114,
    1449,   615,   616,  1409,  1656,   114,   114,   642,   282,   153,
     578,   579,   148,   142,   615,   616,   580,   618,  1983,  1319,
     417,   154,   282,   987,   418,   282,   557,   338,   362,   405,
     428,   429,   114,   407,   735,   556,   736,  1306,  -550,   537,
     737,   282,   282,   560,  1250,   282,   282,   282,   282,  2099,
     424,   425,   426,   427,   752,   175,   753,   120,   282,   155,
     754,   726,   282,   727,   364,   915,   179,   917,   918,  1297,
    1300,   366,  2021,   428,   429,   428,   429,   156,   163,   582,
     164,   583,  1128,  1299,   931,   584,   282,   157,   282,   282,
     282,   428,   429,   441,   761,   582,   753,   583,   428,   429,
     762,   584,   706,  1046,   707,   282,  1047,   163,   708,   164,
     331,   282,   332,   296,   615,   616,   615,   616,   424,   425,
     426,   427,  1046,   801,   967,  1047,   852,   802,   766,   615,
     616,   615,   616,  2094,   764,   424,   425,   426,   427,   716,
     171,   428,   429,  1029,   425,   426,   427,  1338,   610,   611,
    1925,   163,  1926,   164,   159,   802,   619,   624,   428,   429,
    1153,   311,   312,  1003,   172,   615,   616,   561,  1343,   450,
     282,   844,   321,   313,   282,   282,   802,   319,   424,   425,
     426,   427,   424,   425,   426,   427,  1927,   615,   616,  1205,
     282,   862,   615,   616,   436,   802,   521,   738,   522,   436,
     436,   428,   429,   615,   616,   190,  1206,   436,  1928,  1955,
     615,   616,  2122,  2055,  1207,  1208,  1209,  1956,   771,   173,
    1210,  1211,   193,  1396,  2056,  1397,  1205,    65,    66,    67,
      68,  2057,   176,    71,   345,   615,   616,   346,  1407,   177,
     436,   289,    80,  1206,   290,    83,   291,   282,   615,   616,
     347,  1207,  1208,  1209,   191,   436,  2058,  1210,  1211,   615,
     616,   163,   322,   164,   858,  2161,   323,  1152,  1979,  2059,
    1153,   282,   324,   325,   180,   326,   327,   897,   282,   834,
    2106,   802,   726,   450,   727,   450,   282,   450,   450,   456,
    1252,   663,   848,   142,   533,   328,   282,   282,   282,   450,
     450,   282,   162,   181,   162,   898,   282,   540,   541,   802,
     731,   857,   975,   859,   182,   280,   615,   616,   615,   616,
     183,   292,   615,   616,   293,   615,   616,   294,   302,   295,
     195,   303,   282,   196,   192,  1205,   197,  2107,   889,  2116,
      44,  2109,   166,  2145,   184,   167,  2146,  1494,   168,   198,
     169,   185,  1206,  1498,   615,   616,   428,   429,   282,   282,
    1207,  1208,  1209,   726,  -552,   727,  1210,  1211,  1984,  1985,
     625,  1272,  1023,   186,  1986,  1194,  1461,   442,   442,   187,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,  1533,  1639,
    1640,   726,   580,   727,   990,   282,   282,   282,  1806,   802,
     428,   429,  1549,   703,   576,   577,   578,   579,  1554,  1988,
    1989,   188,   580,   941,  1169,  1986,  1170,   802,   456,   450,
     456,   450,   450,   450,   450,   450,   450,   189,   450,   450,
     450,   450,   989,   194,   450,   450,   450,   162,   995,  1041,
    2147,   282,   200,   922,   923,  1328,   994,  1282,   611,  -553,
     286,  1007,  1004,   496,   201,   802,   436,   573,   574,   703,
     576,   577,   578,   579,  1633,   282,   717,   726,   580,   727,
     428,   429,  1732,   283,  1700,   282,  -554,  1479,  -557,   304,
    1013,  1033,  1156,  1189,  1014,  1076,   802,   802,   746,  1190,
    1234,  1235,  1236,   802,   802,   802,   802,   282,  1050,   284,
    1281,  1278,  1282,   282,   280,   802,   456,  1061,   726,  1309,
     727,   728,  1065,   802,   450,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,  1168,  1384,  1391,  1406,   580,   802,   802,
     802,   593,  1414,  1618,  1429,  1619,   802,   282,   802,   282,
     285,  1710,  1545,  1595,   600,  1626,   802,   802,   605,   802,
    1627,  1716,  1739,  1719,   802,  1717,   802,  1741,   305,  1822,
     436,   802,   436,   802,   436,   436,   436,  1868,   330,  1884,
    1906,  1869,   282,  1885,   802,  1907,   436,   436,   333,   802,
    2085,   438,   282,  1699,   310,  1700,   448,   451,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,  1186,  1922,  1952,  1961,   580,  1869,
     759,  1962,  1977,  1981,  1449,  1982,  1869,   802,   282,   802,
    2010,  2030,  1829,   318,   802,   802,  2031,   483,  2038,  2046,
     802,   282,  1869,   802,  2054,  1682,  1683,  2068,   802,  2129,
     282,  1869,   497,  2105,  2112,   615,   616,   802,   802,   282,
    1298,  1301,  2113,   329,   334,  2114,   802,   746,  1449,   802,
     713,   714,   561,  2115,  2117,  2119,  2140,   802,   802,  1962,
    1869,  2148,   343,  2158,  2169,   802,   724,  2159,  1869,   746,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,   748,  2172,
     344,   349,   580,  2173,  1480,   436,   436,   436,   436,   436,
     436,   436,   436,   436,   350,   436,   436,   436,   436,  1302,
     351,   436,   436,   436,   726,   357,   727,   894,  1791,   358,
    1792,   746,   359,   360,   746,   367,   368,  1197,   369,   370,
     371,  1519,   401,   412,  1203,   419,   746,   413,  1214,   445,
     446,   282,   420,   459,   142,   461,   463,   464,   465,   466,
     471,   475,  1308,   489,   529,   532,   494,   282,  -192,  -193,
    -194,   538,   841,   539,  1802,  1803,  1040,   542,  1042,   545,
     282,   546,   593,   547,   853,   450,   558,   590,   415,   591,
     598,   604,   601,   436,   606,   608,   609,  1408,  1410,   628,
     629,   436,   630,   631,   632,   633,   634,   640,   643,   646,
     644,   647,   658,   648,   649,   282,   282,   650,   282,   282,
     282,   651,   282,   282,   282,   282,   282,   282,   282,   282,
     282,   282,   659,   665,   282,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,   709,   653,   654,   655,   580,   666,   674,
     675,   676,   687,  1830,   688,   690,   692,   693,   691,   695,
     696,   694,   282,   450,   450,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   800,   576,
     577,   578,   579,   698,   282,  1359,   699,   580,   718,   726,
     715,   727,   725,   730,   732,   733,   141,   768,   739,   793,
     282,   282,   282,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,   772,   777,   798,   804,   580,   805,   835,   806,   845,
    1703,   580,   865,   866,  1705,   863,   892,   902,   903,   746,
     906,   746,   746,  1715,   909,   914,   913,  1408,  1410,   943,
     916,  1018,   945,   948,   282,   746,   442,   838,   946,   840,
     956,   842,   843,   965,   961,   968,  1006,   997,   282,   282,
     282,  1016,  1012,   854,   855,   999,  1017,  1025,  1026,   282,
    1024,  1037,  1038,  1055,   450,  1377,  1053,  1043,   450,  1058,
    1064,   746,  1066,  1069,  1077,  1070,  1078,   450,  1079,   450,
    1081,  1082,   450,  1083,  1085,  1086,  1087,  1088,  1089,  1463,
    1097,  1090,  1091,  1092,  1093,  1094,  1101,   282,  1126,  1127,
    1135,  1136,  1137,  1142,  1145,  1151,  1157,   746,  1161,  1160,
    1173,  1992,  1165,   364,  1174,  1558,  1175,  1563,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,  1181,   450,  1192,  1183,   580,  1200,   282,
    1201,   282,  1204,  1226,  1227,  1229,  1255,  1241,  1242,  1243,
    1244,  1256,   436,  1268,  1271,  1257,   282,  1260,  1270,  1277,
    1280,  1284,  1285,  1283,  1289,  1290,  1291,  1286,  1295,   282,
    1304,  1305,  1310,   282,  1311,  1325,  1332,  1334,  1331,  1333,
    1342,  1346,  1349,   974,  1362,   976,   977,   978,   979,   980,
     981,  1339,   983,   984,   985,   986,  1504,  1345,   991,   992,
     993,  1366,  1378,  1392,  1400,  -195,  1418,  1419,  1420,  1423,
    1424,  1426,  1711,  1724,  1724,  1425,  1180,  1427,  1428,  1430,
    1437,  1438,  1857,  1441,  1442,   282,   282,   282,  1443,  1452,
    1453,  1464,  1454,   282,  1455,  1456,  1457,  1458,  1459,  1460,
     436,   436,  1465,  1466,  1467,  1468,  1469,  1470,   282,  1471,
    1476,  1484,  1489,  1507,  1496,  1490,  1497,  1510,  1515,   282,
    1516,  1521,  1499,  1501,  1505,  1698,  1518,   282,  1526,  1528,
    1059,   282,  1530,  1529,  1522,  1532,  1546,  1574,  1067,  1544,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,   746,  1550,
    1552,  1565,   580,  1567,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,  1783,  1566,  1751,   364,   580,  1292,  1572,  1578,
    1579,   456,   456,  1580,  1585,  1586,  1587,   282,  1655,  1588,
    1589,  1590,  1591,  1312,  1592,  1593,  1594,  1604,  1615,  1611,
    1621,  1625,  1631,  1632,  1636,   450,  1638,  1637,  1644,  1645,
    1646,   436,  1647,   282,  1652,   436,  1653,  1659,   282,  1660,
    1661,  1662,  1663,  1664,   436,  1668,   436,  1673,  1674,   436,
    1667,  1670,  1679,  1675,  1680,  1681,  1702,  1685,  1500,  1688,
    1707,  1360,  1361,  1695,  1363,  1364,  1365,  1712,  1367,  1368,
    1369,  1370,  1371,  1372,  1373,  1374,  1375,  1376,  1797,  1860,
    1275,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,  1998,
    1713,   436,  1706,   580,  1783,  1726,  1731,  1733,  1734,  1771,
     282,  1736,  1738,   282,   282,  1742,  1743,  1746,  1388,  1752,
    1745,  1834,  1749,  2015,   442,   442,  1758,  1759,  1839,  1842,
    1760,  1761,  1762,  1763,  1764,  1765,   282,  1766,  1769,  1770,
    1401,  1772,   282,  1153,  1784,  1786,  1788,  1851,  1777,  1811,
    1812,  1821,   282,  1785,  1825,  1853,  1415,  1807,  1854,  1844,
    -555,  1856,   450,  1843,  1861,   450,  1873,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
    1845,   759,  1874,   282,   580,  1875,  1876,  1891,  1887,  1897,
    1308,  1882,  2063,  1888,   364,  1889,   282,   456,  1890,  1903,
    1904,  1905,  1910,  1911,  1917,  1913,  1914,  1919,  1920,  1929,
    1462,  1931,  1936,  1937,  1943,  1944,  1946,  1949,  1814,  1815,
    1816,  1817,  1818,  1820,  1473,  1474,  1475,  1945,  1948,  1337,
     456,  1950,  2086,  1951,  1953,  1954,  1959,  1976,  1995,  1996,
    2000,  2011,  2012,  2013,  2014,  2016,  2022,  2042,  2045,   282,
    -556,  2049,  2050,  2061,  2062,  2064,  2025,  2029,  2077,  2026,
    2082,   282,   282,  2047,  2087,  2130,  2035,  2131,  2033,  2132,
    2037,  2039,  2043,  1506,  2053,  2065,  2072,  2133,  2101,  2135,
    2136,  2118,  2137,  1941,   282,  2138,  2141,  2153,  2157,  2165,
    2170,  2152,  2160,  2171,  2174,  2175,   282,  2095,   436,   436,
    1548,   767,  1607,  1940,  1757,  1967,  1576,  1787,   119,   282,
     130,   131,   132,   905,  1990,  1697,   138,  1389,  1390,  1855,
    1840,  1060,   436,  1725,  1971,     0,     0,     0,     0,     0,
       7,     8,  1562,  1517,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1577,   282,     0,     0,  1581,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   282,
       0,     0,     0,     0,  1972,   282,   282,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   282,     0,     0,   282,     0,     0,     0,     0,
       0,   442,     0,     0,     0,     0,   282,     0,     0,     0,
     282,  1628,  1629,  1630,     0,     0,     0,     0,     0,  1635,
       0,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,  1650,   364,     0,    41,    42,    43,
      44,     0,     0,    46,     0,  1658,   282,     0,  1486,     0,
     740,    53,  1488,  1666,    56,   741,     0,   742,   743,     0,
     744,  1492,     0,  1493,     0,     0,  1495,     0,     0,   436,
       0,     0,   436,   282,     0,     0,  2071,    77,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   746,     0,   746,
       0,     0,  2002,  2003,  2004,  2005,  2006,     0,     0,     0,
       0,    91,    92,    93,   436,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1535,     0,
       0,     0,     0,  1730,   282,  1648,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   436,     0,     0,
       0,     0,   282,   282,   282,   282,   282,     0,     0,  1750,
       0,     0,     0,  2040,     0,     0,     0,  2044,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,   928,     0,     0,     0,
     580,     0,     0,     0,     0,     0,   282,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   282,     0,     0,     0,  2078,     0,     0,
       0,  2080,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   450,  1790,  2092,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   282,
       0,     0,     0,     0,     0,     0,   450,   450,     0,     0,
       0,     0,  1730,     0,     0,   410,     0,     0,     0,     0,
     450,     0,     0,     0,     0,     0,     0,     0,  1826,     0,
       0,     0,     0,     0,  2127,  2128,     0,     0,     0,     0,
       0,     0,   437,     0,     0,     0,   447,     0,     0,     0,
     452,     0,     0,     0,     0,     0,     0,     0,     0,  1850,
       0,     0,     0,     0,     0,     0,     0,     0,  1704,     0,
       0,     0,  1859,   468,   469,   470,     0,   472,   473,   474,
       0,   476,   477,   478,   479,   480,   481,   482,     0,   484,
     485,   486,   487,     0,     0,     0,   491,     0,     0,  1744,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,  1730,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   746,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1918,     0,     0,     0,     0,     0,   548,   550,   552,   553,
     491,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   587,   491,   491,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   603,     0,   746,     0,     0,
       0,     0,     0,     0,     0,   613,   614,     0,     0,     0,
       0,   614,  1960,     0,     0,     0,     0,     0,     0,     0,
     746,   491,   636,     0,     0,  1968,  1835,     0,     0,  1837,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   656,   491,     0,  1730,     0,
     660,   661,   662,     0,   664,     0,     0,     0,   667,   668,
     669,     0,  1997,   670,     0,     0,     0,     0,     0,     0,
       0,     0,   436,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,   436,   436,   580,     0,     0,     0,     0,
       0,   702,  1730,     0,     0,   491,   711,   436,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   719,   720,     0,     0,     0,     0,     0,  2041,
       0,     0,     0,  1708,     0,     0,     0,   747,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   773,   774,   775,   776,     0,
       0,     0,     0,   781,     0,     0,     0,   785,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   799,   550,     0,     0,     0,   803,     0,     0,     0,
     807,   808,   809,   810,   811,   812,   813,   814,   815,   816,
     817,   818,   819,   820,   821,   822,   823,   825,   826,   827,
     828,   829,   830,   831,   831,     0,   836,   837,     0,   839,
       0,     0,     0,     0,     0,     0,   846,     0,     0,   850,
     851,     0,     0,     0,     0,     0,   831,     0,     0,     0,
       0,   491,  1730,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   867,   868,   869,   870,   871,   872,   873,
     874,   875,   876,   877,   878,   879,   880,   882,   884,   885,
     886,   887,   888,     0,   890,   891,     0,     0,     0,     0,
       0,     0,   899,   900,   901,     0,     0,     0,     0,     0,
     907,   908,     0,   491,   491,  1730,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   548,   702,   927,     0,     0,     0,   932,   933,
     934,   935,   936,   937,   938,   939,   940,     0,     0,     0,
       0,     0,   947,     0,   949,     0,   950,     0,     0,     0,
     491,   491,   491,   955,     0,   957,   958,   959,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,   982,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   491,     0,     0,     0,
       0,     0,     0,  1001,     0,  1002,     0,     0,   890,   891,
       0,     0,     0,     0,     0,     0,     0,  1015,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     491,     0,  1031,     0,     0,     0,     0,  1034,  1035,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   491,     0,  1051,  1052,     0,     0,   491,  1057,
       0,     0,     0,  1051,     0,     0,     0,     0,     0,     0,
    1068,   160,   372,     0,  1071,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2150,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1095,     0,   882,     0,  1098,     0,     0,     0,     0,     0,
    2166,  2168,     0,   373,   374,   375,   376,   377,   378,   379,
     380,   381,     0,     0,  2176,     0,   382,   383,   384,   385,
       0,  1714,     0,     0,   386,   387,   388,   491,     0,   389,
       0,   390,     0,     0,     0,     0,     0,   491,     0,     0,
       0,     0,     0,     0,   391,     0,     0,   392,     0,     0,
       0,     0,   492,     0,     0,     0,     0,     0,     0,   493,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   491,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   491,     0,     0,     0,
       0,     0,     0,     0,  1179,     0,     0,     0,  1185,     0,
       0,     0,     0,     0,   491,     0,     0,     0,     0,     0,
       0,     0,  1198,  1199,     0,     0,     0,     0,  1202,     0,
       0,     0,     0,   551,     0,     0,   492,     0,     0,     0,
       0,     0,     0,     0,     0,  1228,     0,     0,  1230,     0,
       0,     0,   492,   492,     0,     0,     0,     0,     0,   588,
     589,     0,     0,     0,     0,   393,   394,   395,  1245,     0,
    1247,  1248,     0,     0,     0,     0,   396,     0,     0,     0,
     397,     0,   398,   114,     0,     0,     0,     0,     0,  1266,
       0,     0,     0,     0,     0,     0,     0,   492,     0,     0,
       0,     0,     0,     0,   635,     0,     0,     0,     0,     0,
       0,     0,  1288,     0,     0,     0,     0,     0,     0,  1294,
       0,     0,   492,     0,  1303,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1329,
    1330,     0,     0,     0,     0,   491,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,     0,     0,  1354,  1355,  1356,  1357,     0,     0,
       0,   492,     0,     0,     0,     0,     0,     0,   710,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   491,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1380,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1381,  1382,  1383,     0,     0,
       0,     0,  1385,  1386,     0,  1387,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1394,  1395,     0,
       0,     0,     0,     0,     0,  1398,  1399,  1824,   551,     0,
       0,     0,  1403,     0,  1404,  1405,     0,     0,     0,     0,
       0,  1411,  1412,     0,     0,     0,   491,   491,     0,     0,
       0,     0,   824,     0,     0,     0,     0,     0,     0,   832,
     833,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   402,   372,     0,     0,     0,  1440,     0,
       0,     0,   856,     0,     0,     0,     0,   492,     0,     0,
       0,     0,     0,     0,   860,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1472,   883,     0,   373,   374,   375,   376,   377,
     378,   379,   380,   381,   491,     0,     0,     0,   382,   383,
     384,   385,  1487,     0,     0,     0,   386,   387,   388,   492,
     492,   389,     0,   390,     0,     0,   910,   911,     0,     0,
       0,     0,     0,  1411,  1412,     0,   391,     0,  1503,   392,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1509,
       0,     0,  1514,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   492,   492,   492,     0,
       0,     0,     0,   952,   953,   954,     0,     0,     0,     0,
    1536,     0,     0,  1539,   491,  1542,   491,     0,     0,     0,
       0,     0,     0,     0,  1551,     0,  1553,     0,  1551,  1557,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   492,     0,     0,     0,     0,     0,     0,   996,
       0,     0,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,   492,   393,   394,   395,
       0,     0,     0,  1027,     0,     0,     0,     0,   396,     0,
       0,     0,   397,     0,   398,     0,     0,     0,   492,     0,
       0,     0,     0,     0,   492,  1049,     0,  1634,     0,     0,
       0,  1056,     0,     0,  1641,  1642,  1643,     0,     0,     0,
       0,     0,     0,  1649,     0,  1651,     0,     0,     0,     0,
       0,     0,  1654,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   491,     0,   883,     0,
    1099,  1671,     0,     0,     0,     0,     0,     0,     0,  1676,
    1677,  1678,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1689,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1701,     0,   492,     0,     0,     0,     0,     0,     0,
    1129,     0,     0,   492,     0,     0,     0,     0,     0,     0,
    1139,     0,     0,     0,  1718,     0,     0,     0,     0,     0,
    1727,  1728,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   492,
       0,     0,     0,     0,     0,     0,  1164,     0,     0,     0,
       0,     0,   492,  1755,     0,     0,     0,     0,     0,  1172,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     492,     0,     0,     0,     0,     0,     0,  1191,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1773,  1774,  1775,     0,     0,     0,     0,
       0,     0,     0,     0,  1776,     0,     0,     0,     0,  1779,
    1780,  1781,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1789,     0,     0,     0,     0,   491,   491,
    1795,     0,  1796,     0,     0,     0,     0,     0,  1804,  1805,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   491,     0,     0,
       0,     0,  1823,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1833,     0,     0,  1836,     0,
       0,     0,  1838,     0,     0,     0,     0,     0,     0,     0,
       0,  1780,  1781,     0,  1847,     0,   729,     0,     0,     0,
       0,   492,  1852,   372,     0,     0,     0,     0,  1336,  1909,
       0,     0,     0,     0,     0,     0,     0,     0,  1865,     0,
       0,     0,     0,  1871,  1872,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1877,     0,     0,
    1878,  1877,     0,  1935,   373,   374,   375,   376,   377,   378,
     379,   380,   381,     0,     0,   492,     0,   382,   383,   384,
     385,     0,  1379,     0,     0,   386,   387,   388,     0,     0,
     389,     0,   390,     0,     0,     0,   491,   491,     0,     0,
       0,     0,     0,     0,     0,   391,     0,     0,   392,     0,
       0,     0,     0,  1912,     0,     0,   847,     0,     0,     0,
       0,     0,  1921,     0,     0,     0,     0,     0,     0,     0,
       0,   550,     0,     0,     0,     0,     0,     0,     0,     0,
    1932,  1933,     0,     0,   491,     0,     0,     0,     0,     0,
       0,     0,   492,   492,     0,     0,     0,     0,     0,  1416,
    1417,     0,     0,     0,     0,     0,     0,   895,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1957,  1958,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1970,     0,     0,
     491,   491,  1975,     0,     0,     0,     0,  1978,     0,     0,
    1980,     0,     0,     0,     0,     0,   393,   394,   395,     0,
     491,     0,     0,     0,     0,     0,     0,   396,     0,     0,
     492,   397,     0,   398,     0,   491,     0,  1485,     0,     0,
       0,     0,     0,  2007,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  2018,     0,
       0,     0,     0,     0,  2023,     0,     0,  2024,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,     0,     0,  2048,     0,     0,   165,
     492,   170,   492,     0,     0,     0,     0,  1541,     0,  1543,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2066,  2067,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2073,  2074,     0,  2075,   491,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1947,     0,  2084,     0,     0,     0,     0,   491,   491,   491,
     491,   491,  1421,     0,     0,     0,     0,  1871,     0,     0,
       0,  2100,     0,     0,     0,     0,     0,     0,  1102,  1103,
    1104,  1105,  1106,  1107,  1108,  1109,  1110,  1111,  1112,     0,
    1114,  1115,  1116,  1117,     0,  1119,  1120,  1121,  1122,  2121,
    2060,     0,     0,     0,  2124,  2125,  2126,     0,     0,  1132,
       0,  1134,     0,     0,     0,   435,     0,     0,   491,  1141,
       0,     0,     0,     0,     0,     0,  1147,  1148,     0,     0,
       0,     0,   492,     0,     0,  1159,     0,     0,     0,  1669,
       0,     0,     0,     0,     0,     0,  2154,  2155,  2156,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   498,   499,   500,   501,   502,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,   517,   518,   519,   520,   523,   524,   525,   526,   527,
     528,     0,   530,   531,     0,     0,     0,     0,     0,   534,
     535,     0,     0,     0,     0,     0,     0,     0,     0,  1756,
       0,   543,   544,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1251,  1253,  1254,     0,     0,     0,
    1258,  1259,     0,     0,  1262,  1263,  1264,  1265,     0,  1267,
       0,     0,     0,     0,  1273,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,  2108,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   492,   492,     0,     0,     0,     0,
       0,  1793,  1794,     0,     0,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,   492,     0,     0,     0,   580,     0,     0,
       0,     0,  1353,     0,     0,     0,     0,  1358,   671,   672,
     673,     0,     0,     0,   677,   678,   679,   680,   681,   682,
     683,     0,   684,     0,     0,     0,   685,   686,     0,     0,
     689,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   701,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   202,   160,     0,     0,
       0,     0,     0,   203,   204,   205,     0,     0,   206,   207,
     208,   209,   210,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   408,   224,   225,   226,   227,
     228,   229,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1413,     0,     0,     0,     0,     0,   235,   236,
     237,   238,   492,   492,     0,     0,   239,     0,     0,  1901,
    1902,     0,     0,     0,     0,     0,     0,     0,   241,   242,
     243,     0,  1886,   561,     0,     0,     0,     0,   244,    23,
       0,   245,     0,     0,     0,     0,     0,   551,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     492,     0,     0,     0,     0,     0,     0,  1934,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,   864,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,  2028,     0,     0,     0,   580,   492,   492,     0,     0,
       0,  1800,     0,  1973,  1974,  1801,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   492,     0,     0,     0,
    1523,     0,     0,  1991,   561,     0,     0,     0,     0,     0,
     255,   492,     0,     0,   256,     0,     0,     0,  2001,   258,
     259,   260,     0,   261,   262,   263,     0,     0,     0,   114,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1559,
       0,     0,   561,     0,     0,   264,   409,     0,     0,     0,
       0,     0,   266,     0,     0,     0,     0,   353,     0,     0,
       0,   602,     0,     0,   269,     0,     0,     0,     0,     0,
       0,  1596,  1597,  1598,  1599,  1600,     0,     0,     0,     0,
    1605,  1606,     0,     0,  1608,     0,  1610,     0,     0,     0,
    1614,     0,     0,  1616,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1624,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   800,
     576,   577,   578,   579,     0,   492,     0,     0,   580,     0,
       0,     0,  2076,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   492,   492,   492,   492,   492,     0,     0,
    2088,  2089,  2090,  2091,  2093,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,  1684,   580,  1687,     0,
       0,     0,  1694,   797,  1696,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,   492,     0,     0,   580,     0,     0,
       0,  2134,     0,     0,     0,  1721,     0,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1154,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   800,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1171,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     202,     6,   372,     0,     0,     0,     0,   203,   204,   205,
       0,  1782,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   969,     0,  1560,   374,   375,   376,   377,   378,   379,
     380,   381,   235,   236,   237,   238,   382,   383,   384,   385,
     239,   240,     0,     0,   386,   387,   388,     0,     0,   389,
       0,   390,   241,   242,   243,     0,     0,     0,     0,     0,
       0,     0,   244,    23,   391,   245,     0,   392,     0,     0,
       0,   246,     0,     0,   247,     0,     0,   248,     0,   249,
       0,     0,     0,     0,     0,     0,    40,     0,     0,     0,
     250,     0,     0,  1862,     0,     0,  1863,     0,     0,     0,
     251,     0,    54,    55,     0,   252,     0,   253,     0,     0,
     254,     0,     0,     0,     0,    65,    66,    67,    68,    69,
       0,    71,    72,    73,    74,    75,    76,     0,     0,     0,
      80,     0,     0,    83,     0,     0,     0,     0,     0,     0,
    1348,     0,     0,     0,  1893,  1894,     0,  1896,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,   960,     0,   255,   393,   394,   395,   256,   257,
     920,     0,     0,   258,   259,   260,   396,   261,   262,   263,
     397,     0,   398,   114,     0,     0,     0,     0,  1162,     0,
       0,     0,     0,     0,     0,     0,     0,  1938,  1939,   264,
     265,     0,     0,     0,     0,     0,   266,     0,     0,     0,
       0,   353,     0,     0,     0,  1561,     0,     0,   269,     0,
       0,     0,     0,  1402,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,  1969,     0,   580,   202,     6,   372,
       0,     0,     0,     0,   203,   204,   205,     0,     0,   206,
     207,   208,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,     0,
     373,   374,   375,   376,   377,   378,   379,   380,   381,   235,
     236,   237,   238,   382,   383,   384,   385,   239,   240,     0,
       0,   386,   387,   388,     0,     0,   389,     0,   390,   241,
     242,   243,     0,     0,     0,     0,     0,     0,     0,   244,
      23,   391,   245,     0,   392,     0,     0,     0,   246,     0,
       0,   247,     0,     0,   248,     0,   249,     0,     0,     0,
       0,     0,     0,    40,     0,     0,     0,   250,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   251,     0,    54,
      55,     0,   252,     0,   253,  1524,  1525,   254,     0,     0,
       0,     0,    65,    66,    67,    68,    69,     0,    71,    72,
      73,    74,    75,    76,     0,     0,     0,    80,     0,     0,
      83,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,  1163,
       0,     0,     0,   580,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   255,   393,   394,   395,   256,   257,     0,     0,     0,
     258,   259,   260,   396,   261,   262,   263,   397,     0,   398,
     114,     0,  1612,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1620,     0,     0,     0,   264,   265,     0,     0,
       0,     0,     0,   266,     0,     0,  2143,     0,   353,   202,
     160,   372,   268,   421,     0,   269,   203,   204,   205,     0,
       0,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   408,   224,
     225,   226,   227,   228,   229,     0,     0,     0,     0,     0,
       0,     0,   373,   374,   375,   376,   377,   378,   379,   380,
     381,   235,   236,   237,   238,   382,   383,   384,   385,   239,
       0,     0,     0,   386,   387,   388,     0,     0,   389,     0,
     390,   241,   242,   243,     0,     0,     0,     0,     0,     0,
       0,   244,    23,   391,   245,     0,   392,     0,     0,     0,
     292,     0,     0,   293,     0,     0,   294,     0,   295,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    44,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   422,     0,     0,     0,
       0,     0,     0,     0,     0,  1753,     0,     0,     0,     0,
       0,     0,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,  1478,     0,   423,     0,   580,     0,
     726,     0,   727,   255,   393,   394,   395,   256,  1686,     0,
       0,     0,   258,   259,   260,   396,   261,   262,   263,   397,
       0,   398,   114,   424,   425,   426,   427,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   264,   409,
       0,     0,     0,     0,     0,   266,   428,   429,     0,     0,
     430,     0,   431,   202,     6,   352,   432,   269,     0,     0,
     203,   204,   205,     0,     0,   206,   207,   208,   209,   210,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   235,   236,   237,   238,     0,
       0,     0,     0,   239,   240,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   241,   242,   243,     0,     0,
       0,     0,     0,     0,     0,   244,    23,     0,   245,     0,
       0,     0,     0,     0,   246,     0,     0,   247,     0,     0,
     248,     0,   249,     0,     0,     0,     0,     0,     0,    40,
       0,     0,     0,   250,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   251,     0,    54,    55,     0,   252,     0,
     253,     0,     0,   254,     0,     0,     0,     0,    65,    66,
      67,    68,    69,     0,    71,    72,    73,    74,    75,    76,
       0,     0,     0,    80,     0,     0,    83,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   255,     0,     0,
       0,   256,   257,     0,     0,     0,   258,   259,   260,     0,
     261,   262,   263,     0,     0,     0,   114,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   264,   265,     0,     0,     0,     0,     0,   266,
       0,     0,     0,     0,   353,   202,     6,     0,   354,     0,
     657,   269,   203,   204,   205,     0,     0,   206,   207,   208,
     209,   210,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   235,   236,   237,
     238,     0,     0,     0,     0,   239,   240,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   241,   242,   243,
       0,     0,     0,     0,     0,     0,     0,   244,    23,     0,
     245,     0,     0,     0,     0,     0,   246,     0,     0,   247,
       0,     0,   248,     0,   249,     0,     0,     0,     0,     0,
       0,    40,     0,     0,     0,   250,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   251,     0,    54,    55,     0,
     252,     0,   253,     0,     0,   254,     0,     0,     0,     0,
      65,    66,    67,    68,    69,     0,    71,    72,    73,    74,
      75,    76,     0,     0,     0,    80,     0,     0,    83,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,  1827,     0,     0,
       0,  1828,     0,     0,     0,     0,     0,     0,     0,   255,
       0,     0,     0,   256,   257,     0,     0,     0,   258,   259,
     260,     0,   261,   262,   263,     0,     0,     0,   114,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   264,   490,     0,     0,     0,     0,
       0,   266,     0,     0,     0,     0,   353,   202,     6,  1657,
       0,   554,     0,   269,   203,   204,   205,     0,     0,   206,
     207,   208,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   235,
     236,   237,   238,     0,     0,     0,     0,   239,   240,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   241,
     242,   243,     0,     0,     0,     0,     0,     0,     0,   244,
      23,     0,   245,     0,     0,     0,     0,     0,   246,     0,
       0,   247,     0,     0,   248,     0,   249,     0,     0,     0,
       0,     0,     0,    40,     0,     0,     0,   250,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   251,     0,    54,
      55,     0,   252,     0,   253,     0,     0,   254,     0,     0,
       0,     0,    65,    66,    67,    68,    69,     0,    71,    72,
      73,    74,    75,    76,     0,     0,     0,    80,     0,     0,
      83,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,     0,  1993,
       0,     0,     0,  1994,     0,     0,     0,     0,     0,     0,
       0,   255,     0,     0,     0,   256,   257,     0,     0,     0,
     258,   259,   260,     0,   261,   262,   263,     0,     0,     0,
     114,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   264,   265,     0,     0,
       0,     0,     0,   266,     0,     0,     0,     0,   353,   202,
       6,     0,   268,     0,     0,   269,   203,   204,   205,     0,
       0,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   235,   236,   237,   238,     0,     0,     0,     0,   239,
     240,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   241,   242,   243,     0,     0,     0,     0,     0,     0,
       0,   244,    23,     0,   245,     0,     0,     0,     0,     0,
     246,     0,     0,   247,     0,     0,   248,     0,   249,     0,
       0,     0,     0,     0,     0,    40,     0,     0,     0,   250,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   251,
       0,    54,    55,     0,   252,     0,   253,     0,     0,   254,
       0,     0,     0,     0,    65,    66,    67,    68,    69,     0,
      71,    72,    73,    74,    75,    76,     0,     0,     0,    80,
       0,     0,    83,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   757,     0,     0,     0,     0,
       0,     0,     0,   255,     0,     0,     0,   256,   257,     0,
       0,     0,   258,   259,   260,     0,   261,   262,   263,     0,
       0,     0,   114,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   620,  1819,
       0,     0,     0,     0,     0,   266,     0,     0,     0,     0,
     622,   202,     6,     0,   314,   554,     0,   269,   203,   204,
     205,     0,     0,   206,   207,   208,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   220,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   235,   236,   237,   238,     0,     0,     0,
       0,   239,   240,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   241,   242,   243,     0,     0,     0,     0,
       0,     0,     0,   244,    23,     0,   245,     0,     0,     0,
       0,     0,   246,     0,     0,   247,     0,     0,   248,     0,
     249,     0,     0,     0,     0,     0,     0,    40,     0,     0,
       0,   250,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   251,     0,    54,    55,     0,   252,     0,   253,     0,
       0,   254,     0,     0,     0,     0,    65,    66,    67,    68,
      69,     0,    71,    72,    73,    74,    75,    76,     0,     0,
       0,    80,     0,     0,    83,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   920,     0,     0,
       0,     0,     0,     0,     0,   255,     0,     0,     0,   256,
     257,     0,     0,     0,   258,   259,   260,     0,   261,   262,
     263,     0,     0,     0,   114,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     264,   265,     0,     0,     0,     0,     0,   266,     0,     0,
       0,     0,   267,   202,     6,     0,   268,     0,     0,   269,
     203,   204,   205,     0,     0,   206,   207,   208,   209,   210,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   235,   236,   237,   238,     0,
       0,     0,     0,   239,   240,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   241,   242,   243,     0,     0,
       0,     0,     0,     0,     0,   244,    23,     0,   245,     0,
       0,     0,     0,     0,   246,     0,     0,   247,     0,     0,
     248,     0,   249,     0,     0,     0,     0,     0,     0,    40,
       0,     0,     0,   250,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   251,     0,    54,    55,     0,   252,     0,
     253,     0,     0,   254,     0,     0,     0,     0,    65,    66,
      67,    68,    69,     0,    71,    72,    73,    74,    75,    76,
       0,     0,     0,    80,     0,     0,    83,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   921,
       0,     0,     0,     0,     0,     0,     0,   255,     0,     0,
       0,   256,   257,     0,     0,     0,   258,   259,   260,     0,
     261,   262,   263,     0,     0,     0,   114,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   264,   265,     0,     0,     0,     0,     0,   266,
       0,     0,     0,     0,   353,   202,     6,     0,   268,     0,
       0,   269,   203,   204,   205,     0,     0,   206,   207,   208,
     209,   210,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   235,   236,   237,
     238,     0,     0,     0,     0,   239,   240,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   241,   242,   243,
       0,     0,     0,     0,     0,     0,     0,   244,    23,     0,
     245,     0,     0,     0,     0,     0,   246,     0,     0,   247,
       0,     0,   248,     0,   249,     0,     0,     0,     0,     0,
       0,    40,     0,     0,     0,   250,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   251,     0,    54,    55,     0,
     252,     0,   253,     0,     0,   254,     0,     0,     0,     0,
      65,    66,    67,    68,    69,     0,    71,    72,    73,    74,
      75,    76,     0,     0,     0,    80,     0,     0,    83,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   966,     0,     0,     0,     0,     0,     0,     0,   255,
       0,     0,     0,   256,   257,     0,     0,     0,   258,   259,
     260,     0,   261,   262,   263,     0,     0,     0,   114,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   264,   490,     0,     0,     0,     0,
       0,   266,     0,     0,     0,     0,   353,   202,     6,     0,
     549,     0,     0,   269,   203,   204,   205,     0,     0,   206,
     207,   208,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   235,
     236,   237,   238,     0,     0,     0,     0,   239,   240,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   241,
     242,   243,     0,     0,     0,     0,     0,     0,     0,   244,
      23,     0,   245,     0,     0,     0,     0,     0,   246,     0,
       0,   247,     0,     0,   248,     0,   249,     0,     0,     0,
       0,     0,     0,    40,     0,     0,     0,   250,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   251,     0,    54,
      55,     0,   252,     0,   253,     0,     0,   254,     0,     0,
       0,     0,    65,    66,    67,    68,    69,     0,    71,    72,
      73,    74,    75,    76,     0,     0,     0,    80,     0,     0,
      83,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1072,     0,     0,     0,     0,     0,     0,
       0,   255,     0,     0,     0,   256,   257,     0,     0,     0,
     258,   259,   260,     0,   261,   262,   263,     0,     0,     0,
     114,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   264,   490,     0,     0,
       0,     0,     0,   266,     0,     0,     0,     0,   353,   202,
       6,     0,     0,   554,     0,   269,   203,   204,   205,     0,
       0,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   235,   236,   237,   238,     0,     0,     0,     0,   239,
     240,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   241,   242,   243,     0,     0,     0,     0,     0,     0,
       0,   244,    23,     0,   245,     0,     0,     0,     0,     0,
     246,     0,     0,   247,     0,     0,   248,     0,   249,     0,
       0,     0,     0,     0,     0,    40,     0,     0,     0,   250,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   251,
       0,    54,    55,     0,   252,     0,   253,     0,     0,   254,
       0,     0,     0,     0,    65,    66,    67,    68,    69,     0,
      71,    72,    73,    74,    75,    76,     0,     0,     0,    80,
       0,     0,    83,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1073,     0,     0,     0,     0,
       0,     0,     0,   255,     0,     0,     0,   256,   257,     0,
       0,     0,   258,   259,   260,     0,   261,   262,   263,     0,
       0,     0,   114,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   264,   265,
       0,     0,     0,     0,     0,   266,     0,     0,     0,     0,
     592,   202,     6,     0,   268,     0,     0,   269,   203,   204,
     205,     0,     0,   206,   207,   208,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   220,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   235,   236,   237,   238,     0,     0,     0,
       0,   239,   240,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   241,   242,   243,     0,     0,     0,     0,
       0,     0,     0,   244,    23,     0,   245,     0,     0,     0,
       0,     0,   246,     0,     0,   247,     0,     0,   248,     0,
     249,     0,     0,     0,     0,     0,     0,    40,     0,     0,
       0,   250,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   251,     0,    54,    55,     0,   252,     0,   253,     0,
       0,   254,     0,     0,     0,     0,    65,    66,    67,    68,
      69,     0,    71,    72,    73,    74,    75,    76,     0,     0,
       0,    80,     0,     0,    83,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1074,     0,     0,
       0,     0,     0,     0,     0,   255,     0,     0,     0,   256,
     257,     0,     0,     0,   258,   259,   260,     0,   261,   262,
     263,     0,     0,     0,   114,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     264,   265,     0,     0,     0,     0,     0,   266,     0,     0,
       0,     0,   599,   202,     6,     0,   268,     0,     0,   269,
     203,   204,   205,     0,     0,   206,   207,   208,   209,   210,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   235,   236,   237,   238,     0,
       0,     0,     0,   239,   240,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   241,   242,   243,     0,     0,
       0,     0,     0,     0,     0,   244,    23,     0,   245,     0,
       0,     0,     0,     0,   246,     0,     0,   247,     0,     0,
     248,     0,   249,     0,     0,     0,     0,     0,     0,    40,
       0,     0,     0,   250,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   251,     0,    54,    55,     0,   252,     0,
     253,     0,     0,   254,     0,     0,     0,     0,    65,    66,
      67,    68,    69,     0,    71,    72,    73,    74,    75,    76,
       0,     0,     0,    80,     0,     0,    83,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1075,
       0,     0,     0,     0,     0,     0,     0,   255,     0,     0,
       0,   256,   257,     0,     0,     0,   258,   259,   260,     0,
     261,   262,   263,     0,     0,     0,   114,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   264,   490,     0,     0,     0,     0,     0,   266,
       0,     0,     0,     0,   353,   202,     6,     0,   881,     0,
    1335,   269,   203,   204,   205,     0,     0,   206,   207,   208,
     209,   210,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   235,   236,   237,
     238,     0,     0,     0,     0,   239,   240,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   241,   242,   243,
       0,     0,     0,     0,     0,     0,     0,   244,    23,     0,
     245,     0,     0,     0,     0,     0,   246,     0,     0,   247,
       0,     0,   248,     0,   249,     0,     0,     0,     0,     0,
       0,    40,     0,     0,     0,   250,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   251,     0,    54,    55,     0,
     252,     0,   253,     0,     0,   254,     0,     0,     0,     0,
      65,    66,    67,    68,    69,     0,    71,    72,    73,    74,
      75,    76,     0,     0,     0,    80,     0,     0,    83,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1113,     0,     0,     0,     0,     0,     0,     0,   255,
       0,     0,     0,   256,   257,     0,     0,     0,   258,   259,
     260,     0,   261,   262,   263,     0,     0,     0,   114,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   264,   490,     0,     0,     0,     0,
       0,   266,   202,     6,     0,     0,   353,  1540,     0,   203,
     204,   205,     0,   269,   206,   207,   208,   209,   210,   211,
     212,   213,   214,   215,   216,   217,   218,   219,   220,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   235,   236,   237,   238,     0,     0,
       0,     0,   239,   240,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   241,   242,   243,     0,     0,     0,
       0,     0,     0,     0,   244,    23,     0,   245,     0,     0,
       0,     0,     0,   246,     0,     0,   247,     0,     0,   248,
       0,   249,     0,     0,     0,     0,     0,     0,    40,     0,
       0,     0,   250,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   251,     0,    54,    55,     0,   252,     0,   253,
       0,     0,   254,     0,     0,     0,     0,    65,    66,    67,
      68,    69,     0,    71,    72,    73,    74,    75,    76,     0,
       0,     0,    80,     0,     0,    83,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1123,     0,
       0,     0,     0,     0,     0,     0,   255,     0,     0,     0,
     256,   257,     0,     0,     0,   258,   259,   260,     0,   261,
     262,   263,     0,     0,     0,   114,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   264,   490,     0,     0,     0,     0,     0,   266,   202,
       6,     0,     0,   353,     0,     0,   203,   204,   205,     0,
     269,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   235,   236,   237,   238,     0,     0,     0,     0,   239,
     240,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   241,   242,   243,     0,     0,     0,     0,     0,     0,
       0,   244,    23,     0,   245,     0,     0,     0,     0,     0,
     246,     0,     0,   247,     0,     0,   248,     0,   249,     0,
       0,     0,     0,     0,     0,    40,     0,     0,     0,   250,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   251,
       0,    54,    55,     0,   252,     0,   253,     0,     0,   254,
       0,     0,     0,     0,    65,    66,    67,    68,    69,     0,
      71,    72,    73,    74,    75,    76,     0,     0,     0,    80,
       0,     0,    83,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1124,     0,     0,     0,     0,
       0,     0,     0,   255,     0,     0,     0,   256,   257,     0,
       0,     0,   258,   259,   260,     0,   261,   262,   263,     0,
       0,     0,   114,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   264,   265,
       0,     0,     0,     0,     0,   266,     0,     0,     0,     0,
     353,   202,     6,     0,  1665,     0,     0,   269,   203,   204,
     205,     0,     0,   206,   207,   208,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   220,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   235,   236,   237,   238,     0,     0,     0,
       0,   239,   240,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   241,   242,   243,     0,     0,     0,     0,
       0,     0,     0,   244,    23,     0,   245,     0,     0,     0,
       0,     0,   246,     0,     0,   247,     0,     0,   248,     0,
     249,     0,     0,     0,     0,     0,     0,    40,     0,     0,
       0,   250,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   251,     0,    54,    55,     0,   252,     0,   253,     0,
       0,   254,     0,     0,     0,     0,    65,    66,    67,    68,
      69,     0,    71,    72,    73,    74,    75,    76,     0,     0,
       0,    80,     0,     0,    83,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1125,     0,     0,
       0,     0,     0,     0,     0,   255,     0,     0,     0,   256,
     257,     0,     0,     0,   258,   259,   260,     0,   261,   262,
     263,     0,     0,     0,   114,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     620,  1819,     0,     0,     0,     0,     0,   266,     0,     0,
       0,     0,   622,   202,     6,     0,   314,     0,     0,   269,
     203,   204,   205,     0,     0,   206,   207,   208,   209,   210,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   235,   236,   237,   238,     0,
       0,     0,     0,   239,   240,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   241,   242,   243,     0,     0,
       0,     0,     0,     0,     0,   244,    23,     0,   245,     0,
       0,     0,     0,     0,   246,     0,     0,   247,     0,     0,
     248,     0,   249,     0,     0,     0,     0,     0,     0,    40,
       0,     0,     0,   250,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   251,     0,    54,    55,     0,   252,     0,
     253,     0,     0,   254,     0,     0,     0,     0,    65,    66,
      67,    68,    69,     0,    71,    72,    73,    74,    75,    76,
       0,     0,     0,    80,     0,     0,    83,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1130,
       0,     0,     0,     0,     0,     0,     0,   255,     0,     0,
       0,   256,   257,     0,     0,     0,   258,   259,   260,     0,
     261,   262,   263,     0,     0,     0,   114,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   264,   490,     0,     0,     0,     0,     0,   266,
     202,   160,   372,     0,   353,     0,     0,   203,   204,   205,
       0,   269,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   408,
     224,   225,   226,   227,   228,   229,     0,     0,     0,     0,
       0,     0,     0,   373,   374,   375,   376,   377,   378,   379,
     380,   381,   235,   236,   237,   238,   382,   383,   384,   385,
     239,     0,     0,     0,   386,   387,   388,     0,     0,   389,
       0,   390,   241,   242,   243,     0,     0,     0,     0,     0,
       0,     0,   244,    23,   391,   245,     0,   392,   202,   160,
     372,     0,     0,     0,     0,   203,   204,   205,     0,     0,
     206,   207,   208,   209,   210,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   221,   222,   408,   224,   225,
     226,   227,   228,   229,     0,     0,     0,     0,     0,     0,
       0,   373,   374,   375,   376,   377,   378,   379,   380,   381,
     235,   236,   237,   238,   382,   383,   384,   385,   239,     0,
       0,     0,   386,   387,   388,     0,     0,   389,     0,   390,
     241,   242,   243,     0,     0,     0,     0,     0,     0,     0,
     244,    23,   391,   245,     0,   392,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   255,   393,   394,   395,   256,     0,
       0,     0,     0,   258,   259,   260,   396,   261,   262,   263,
     397,     0,   398,   114,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   264,
     409,     0,     0,     0,     0,     0,   266,     0,     0,     0,
       0,   353,     0,     0,     0,  1556,     0,     0,   269,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1131,   255,   393,   394,   395,   256,     0,     0,     0,
       0,   258,   259,   260,   396,   261,   262,   263,   397,     0,
     398,   114,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   264,   409,     0,
       0,     0,     0,     0,   266,     0,     0,     0,     0,   353,
     202,   160,   372,  1832,     0,     0,   269,   203,   204,   205,
       0,     0,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   408,
     224,   225,   226,   227,   228,   229,     0,     0,     0,     0,
       0,     0,     0,  1560,   374,   375,   376,   377,   378,   379,
     380,   381,   235,   236,   237,   238,   382,   383,   384,   385,
     239,     0,     0,     0,   386,   387,   388,     0,     0,   389,
       0,   390,   241,   242,   243,     0,     0,     0,     0,     0,
       0,     0,   244,    23,   391,   245,     0,   392,   202,   160,
     372,     0,     0,     0,     0,   203,   204,   205,     0,     0,
     206,   207,   208,   209,   210,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   221,   222,   408,   224,   225,
     226,   227,   228,   229,     0,     0,     0,     0,     0,     0,
       0,   373,   374,   375,   376,   377,   378,   379,   380,   381,
     235,   236,   237,   238,   382,   383,   384,   385,   239,     0,
       0,     0,   386,   387,   388,     0,     0,   389,     0,   390,
     241,   242,   243,     0,     0,     0,     0,     0,     0,     0,
     244,    23,   391,   245,     0,   392,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   255,   393,   394,   395,   256,     0,
       0,     0,     0,   258,   259,   260,   396,   261,   262,   263,
     397,     0,   398,   114,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   264,
     409,     0,     0,     0,     0,     0,   266,     0,     0,     0,
       0,   353,     0,     0,     0,  1720,     0,     0,   269,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1133,   255,   393,   394,   395,   256,     0,     0,     0,
       0,   258,   259,   260,   396,   261,   262,   263,   397,     0,
     398,   114,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   264,   409,     0,
       0,     0,     0,     0,   266,   202,   160,     0,     0,   353,
       0,     0,   203,   204,   205,     0,   269,   206,   207,   208,
     209,   210,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   408,   224,   225,   226,   227,   228,
     229,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   235,   236,   237,
     238,     0,     0,     0,     0,   239,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   241,   242,   243,
       0,     0,     0,     0,   202,   160,     0,   244,    23,   849,
     245,   203,   204,   205,     0,     0,   206,   207,   208,   209,
     210,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   408,   224,   225,   226,   227,   228,   229,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   235,   236,   237,   238,
       0,     0,     0,     0,   239,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   241,   242,   243,     0,
       0,     0,     0,     0,     0,     0,   244,    23,     0,   245,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,     0,     0,     0,     0,   255,
       0,     0,  1240,   256,     0,     0,     0,     0,   258,   259,
     260,     0,   261,   262,   263,     0,     0,     0,   114,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   620,   621,     0,     0,     0,     0,
       0,   266,     0,     0,     0,     0,   622,     0,     0,     0,
     314,     0,     0,   269,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   255,     0,
       0,     0,   256,     0,     0,     0,  1261,   258,   259,   260,
       0,   261,   262,   263,     0,     0,     0,   114,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   264,   409,     0,     0,     0,     0,     0,
     266,   202,   160,     0,     0,   353,     0,     0,   203,   204,
     205,     0,   269,   206,   207,   208,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   220,   221,   222,
     408,   224,   225,   226,   227,   228,   229,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   235,   236,   237,   238,     0,     0,     0,
       0,   239,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   241,   242,   243,     0,     0,     0,     0,
       0,     0,     0,   244,    23,     0,   245,   202,   160,     0,
    1287,     0,     0,     0,   203,   204,   205,     0,     0,   206,
     207,   208,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   408,   224,   225,   226,
     227,   228,   229,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     7,     8,     0,   235,
     236,   237,   238,     0,     0,     0,     0,   239,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   241,
     242,   243,     0,     0,     0,     0,     0,     0,     0,   244,
      23,     0,   245,     0,     0,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,   255,     0,   580,     0,   256,
       0,     0,     0,     0,   258,   259,   260,  1344,   261,   262,
     263,     0,     0,     0,   114,     0,     0,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
     264,   409,     0,    41,    42,    43,    44,   266,     0,    46,
       0,     0,   353,  1062,     0,     0,   740,    53,     0,   269,
      56,   741,     0,   742,   743,     0,   744,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    77,     0,     0,     0,     0,     0,     0,
       0,   255,     0,     0,     0,   256,     0,     0,     0,     0,
     258,   259,   260,     0,   261,   262,   263,    91,    92,    93,
     114,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   264,   409,     0,     0,
       0,     0,     0,   266,   202,   160,     0,  1439,   353,     0,
       0,   203,   204,   205,     0,   269,   206,   207,   208,   209,
     210,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   408,   224,   225,   226,   227,   228,   229,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   951,     0,     0,     0,   235,   236,   237,   238,
       0,     0,     0,     0,   239,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   241,   242,   243,     0,
       0,     0,     0,     0,     0,     0,   244,    23,     0,   245,
     202,   160,     0,  1502,     0,     0,     0,   203,   204,   205,
       0,     0,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   408,
     224,   225,   226,   227,   228,   229,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     7,
       8,     0,   235,   236,   237,   238,     0,     0,     0,     0,
     239,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   241,   242,   243,     0,     0,     0,     0,     0,
       0,     0,   244,    23,     0,   245,     0,     0,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,   255,     0,
     580,     0,   256,     0,     0,     0,     0,   258,   259,   260,
    1433,   261,   262,   263,     0,     0,     0,   114,     0,     0,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,   264,   409,     0,    41,    42,    43,    44,
     266,     0,    46,     0,     0,   353,     0,     0,     0,   740,
      53,     0,   269,    56,   741,     0,   742,   743,     0,   744,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    77,     0,     0,     0,
       0,     0,     0,     0,   255,     0,     0,     0,   256,     0,
       0,     0,     0,   258,   259,   260,     0,   261,   262,   263,
      91,    92,    93,   114,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   264,
     409,     0,     0,     0,     0,     0,   266,   202,   160,     0,
       0,   353,     0,     0,   203,   204,   205,     0,   269,   206,
     207,   208,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   408,   224,   225,   226,
     227,   228,   229,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   998,     0,     0,     0,   235,
     236,   237,   238,     0,     0,     0,     0,   239,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   241,
     242,   243,     0,     0,     0,     0,   202,   160,     0,   244,
      23,     0,   245,   203,   204,   205,     0,     0,   206,   207,
     208,   209,   210,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   408,   224,   225,   226,   227,
     228,   229,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   235,   236,
     237,   238,     7,     8,     0,     0,   239,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   241,   242,
     243,     0,     0,     0,     0,     0,     0,     0,   244,    23,
       0,   245,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,     0,     0,
       0,   255,     0,     0,  1434,   256,     0,     0,     0,     0,
     258,   259,   260,     0,   261,   262,   263,     0,     0,     0,
     114,     0,     0,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,   264,   409,     0,    41,
      42,    43,    44,   266,     0,    46,     0,     0,   353,  1555,
       0,     0,   740,    53,     0,   269,    56,   741,     0,   742,
     743,     0,   744,     0,     0,     0,  1009,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    77,
     255,     0,     0,  1010,   256,     0,     0,     0,     0,   258,
     259,   260,     0,   261,   262,   263,     0,     0,     0,   114,
       0,     0,     0,    91,    92,    93,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   264,   409,     0,     0,     0,
       0,     0,   266,   202,   160,     0,     0,   353,  1879,     0,
     203,   204,   205,     0,   269,   206,   207,   208,   209,   210,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   408,   224,   225,   226,   227,   228,   229,     0,
       0,     0,     0,     0,     0,     0,     0,   160,   372,     0,
       0,     0,     0,     0,     0,   235,   236,   237,   238,     0,
       0,     0,     0,   239,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   241,   242,   243,     0,     0,
       0,     0,     0,     0,     0,   244,    23,     0,   245,   373,
     374,   375,   376,   377,   378,   379,   380,   381,     0,     0,
       0,     0,   382,   383,   384,   385,     0,     0,     0,     0,
     386,   387,   388,     0,     0,   389,     0,   390,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     391,     0,     0,   392,     0,     0,     0,   292,     0,     0,
     293,     0,     0,   294,     0,   295,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    44,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   422,     0,     0,     0,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,   255,     0,   580,
       0,   256,     0,     0,     0,     0,   258,   259,   260,  1477,
     261,   262,   263,     0,     0,     0,   114,     0,     0,     0,
       0,     0,   160,   372,     0,     0,     0,     0,     0,     0,
       0,     0,   264,   409,     0,     0,     0,     0,     0,   266,
       0,   393,   394,   395,   353,     0,     0,     0,     0,     0,
       0,   269,   396,     0,     0,     0,   397,     0,   398,   114,
       0,     0,     0,     0,   373,   374,   375,   376,   377,   378,
     379,   380,   381,     0,     0,     0,     0,   382,   383,   384,
     385,     0,     0,   428,   429,   386,   387,   388,     0,     0,
     389,  -551,   390,  1534,   160,   372,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   391,     0,     0,   392,     0,
       0,     0,   292,     0,     0,   293,     0,     0,   294,     0,
     295,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    44,     0,     0,     0,     0,   373,   374,   375,   376,
     377,   378,   379,   380,   381,     0,     0,     0,   422,   382,
     383,   384,   385,     0,     0,     0,     0,   386,   387,   388,
       0,     0,   389,     0,   390,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   391,     0,     0,
     392,     0,     0,     0,   292,     0,     0,   293,     0,     0,
     294,     0,   295,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    44,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   393,   394,   395,     0,
     422,     0,     0,     0,     0,     0,     0,   396,     0,     7,
       8,   397,     0,   398,   114,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,   726,     0,   727,  1720,  1601,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     7,     8,   393,   394,
     395,     0,     0,     0,     0,     0,     0,     0,     0,   396,
       0,     0,     0,   397,     0,   398,   114,     0,     0,     0,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,     0,     0,     0,    41,    42,    43,    44,
       0,     0,    46,     0,     0,     0,     0,     0,  1720,   740,
      53,     0,     0,    56,   741,     0,   742,   743,     0,   744,
       0,     0,     0,     7,     8,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    77,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
       0,     0,     0,    41,    42,    43,    44,     0,     0,    46,
      91,    92,    93,     0,     0,     0,   740,    53,     0,     0,
      56,   741,     0,   742,   743,     0,   744,     0,     0,     0,
       7,     8,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    77,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    91,    92,    93,
      41,    42,    43,    44,     0,     0,    46,     0,     0,     0,
       0,     0,     0,   740,    53,  1000,     0,    56,   741,     0,
     742,   743,     0,   744,     0,     0,     0,     7,     8,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      77,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,     0,     0,     0,    41,    42,    43,
      44,     0,     0,    46,    91,    92,    93,     0,     0,     0,
     740,    53,  1193,     0,    56,   741,     0,   742,   743,     0,
     744,     0,     0,     0,     7,     8,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    77,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    91,    92,    93,    41,    42,    43,    44,     0,     0,
      46,     0,     0,     0,     0,     0,     0,   740,    53,  1195,
       0,    56,   741,     0,   742,   743,     0,   744,     0,     0,
       0,     7,     8,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    77,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,     0,     0,
       0,    41,    42,    43,    44,     0,     0,    46,    91,    92,
      93,     0,     0,     0,   740,    53,  1196,     0,    56,   741,
       0,   742,   743,     0,   744,     0,     0,     0,     7,     8,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    77,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    91,    92,    93,    41,    42,
      43,    44,     0,     0,    46,     0,     0,     0,     0,     0,
       0,   740,    53,  1215,     0,    56,   741,     0,   742,   743,
       0,   744,     0,     0,     0,     7,     8,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    77,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,     0,     0,     0,    41,    42,    43,    44,     0,
       0,    46,    91,    92,    93,     0,     0,     0,   740,    53,
    1246,     0,    56,   741,     0,   742,   743,     0,   744,     0,
       0,     0,     7,     8,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    77,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    91,
      92,    93,    41,    42,    43,    44,     0,     0,    46,     0,
       0,     0,     0,     0,     0,   740,    53,  1279,     0,    56,
     741,     0,   742,   743,     0,   744,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    77,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,     0,     0,     0,    41,
      42,    43,    44,     0,     0,    46,    91,    92,    93,     0,
       0,     0,   740,    53,  1431,     0,    56,   741,     0,   742,
     743,     0,   744,    -4,     1,     0,     0,    -4,     0,     0,
       0,     0,     0,     0,     0,     0,    -4,    -4,     0,    77,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,    91,    92,    93,     0,    -4,    -4,     0,
       0,     0,  1602,     0,     0,     0,     0,     0,     0,     0,
       0,  2051,     0,     0,     0,     0,     0,    -4,    -4,    -4,
       0,     0,     0,    -4,    -4,     0,    -4,     0,     0,     0,
       0,    -4,    -4,    -4,    -4,     0,    -4,    -4,     0,    -4,
       0,     0,     0,     0,    -4,    -4,    -4,    -4,    -4,    -4,
      -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,
       0,     0,    -4,    -4,    -4,    -4,    -4,    -4,  2110,    -4,
       0,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,
      -4,    -4,    -4,    -4,    -4,    -4,     0,     0,    -4,    -4,
      -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,
      -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,
      -4,    -4,     0,     6,     0,    -4,    -4,    -4,     0,     0,
       0,    -4,     7,     8,     0,     0,    -4,    -4,    -4,    -4,
       0,     0,    -4,     0,    -4,     0,    -4,    -4,    -4,    -4,
      -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,    -4,
       0,     0,     0,     9,    10,     0,    -4,    -4,    -4,     0,
       0,     0,     0,     0,     0,     0,     0,    -4,     0,    -4,
       0,     0,     0,    11,    12,    13,     0,     0,     0,    14,
      15,     0,    16,     0,     0,     0,     0,    17,    18,    19,
      20,     0,    21,    22,     0,    23,     0,     0,     0,     0,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,     0,     0,    40,    41,
      42,    43,    44,    45,     0,    46,     0,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,     0,     0,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,    82,    83,    84,    85,     0,     0,
       0,    86,    87,    88,     0,     0,     0,    89,     0,     0,
       0,     0,    90,    91,    92,    93,   160,   372,    94,     0,
      95,     0,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,     0,     0,     0,     0,
       0,     0,   110,   111,   112,     0,     7,     8,     0,     0,
       0,     0,     0,   113,     0,   114,     0,     0,   373,   374,
     375,   376,   377,   378,   379,   380,   381,     0,     0,     0,
       0,   382,   383,   384,   385,     0,     0,     0,     0,   386,
     387,   388,     0,     0,   389,     0,   390,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   391,
       0,     0,   392,     0,     0,     0,   292,     0,     0,   293,
       0,     0,   294,     0,   295,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    44,     0,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
       0,     0,   422,    41,    42,    43,    44,     0,     0,    46,
       0,     0,     0,     0,     0,     0,   740,    53,     0,     0,
      56,   741,     0,   742,   743,     0,   744,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    77,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,    91,    92,    93,
     393,   394,   395,     0,     0,     0,  1603,     0,     0,     0,
       0,   396,     0,     0,     0,   397,     0,   398,   114,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1617,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1748,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1866,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1867,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1895,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1898,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1965,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2019,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2020,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2034,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2036,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2069,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2097,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2102,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2103,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2104,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2144,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,     0,     0,     0,  1044,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,  1100,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,     0,     0,  1143,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,  1144,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,     0,     0,  1187,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,  1218,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
       0,     0,  1237,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
       0,  1293,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,     0,     0,
    1313,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,     0,  1435,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,     0,     0,     0,  1508,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,  1798,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,     0,     0,  1809,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,  1846,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,     0,     0,  1915,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,  1930,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
       0,     0,  1942,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
       0,  1999,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,     0,     0,
    2008,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,     0,  2009,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,     0,     0,     0,  2032,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,     0,  2083,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,     0,     0,  2123,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,     0,     0,  2142,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,     0,     0,  2162,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,     0,     0,  2163,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
       0,     0,  2164,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,     0,
     704,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,     0,   912,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,     0,  1908,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,   726,     0,   727,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
     972,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,  1020,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,  1167,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,  1232,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,  1233,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,     0,     0,  1238,   562,   563,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   703,
     576,   577,   578,   579,     0,     0,     0,     0,   580,     0,
       0,     0,  1239,   562,   563,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   703,   576,   577,   578,
     579,     0,     0,     0,     0,   580,     0,     0,     0,  1327,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   703,   576,   577,   578,   579,     0,     0,
       0,     0,   580,     0,     0,     0,  1341,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     703,   576,   577,   578,   579,     0,     0,     0,     0,   580,
       0,     0,     0,  1537,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,     0,     0,
    1623,   562,   563,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   703,   576,   577,   578,   579,     0,
       0,     0,     0,   580,     0,     0,     0,  1672,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   703,   576,   577,   578,   579,     0,     0,     0,     0,
     580,     0,     0,     0,  1858,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   703,   576,
     577,   578,   579,     0,     0,     0,     0,   580,     0,     0,
       0,  1900,   562,   563,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   703,   576,   577,   578,   579,
       0,     0,     0,     0,   580,     0,     0,     0,  1916,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   778,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   779,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   780,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   782,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   783,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   784,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   786,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   787,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   788,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   789,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   790,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   791,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   792,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   794,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   795,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   796,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   861,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   896,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   942,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   960,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   962,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   963,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   964,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,   970,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,   971,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1008,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1019,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1080,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1084,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1096,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1166,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1176,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1177,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1178,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1188,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1217,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1219,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1220,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1221,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1222,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1223,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1224,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1225,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1231,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1326,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1340,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1538,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1568,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1569,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1570,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1571,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1609,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1622,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1735,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1737,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1740,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1747,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1799,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1808,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1831,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1899,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  1963,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  1964,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580,     0,  2096,   562,
     563,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   703,   576,   577,   578,   579,     0,     0,     0,
       0,   580,     0,  2139,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   703,   576,   577,
     578,   579,     0,     0,     0,     0,   580
};

static const yytype_int16 yycheck[] =
{
       3,  1212,   307,     3,  1216,     4,    60,   617,   618,  1397,
       6,   721,  1659,  1703,     4,  1705,     4,   191,   319,   729,
      23,   708,     6,     4,     4,   199,     5,     4,     6,    23,
       5,   140,   153,     4,     4,   134,     5,   338,     6,     9,
       5,     7,     4,     4,     0,  1880,     4,     6,  1690,    52,
     737,    94,     6,     4,   100,     4,   232,   233,   247,     9,
       4,     6,   134,   105,   253,   107,     4,   754,   245,   100,
     113,     4,     4,   249,     6,   762,   253,   123,    81,   151,
       7,     7,   245,   182,     9,   154,    89,   159,   160,   161,
     253,    94,   123,   165,   166,   245,   245,   247,   140,     9,
       9,   250,   105,   253,   405,    86,   407,     9,   232,   233,
     179,  1753,     9,    14,   415,   289,   290,   291,   232,   233,
     244,     6,   244,   165,   248,  1772,   195,   196,     6,   134,
       6,   248,   244,    97,   246,   252,   100,   140,   102,   253,
       6,   141,   145,   146,    66,    67,   151,    69,     7,   113,
     153,   244,   103,   104,   159,   160,   161,   331,   332,   333,
     165,   166,    94,   100,   245,     7,   247,   170,   232,   233,
     232,   233,   253,     7,   348,   250,   863,   249,   253,   169,
     170,   171,   172,   186,   232,   233,   248,   245,   191,   253,
       6,   194,  2027,   232,   233,   253,   199,   200,   201,   167,
     244,   246,   232,   233,   248,   253,   232,   233,     7,   190,
     232,   233,     7,   167,   253,     6,   232,   233,   179,   180,
     244,   179,   180,   253,   248,     6,  1614,   253,   216,     6,
     245,   253,   194,   240,   241,   216,   216,   253,   253,   216,
    1451,   232,   233,   250,   249,   216,   216,   246,   251,   244,
     237,   238,   248,   252,   232,   233,   243,   248,  1905,   249,
     244,     6,   265,   251,   248,   268,   269,   248,   248,   248,
     240,   241,   216,   248,   244,   269,   246,   248,   248,   248,
     250,   284,   285,   248,   250,   288,   289,   290,   291,   248,
     217,   218,   219,   220,   244,   244,   246,     3,   301,     6,
     250,   245,   305,   247,   307,   606,   244,   608,   609,  1019,
    1020,   244,  1959,   240,   241,   240,   241,     6,   244,   244,
     246,   246,   247,   250,   625,   250,   329,     6,   331,   332,
     333,   240,   241,     4,   244,   244,   246,   246,   240,   241,
     250,   250,   244,   228,   246,   348,   231,   244,   250,   246,
     100,   354,   102,    59,   232,   233,   232,   233,   217,   218,
     219,   220,   228,   249,   665,   231,   540,   253,   247,   232,
     233,   232,   233,  2063,   253,   217,   218,   219,   220,   421,
     244,   240,   241,   217,   218,   219,   220,   245,   310,   311,
     253,   244,   253,   246,   246,   253,   318,   319,   240,   241,
     253,   232,   233,   704,   244,   232,   233,     8,   245,   412,
     413,   532,     4,   244,   417,   418,   253,   248,   217,   218,
     219,   220,   217,   218,   219,   220,   253,   232,   233,   134,
     433,   249,   232,   233,   140,   253,   244,   440,   246,   145,
     146,   240,   241,   232,   233,   248,   151,   153,   253,   245,
     232,   233,  2099,   253,   159,   160,   161,   253,   461,   244,
     165,   166,     6,  1150,   253,  1152,   134,   138,   139,   140,
     141,   253,   244,   144,    97,   232,   233,   100,  1165,   244,
     186,    97,   153,   151,   100,   156,   102,   490,   232,   233,
     113,   159,   160,   161,   115,   201,   253,   165,   166,   232,
     233,   244,    94,   246,   546,  2152,    98,   250,  1896,   253,
     253,   514,   104,   105,   244,   107,   108,   249,   521,   522,
     253,   253,   245,   526,   247,   528,   529,   530,   531,   532,
     253,   250,   535,   252,   240,   127,   539,   540,   541,   542,
     543,   544,   545,   244,   547,   249,   549,   253,   254,   253,
       6,   545,   673,   547,   244,    52,   232,   233,   232,   233,
     244,    94,   232,   233,    97,   232,   233,   100,    97,   102,
      94,   100,   575,    97,   100,   134,   100,   253,   581,   253,
     113,   249,    94,   253,   244,    97,   253,  1274,   100,   113,
     102,   244,   151,  1280,   232,   233,   240,   241,   601,   602,
     159,   160,   161,   245,   248,   247,   165,   166,   138,   139,
     248,   253,   721,   244,   144,   916,  1226,   617,   618,   244,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,  1325,     6,
       7,   245,   243,   247,   686,   648,   649,   650,   249,   253,
     240,   241,  1339,   234,   235,   236,   237,   238,  1345,   138,
     139,   244,   243,   249,   244,   144,   246,   253,   671,   672,
     673,   674,   675,   676,   677,   678,   679,   244,   681,   682,
     683,   684,   685,   248,   687,   688,   689,   690,   691,   743,
     249,   694,   102,   615,   616,   244,   690,   246,   620,   248,
     100,   249,   705,   200,   244,   253,   412,   232,   233,   234,
     235,   236,   237,   238,  1424,   718,   422,   245,   243,   247,
     240,   241,   244,   244,   246,   728,   248,     8,   248,   244,
     249,   734,   249,   249,   253,   777,   253,   253,   444,   249,
     249,   249,   249,   253,   253,   253,   253,   750,   751,   248,
     244,   249,   246,   756,   251,   253,   759,   760,   245,   249,
     247,   248,   765,   253,   767,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,   892,   249,   249,   249,   243,   253,   253,
     253,   288,   249,   244,   249,   246,   253,   800,   253,   802,
     248,  1511,   249,   249,   301,   249,   253,   253,   305,   253,
     249,   249,   249,  1523,   253,   253,   253,   249,    97,   249,
     526,   253,   528,   253,   530,   531,   532,   249,   244,   249,
     249,   253,   835,   253,   253,   249,   542,   543,   100,   253,
    2052,   140,   845,   244,   248,   246,   145,   146,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   906,   249,   249,   249,   243,   253,
     253,   253,   249,   249,  2085,   249,   253,   253,   881,   253,
     249,   249,     8,   248,   253,   253,   249,   186,   249,   249,
     253,   894,   253,   253,   249,     6,     7,   249,   253,  2111,
     903,   253,   201,   249,   249,   232,   233,   253,   253,   912,
    1019,  1020,   249,   248,    97,   249,   253,   623,  2129,   253,
     417,   418,     8,   249,   249,   249,   249,   253,   253,   253,
     253,   249,    97,   249,   249,   253,   433,   253,   253,   645,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,     6,   249,
     100,     4,   243,   253,   245,   671,   672,   673,   674,   675,
     676,   677,   678,   679,     4,   681,   682,   683,   684,  1021,
       4,   687,   688,   689,   245,     4,   247,   248,     4,   244,
       6,   697,   244,   244,   700,   244,   244,   919,     6,     6,
     246,  1306,   244,   246,   926,     6,   712,   248,   930,   248,
       6,  1014,   252,   244,   252,     9,   244,   244,   244,   244,
     244,   244,  1025,   244,   244,   244,   248,  1030,   182,   182,
     182,   248,   529,   123,  1644,  1645,   742,   244,   744,   244,
    1043,   244,   539,   244,   541,  1048,   182,   244,   248,   244,
     248,     4,   248,   759,   248,   248,   248,  1166,  1167,   244,
     244,   767,   244,   244,   244,   244,   244,   244,     6,     6,
     246,     6,     6,   248,   248,  1078,  1079,   248,  1081,  1082,
    1083,   248,  1085,  1086,  1087,  1088,  1089,  1090,  1091,  1092,
    1093,  1094,     6,   248,  1097,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,   412,   246,   246,   246,   243,   182,   244,
     244,   244,   244,   249,   244,   244,     6,     6,   246,     6,
       8,   248,  1135,  1136,  1137,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,     6,  1157,  1077,     8,   243,     7,   245,
     248,   247,     6,   248,   248,     6,   248,    89,   249,   245,
    1173,  1174,  1175,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,   253,   253,     7,     6,   243,     6,   248,   245,   248,
    1501,   243,    64,     8,  1505,   250,     4,     7,     7,   915,
     244,   917,   918,  1518,     6,     6,   245,  1326,  1327,     7,
     248,   718,     6,     6,  1227,   931,  1226,   526,     7,   528,
       6,   530,   531,   179,   249,   248,   247,   249,  1241,  1242,
    1243,     6,   245,   542,   543,   249,   248,     7,     6,  1252,
     250,   246,   244,     6,  1257,     6,     4,   248,  1261,     6,
       6,   967,     6,   245,     7,   246,     7,  1270,     7,  1272,
       7,     7,  1275,     7,     7,     7,     7,     7,     7,     6,
     248,     7,     7,     7,     7,     7,     6,  1290,   245,   247,
     253,   253,   253,   249,   253,   245,     7,  1003,   249,   248,
       7,  1911,   250,  1306,   248,  1347,   248,  1349,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,     4,  1327,   250,     6,   243,   249,  1332,
     249,  1334,   134,   248,     7,     6,   245,     7,     7,     7,
     250,   245,  1048,     9,   245,   253,  1349,   253,   253,   247,
     250,   182,     7,   252,   249,   248,     6,   154,     6,  1362,
       6,     4,    46,  1366,    46,   250,   248,   248,   244,   244,
     244,   244,     4,   672,     7,   674,   675,   676,   677,   678,
     679,   250,   681,   682,   683,   684,     6,   250,   687,   688,
     689,     7,     7,   253,   245,   182,     7,     7,     7,     6,
     245,     7,  1511,  1524,  1525,   253,   903,     7,     7,     4,
     112,     4,  1717,   248,     6,  1418,  1419,  1420,   244,     7,
       6,   248,     7,  1426,     7,     7,     7,     7,     7,     7,
    1136,  1137,     6,     6,     6,   100,     7,     6,  1441,     6,
       4,     4,   251,   249,   245,   253,   253,     6,     6,  1452,
       6,     6,   248,   248,   248,  1497,     7,  1460,     6,   244,
     759,  1464,   248,   244,   246,     6,     6,   127,   767,   249,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,  1194,   250,
       6,   247,   243,     6,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,  1621,   253,  1556,  1518,   243,  1014,     6,     6,
       6,  1524,  1525,     6,     6,     6,     6,  1530,  1450,     6,
       6,     6,     6,  1030,     6,     6,     5,   245,     6,   245,
       4,     6,     4,     6,   248,  1548,     7,   248,   248,   248,
     248,  1257,   246,  1556,     6,  1261,     6,   248,  1561,   248,
     248,   248,   248,   248,  1270,     6,  1272,     6,     6,  1275,
     249,   248,     6,   178,     6,     6,   248,   253,  1284,   253,
       6,  1078,  1079,   245,  1081,  1082,  1083,   250,  1085,  1086,
    1087,  1088,  1089,  1090,  1091,  1092,  1093,  1094,  1640,  1720,
     253,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,  1920,
       7,  1327,   253,   243,  1733,   244,   248,     4,     6,     5,
    1633,     6,   249,  1636,  1637,     6,     6,     4,  1135,     7,
     244,  1683,   244,  1948,  1644,  1645,     6,     6,  1690,  1691,
       6,     6,  1574,     6,     6,     6,  1659,     6,    98,     6,
    1157,   248,  1665,   253,   248,     6,     6,  1709,   245,     6,
       6,     6,  1675,   248,     6,     6,  1173,   253,     4,   245,
     248,     6,  1685,   253,     6,  1688,     6,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
     253,   253,     6,  1706,   243,     6,     6,     5,     7,   245,
    1713,   248,  2013,   248,  1717,   248,  1719,  1720,   248,     6,
       6,   248,     6,   248,     6,   249,   249,     7,   248,     6,
    1227,   177,     6,   245,   249,   249,     6,     6,  1660,  1661,
    1662,  1663,  1664,  1665,  1241,  1242,  1243,   249,     7,  1048,
    1753,   250,  2053,     6,     6,     6,   248,     6,     6,   249,
       6,   245,     6,   248,     6,   180,     6,     6,     6,  1772,
     248,     6,   134,     6,     6,     6,   249,   244,     6,   249,
       6,  1784,  1785,   249,     6,     6,   245,     6,   248,     6,
     248,   248,   248,  1290,   248,   248,   245,     6,   249,     6,
       6,   249,     6,  1845,  1807,   248,     6,     6,     6,     6,
       6,   248,   248,     6,     6,     6,  1819,  2065,  1524,  1525,
    1338,   458,  1392,  1843,  1566,  1881,  1361,  1629,     3,  1832,
       3,     3,     3,   596,  1907,  1494,     3,  1136,  1137,  1713,
    1690,   759,  1548,  1525,  1886,    -1,    -1,    -1,    -1,    -1,
      13,    14,  1349,     6,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1362,  1869,    -1,    -1,  1366,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1882,
      -1,    -1,    -1,    -1,  1887,  1888,  1889,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1905,    -1,    -1,  1908,    -1,    -1,    -1,    -1,
      -1,  1911,    -1,    -1,    -1,    -1,  1919,    -1,    -1,    -1,
    1923,  1418,  1419,  1420,    -1,    -1,    -1,    -1,    -1,  1426,
      -1,    94,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,  1441,  1948,    -1,   110,   111,   112,
     113,    -1,    -1,   116,    -1,  1452,  1959,    -1,  1257,    -1,
     123,   124,  1261,  1460,   127,   128,    -1,   130,   131,    -1,
     133,  1270,    -1,  1272,    -1,    -1,  1275,    -1,    -1,  1685,
      -1,    -1,  1688,  1986,    -1,    -1,  2028,   150,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1703,    -1,  1705,
      -1,    -1,  1924,  1925,  1926,  1927,  1928,    -1,    -1,    -1,
      -1,   174,   175,   176,  1720,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1327,    -1,
      -1,    -1,    -1,  1530,  2037,     6,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1753,    -1,    -1,
      -1,    -1,  2055,  2056,  2057,  2058,  2059,    -1,    -1,  1556,
      -1,    -1,    -1,  1985,    -1,    -1,    -1,  1989,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,   249,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,  2099,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2116,    -1,    -1,    -1,  2039,    -1,    -1,
      -1,  2043,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2138,  1633,  2059,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2152,
      -1,    -1,    -1,    -1,    -1,    -1,  2159,  2160,    -1,    -1,
      -1,    -1,  1659,    -1,    -1,   113,    -1,    -1,    -1,    -1,
    2173,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1675,    -1,
      -1,    -1,    -1,    -1,  2106,  2107,    -1,    -1,    -1,    -1,
      -1,    -1,   140,    -1,    -1,    -1,   144,    -1,    -1,    -1,
     148,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1706,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,     6,    -1,
      -1,    -1,  1719,   171,   172,   173,    -1,   175,   176,   177,
      -1,   179,   180,   181,   182,   183,   184,   185,    -1,   187,
     188,   189,   190,    -1,    -1,    -1,   194,    -1,    -1,  1548,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,  1772,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1998,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1807,    -1,    -1,    -1,    -1,    -1,   264,   265,   266,   267,
     268,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   283,   284,   285,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   303,    -1,  2063,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   313,   314,    -1,    -1,    -1,
      -1,   319,  1869,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2086,   329,   330,    -1,    -1,  1882,  1685,    -1,    -1,  1688,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   353,   354,    -1,  1905,    -1,
     358,   359,   360,    -1,   362,    -1,    -1,    -1,   366,   367,
     368,    -1,  1919,   371,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2138,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,  2159,  2160,   243,    -1,    -1,    -1,    -1,
      -1,   409,  1959,    -1,    -1,   413,   414,  2173,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   430,   431,    -1,    -1,    -1,    -1,    -1,  1986,
      -1,    -1,    -1,     6,    -1,    -1,    -1,   445,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   463,   464,   465,   466,    -1,
      -1,    -1,    -1,   471,    -1,    -1,    -1,   475,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   489,   490,    -1,    -1,    -1,   494,    -1,    -1,    -1,
     498,   499,   500,   501,   502,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,   517,
     518,   519,   520,   521,   522,    -1,   524,   525,    -1,   527,
      -1,    -1,    -1,    -1,    -1,    -1,   534,    -1,    -1,   537,
     538,    -1,    -1,    -1,    -1,    -1,   544,    -1,    -1,    -1,
      -1,   549,  2099,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   561,   562,   563,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
     578,   579,   580,    -1,   582,   583,    -1,    -1,    -1,    -1,
      -1,    -1,   590,   591,   592,    -1,    -1,    -1,    -1,    -1,
     598,   599,    -1,   601,   602,  2152,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   620,   621,   622,    -1,    -1,    -1,   626,   627,
     628,   629,   630,   631,   632,   633,   634,    -1,    -1,    -1,
      -1,    -1,   640,    -1,   642,    -1,   644,    -1,    -1,    -1,
     648,   649,   650,   651,    -1,   653,   654,   655,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,   680,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   694,    -1,    -1,    -1,
      -1,    -1,    -1,   701,    -1,   703,    -1,    -1,   706,   707,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   715,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     728,    -1,   730,    -1,    -1,    -1,    -1,   735,   736,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   750,    -1,   752,   753,    -1,    -1,   756,   757,
      -1,    -1,    -1,   761,    -1,    -1,    -1,    -1,    -1,    -1,
     768,     4,     5,    -1,   772,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2138,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     798,    -1,   800,    -1,   802,    -1,    -1,    -1,    -1,    -1,
    2159,  2160,    -1,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    -1,    -1,  2173,    -1,    59,    60,    61,    62,
      -1,     6,    -1,    -1,    67,    68,    69,   835,    -1,    72,
      -1,    74,    -1,    -1,    -1,    -1,    -1,   845,    -1,    -1,
      -1,    -1,    -1,    -1,    87,    -1,    -1,    90,    -1,    -1,
      -1,    -1,   194,    -1,    -1,    -1,    -1,    -1,    -1,   194,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   881,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   894,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   902,    -1,    -1,    -1,   906,    -1,
      -1,    -1,    -1,    -1,   912,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   920,   921,    -1,    -1,    -1,    -1,   926,    -1,
      -1,    -1,    -1,   265,    -1,    -1,   268,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   943,    -1,    -1,   946,    -1,
      -1,    -1,   284,   285,    -1,    -1,    -1,    -1,    -1,   284,
     285,    -1,    -1,    -1,    -1,   198,   199,   200,   966,    -1,
     968,   969,    -1,    -1,    -1,    -1,   209,    -1,    -1,    -1,
     213,    -1,   215,   216,    -1,    -1,    -1,    -1,    -1,   987,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   329,    -1,    -1,
      -1,    -1,    -1,    -1,   329,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1010,    -1,    -1,    -1,    -1,    -1,    -1,  1017,
      -1,    -1,   354,    -1,  1022,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1037,
    1038,    -1,    -1,    -1,    -1,  1043,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,    -1,    -1,  1072,  1073,  1074,  1075,    -1,    -1,
      -1,   413,    -1,    -1,    -1,    -1,    -1,    -1,   413,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1097,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1113,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1123,  1124,  1125,    -1,    -1,
      -1,    -1,  1130,  1131,    -1,  1133,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1145,  1146,    -1,
      -1,    -1,    -1,    -1,    -1,  1153,  1154,     6,   490,    -1,
      -1,    -1,  1160,    -1,  1162,  1163,    -1,    -1,    -1,    -1,
      -1,  1169,  1170,    -1,    -1,    -1,  1174,  1175,    -1,    -1,
      -1,    -1,   514,    -1,    -1,    -1,    -1,    -1,    -1,   521,
     522,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,     4,     5,    -1,    -1,    -1,  1206,    -1,
      -1,    -1,   544,    -1,    -1,    -1,    -1,   549,    -1,    -1,
      -1,    -1,    -1,    -1,   549,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1240,   575,    -1,    46,    47,    48,    49,    50,
      51,    52,    53,    54,  1252,    -1,    -1,    -1,    59,    60,
      61,    62,  1260,    -1,    -1,    -1,    67,    68,    69,   601,
     602,    72,    -1,    74,    -1,    -1,   601,   602,    -1,    -1,
      -1,    -1,    -1,  1281,  1282,    -1,    87,    -1,  1286,    90,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1297,
      -1,    -1,  1300,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   648,   649,   650,    -1,
      -1,    -1,    -1,   648,   649,   650,    -1,    -1,    -1,    -1,
    1328,    -1,    -1,  1331,  1332,  1333,  1334,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1342,    -1,  1344,    -1,  1346,  1347,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   694,    -1,    -1,    -1,    -1,    -1,    -1,   694,
      -1,    -1,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,   728,   198,   199,   200,
      -1,    -1,    -1,   728,    -1,    -1,    -1,    -1,   209,    -1,
      -1,    -1,   213,    -1,   215,    -1,    -1,    -1,   750,    -1,
      -1,    -1,    -1,    -1,   756,   750,    -1,  1425,    -1,    -1,
      -1,   756,    -1,    -1,  1432,  1433,  1434,    -1,    -1,    -1,
      -1,    -1,    -1,  1441,    -1,  1443,    -1,    -1,    -1,    -1,
      -1,    -1,  1450,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1464,    -1,   800,    -1,
     802,  1469,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1477,
    1478,  1479,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1489,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1499,    -1,   835,    -1,    -1,    -1,    -1,    -1,    -1,
     835,    -1,    -1,   845,    -1,    -1,    -1,    -1,    -1,    -1,
     845,    -1,    -1,    -1,  1522,    -1,    -1,    -1,    -1,    -1,
    1528,  1529,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   881,
      -1,    -1,    -1,    -1,    -1,    -1,   881,    -1,    -1,    -1,
      -1,    -1,   894,  1561,    -1,    -1,    -1,    -1,    -1,   894,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     912,    -1,    -1,    -1,    -1,    -1,    -1,   912,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1601,  1602,  1603,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1612,    -1,    -1,    -1,    -1,  1617,
    1618,  1619,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1631,    -1,    -1,    -1,    -1,  1636,  1637,
    1638,    -1,  1640,    -1,    -1,    -1,    -1,    -1,  1646,  1647,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1665,    -1,    -1,
      -1,    -1,  1670,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1683,    -1,    -1,  1686,    -1,
      -1,    -1,  1690,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1699,  1700,    -1,  1702,    -1,   435,    -1,    -1,    -1,
      -1,  1043,  1710,     5,    -1,    -1,    -1,    -1,  1043,     6,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1726,    -1,
      -1,    -1,    -1,  1731,  1732,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1745,    -1,    -1,
    1748,  1749,    -1,     6,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    -1,    -1,  1097,    -1,    59,    60,    61,
      62,    -1,  1097,    -1,    -1,    67,    68,    69,    -1,    -1,
      72,    -1,    74,    -1,    -1,    -1,  1784,  1785,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    87,    -1,    -1,    90,    -1,
      -1,    -1,    -1,  1801,    -1,    -1,   535,    -1,    -1,    -1,
      -1,    -1,  1810,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1819,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1828,  1829,    -1,    -1,  1832,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1174,  1175,    -1,    -1,    -1,    -1,    -1,  1174,
    1175,    -1,    -1,    -1,    -1,    -1,    -1,   586,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1866,  1867,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1885,    -1,    -1,
    1888,  1889,  1890,    -1,    -1,    -1,    -1,  1895,    -1,    -1,
    1898,    -1,    -1,    -1,    -1,    -1,   198,   199,   200,    -1,
    1908,    -1,    -1,    -1,    -1,    -1,    -1,   209,    -1,    -1,
    1252,   213,    -1,   215,    -1,  1923,    -1,  1252,    -1,    -1,
      -1,    -1,    -1,  1931,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,  1956,    -1,
      -1,    -1,    -1,    -1,  1962,    -1,    -1,  1965,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,  1994,    -1,    -1,    24,
    1332,    26,  1334,    -1,    -1,    -1,    -1,  1332,    -1,  1334,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2019,  2020,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2033,  2034,    -1,  2036,  2037,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       6,    -1,  2050,    -1,    -1,    -1,    -1,  2055,  2056,  2057,
    2058,  2059,   134,    -1,    -1,    -1,    -1,  2065,    -1,    -1,
      -1,  2069,    -1,    -1,    -1,    -1,    -1,    -1,   807,   808,
     809,   810,   811,   812,   813,   814,   815,   816,   817,    -1,
     819,   820,   821,   822,    -1,   824,   825,   826,   827,  2097,
       6,    -1,    -1,    -1,  2102,  2103,  2104,    -1,    -1,   838,
      -1,   840,    -1,    -1,    -1,   140,    -1,    -1,  2116,   848,
      -1,    -1,    -1,    -1,    -1,    -1,   855,   856,    -1,    -1,
      -1,    -1,  1464,    -1,    -1,   864,    -1,    -1,    -1,  1464,
      -1,    -1,    -1,    -1,    -1,    -1,  2144,  2145,  2146,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,    -1,   237,   238,    -1,    -1,    -1,    -1,    -1,   244,
     245,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1561,
      -1,   256,   257,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   973,   974,   975,    -1,    -1,    -1,
     979,   980,    -1,    -1,   983,   984,   985,   986,    -1,   988,
      -1,    -1,    -1,    -1,   993,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,     6,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1636,  1637,    -1,    -1,    -1,    -1,
      -1,  1636,  1637,    -1,    -1,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,  1665,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,  1071,    -1,    -1,    -1,    -1,  1076,   373,   374,
     375,    -1,    -1,    -1,   379,   380,   381,   382,   383,   384,
     385,    -1,   387,    -1,    -1,    -1,   391,   392,    -1,    -1,
     395,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   408,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,     3,     4,    -1,    -1,
      -1,    -1,    -1,    10,    11,    12,    -1,    -1,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1171,    -1,    -1,    -1,    -1,    -1,    55,    56,
      57,    58,  1784,  1785,    -1,    -1,    63,    -1,    -1,  1784,
    1785,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    76,
      77,    -1,     7,     8,    -1,    -1,    -1,    -1,    85,    86,
      -1,    88,    -1,    -1,    -1,    -1,    -1,  1819,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1832,    -1,    -1,    -1,    -1,    -1,    -1,  1832,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,   557,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,     7,    -1,    -1,    -1,   243,  1888,  1889,    -1,    -1,
      -1,   249,    -1,  1888,  1889,   253,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1908,    -1,    -1,    -1,
    1309,    -1,    -1,  1908,     8,    -1,    -1,    -1,    -1,    -1,
     197,  1923,    -1,    -1,   201,    -1,    -1,    -1,  1923,   206,
     207,   208,    -1,   210,   211,   212,    -1,    -1,    -1,   216,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1348,
      -1,    -1,     8,    -1,    -1,   232,   233,    -1,    -1,    -1,
      -1,    -1,   239,    -1,    -1,    -1,    -1,   244,    -1,    -1,
      -1,   248,    -1,    -1,   251,    -1,    -1,    -1,    -1,    -1,
      -1,  1380,  1381,  1382,  1383,  1384,    -1,    -1,    -1,    -1,
    1389,  1390,    -1,    -1,  1393,    -1,  1395,    -1,    -1,    -1,
    1399,    -1,    -1,  1402,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1414,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,  2037,    -1,    -1,   243,    -1,
      -1,    -1,  2037,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2055,  2056,  2057,  2058,  2059,    -1,    -1,
    2055,  2056,  2057,  2058,  2059,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,  1485,   243,  1487,    -1,
      -1,    -1,  1491,   249,  1493,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,  2116,    -1,    -1,   243,    -1,    -1,
      -1,  2116,    -1,    -1,    -1,  1524,    -1,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   859,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   893,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       3,     4,     5,    -1,    -1,    -1,    -1,    10,    11,    12,
      -1,  1620,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,     8,    -1,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    64,    -1,    -1,    67,    68,    69,    -1,    -1,    72,
      -1,    74,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    85,    86,    87,    88,    -1,    90,    -1,    -1,
      -1,    94,    -1,    -1,    97,    -1,    -1,   100,    -1,   102,
      -1,    -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,    -1,
     113,    -1,    -1,  1722,    -1,    -1,  1725,    -1,    -1,    -1,
     123,    -1,   125,   126,    -1,   128,    -1,   130,    -1,    -1,
     133,    -1,    -1,    -1,    -1,   138,   139,   140,   141,   142,
      -1,   144,   145,   146,   147,   148,   149,    -1,    -1,    -1,
     153,    -1,    -1,   156,    -1,    -1,    -1,    -1,    -1,    -1,
    1065,    -1,    -1,    -1,  1773,  1774,    -1,  1776,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,   245,    -1,   197,   198,   199,   200,   201,   202,
     253,    -1,    -1,   206,   207,   208,   209,   210,   211,   212,
     213,    -1,   215,   216,    -1,    -1,    -1,    -1,     8,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1836,  1837,   232,
     233,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,
      -1,   244,    -1,    -1,    -1,   248,    -1,    -1,   251,    -1,
      -1,    -1,    -1,  1158,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,  1883,    -1,   243,     3,     4,     5,
      -1,    -1,    -1,    -1,    10,    11,    12,    -1,    -1,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    -1,    -1,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,    60,    61,    62,    63,    64,    -1,
      -1,    67,    68,    69,    -1,    -1,    72,    -1,    74,    75,
      76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,
      86,    87,    88,    -1,    90,    -1,    -1,    -1,    94,    -1,
      -1,    97,    -1,    -1,   100,    -1,   102,    -1,    -1,    -1,
      -1,    -1,    -1,   109,    -1,    -1,    -1,   113,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   123,    -1,   125,
     126,    -1,   128,    -1,   130,  1310,  1311,   133,    -1,    -1,
      -1,    -1,   138,   139,   140,   141,   142,    -1,   144,   145,
     146,   147,   148,   149,    -1,    -1,    -1,   153,    -1,    -1,
     156,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,     8,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   197,   198,   199,   200,   201,   202,    -1,    -1,    -1,
     206,   207,   208,   209,   210,   211,   212,   213,    -1,   215,
     216,    -1,  1397,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1407,    -1,    -1,    -1,   232,   233,    -1,    -1,
      -1,    -1,    -1,   239,    -1,    -1,  2125,    -1,   244,     3,
       4,     5,   248,     7,    -1,   251,    10,    11,    12,    -1,
      -1,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    56,    57,    58,    59,    60,    61,    62,    63,
      -1,    -1,    -1,    67,    68,    69,    -1,    -1,    72,    -1,
      74,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    85,    86,    87,    88,    -1,    90,    -1,    -1,    -1,
      94,    -1,    -1,    97,    -1,    -1,   100,    -1,   102,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   113,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   130,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1560,    -1,    -1,    -1,    -1,
      -1,    -1,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,     8,    -1,   190,    -1,   243,    -1,
     245,    -1,   247,   197,   198,   199,   200,   201,   253,    -1,
      -1,    -1,   206,   207,   208,   209,   210,   211,   212,   213,
      -1,   215,   216,   217,   218,   219,   220,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,   233,
      -1,    -1,    -1,    -1,    -1,   239,   240,   241,    -1,    -1,
     244,    -1,   246,     3,     4,     5,   250,   251,    -1,    -1,
      10,    11,    12,    -1,    -1,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    55,    56,    57,    58,    -1,
      -1,    -1,    -1,    63,    64,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,    -1,
      -1,    -1,    -1,    -1,    94,    -1,    -1,    97,    -1,    -1,
     100,    -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,   109,
      -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   123,    -1,   125,   126,    -1,   128,    -1,
     130,    -1,    -1,   133,    -1,    -1,    -1,    -1,   138,   139,
     140,   141,   142,    -1,   144,   145,   146,   147,   148,   149,
      -1,    -1,    -1,   153,    -1,    -1,   156,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,
      -1,   201,   202,    -1,    -1,    -1,   206,   207,   208,    -1,
     210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   232,   233,    -1,    -1,    -1,    -1,    -1,   239,
      -1,    -1,    -1,    -1,   244,     3,     4,    -1,   248,    -1,
       8,   251,    10,    11,    12,    -1,    -1,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    56,    57,
      58,    -1,    -1,    -1,    -1,    63,    64,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    76,    77,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,    -1,
      88,    -1,    -1,    -1,    -1,    -1,    94,    -1,    -1,    97,
      -1,    -1,   100,    -1,   102,    -1,    -1,    -1,    -1,    -1,
      -1,   109,    -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   123,    -1,   125,   126,    -1,
     128,    -1,   130,    -1,    -1,   133,    -1,    -1,    -1,    -1,
     138,   139,   140,   141,   142,    -1,   144,   145,   146,   147,
     148,   149,    -1,    -1,    -1,   153,    -1,    -1,   156,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,   249,    -1,    -1,
      -1,   253,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,
      -1,    -1,    -1,   201,   202,    -1,    -1,    -1,   206,   207,
     208,    -1,   210,   211,   212,    -1,    -1,    -1,   216,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   232,   233,    -1,    -1,    -1,    -1,
      -1,   239,    -1,    -1,    -1,    -1,   244,     3,     4,     5,
      -1,   249,    -1,   251,    10,    11,    12,    -1,    -1,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,
      56,    57,    58,    -1,    -1,    -1,    -1,    63,    64,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,
      76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,
      86,    -1,    88,    -1,    -1,    -1,    -1,    -1,    94,    -1,
      -1,    97,    -1,    -1,   100,    -1,   102,    -1,    -1,    -1,
      -1,    -1,    -1,   109,    -1,    -1,    -1,   113,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   123,    -1,   125,
     126,    -1,   128,    -1,   130,    -1,    -1,   133,    -1,    -1,
      -1,    -1,   138,   139,   140,   141,   142,    -1,   144,   145,
     146,   147,   148,   149,    -1,    -1,    -1,   153,    -1,    -1,
     156,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,   249,
      -1,    -1,    -1,   253,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   197,    -1,    -1,    -1,   201,   202,    -1,    -1,    -1,
     206,   207,   208,    -1,   210,   211,   212,    -1,    -1,    -1,
     216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   232,   233,    -1,    -1,
      -1,    -1,    -1,   239,    -1,    -1,    -1,    -1,   244,     3,
       4,    -1,   248,    -1,    -1,   251,    10,    11,    12,    -1,
      -1,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    55,    56,    57,    58,    -1,    -1,    -1,    -1,    63,
      64,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    85,    86,    -1,    88,    -1,    -1,    -1,    -1,    -1,
      94,    -1,    -1,    97,    -1,    -1,   100,    -1,   102,    -1,
      -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,    -1,   113,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   123,
      -1,   125,   126,    -1,   128,    -1,   130,    -1,    -1,   133,
      -1,    -1,    -1,    -1,   138,   139,   140,   141,   142,    -1,
     144,   145,   146,   147,   148,   149,    -1,    -1,    -1,   153,
      -1,    -1,   156,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   253,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   197,    -1,    -1,    -1,   201,   202,    -1,
      -1,    -1,   206,   207,   208,    -1,   210,   211,   212,    -1,
      -1,    -1,   216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,   233,
      -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,    -1,
     244,     3,     4,    -1,   248,   249,    -1,   251,    10,    11,
      12,    -1,    -1,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    55,    56,    57,    58,    -1,    -1,    -1,
      -1,    63,    64,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    75,    76,    77,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    85,    86,    -1,    88,    -1,    -1,    -1,
      -1,    -1,    94,    -1,    -1,    97,    -1,    -1,   100,    -1,
     102,    -1,    -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,
      -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   123,    -1,   125,   126,    -1,   128,    -1,   130,    -1,
      -1,   133,    -1,    -1,    -1,    -1,   138,   139,   140,   141,
     142,    -1,   144,   145,   146,   147,   148,   149,    -1,    -1,
      -1,   153,    -1,    -1,   156,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,   201,
     202,    -1,    -1,    -1,   206,   207,   208,    -1,   210,   211,
     212,    -1,    -1,    -1,   216,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     232,   233,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,
      -1,    -1,   244,     3,     4,    -1,   248,    -1,    -1,   251,
      10,    11,    12,    -1,    -1,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    55,    56,    57,    58,    -1,
      -1,    -1,    -1,    63,    64,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,    -1,
      -1,    -1,    -1,    -1,    94,    -1,    -1,    97,    -1,    -1,
     100,    -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,   109,
      -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   123,    -1,   125,   126,    -1,   128,    -1,
     130,    -1,    -1,   133,    -1,    -1,    -1,    -1,   138,   139,
     140,   141,   142,    -1,   144,   145,   146,   147,   148,   149,
      -1,    -1,    -1,   153,    -1,    -1,   156,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,
      -1,   201,   202,    -1,    -1,    -1,   206,   207,   208,    -1,
     210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   232,   233,    -1,    -1,    -1,    -1,    -1,   239,
      -1,    -1,    -1,    -1,   244,     3,     4,    -1,   248,    -1,
      -1,   251,    10,    11,    12,    -1,    -1,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    56,    57,
      58,    -1,    -1,    -1,    -1,    63,    64,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    76,    77,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,    -1,
      88,    -1,    -1,    -1,    -1,    -1,    94,    -1,    -1,    97,
      -1,    -1,   100,    -1,   102,    -1,    -1,    -1,    -1,    -1,
      -1,   109,    -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   123,    -1,   125,   126,    -1,
     128,    -1,   130,    -1,    -1,   133,    -1,    -1,    -1,    -1,
     138,   139,   140,   141,   142,    -1,   144,   145,   146,   147,
     148,   149,    -1,    -1,    -1,   153,    -1,    -1,   156,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   253,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,
      -1,    -1,    -1,   201,   202,    -1,    -1,    -1,   206,   207,
     208,    -1,   210,   211,   212,    -1,    -1,    -1,   216,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   232,   233,    -1,    -1,    -1,    -1,
      -1,   239,    -1,    -1,    -1,    -1,   244,     3,     4,    -1,
     248,    -1,    -1,   251,    10,    11,    12,    -1,    -1,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,
      56,    57,    58,    -1,    -1,    -1,    -1,    63,    64,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,
      76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,
      86,    -1,    88,    -1,    -1,    -1,    -1,    -1,    94,    -1,
      -1,    97,    -1,    -1,   100,    -1,   102,    -1,    -1,    -1,
      -1,    -1,    -1,   109,    -1,    -1,    -1,   113,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   123,    -1,   125,
     126,    -1,   128,    -1,   130,    -1,    -1,   133,    -1,    -1,
      -1,    -1,   138,   139,   140,   141,   142,    -1,   144,   145,
     146,   147,   148,   149,    -1,    -1,    -1,   153,    -1,    -1,
     156,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   253,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   197,    -1,    -1,    -1,   201,   202,    -1,    -1,    -1,
     206,   207,   208,    -1,   210,   211,   212,    -1,    -1,    -1,
     216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   232,   233,    -1,    -1,
      -1,    -1,    -1,   239,    -1,    -1,    -1,    -1,   244,     3,
       4,    -1,    -1,   249,    -1,   251,    10,    11,    12,    -1,
      -1,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    55,    56,    57,    58,    -1,    -1,    -1,    -1,    63,
      64,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    85,    86,    -1,    88,    -1,    -1,    -1,    -1,    -1,
      94,    -1,    -1,    97,    -1,    -1,   100,    -1,   102,    -1,
      -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,    -1,   113,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   123,
      -1,   125,   126,    -1,   128,    -1,   130,    -1,    -1,   133,
      -1,    -1,    -1,    -1,   138,   139,   140,   141,   142,    -1,
     144,   145,   146,   147,   148,   149,    -1,    -1,    -1,   153,
      -1,    -1,   156,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   253,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   197,    -1,    -1,    -1,   201,   202,    -1,
      -1,    -1,   206,   207,   208,    -1,   210,   211,   212,    -1,
      -1,    -1,   216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,   233,
      -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,    -1,
     244,     3,     4,    -1,   248,    -1,    -1,   251,    10,    11,
      12,    -1,    -1,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    55,    56,    57,    58,    -1,    -1,    -1,
      -1,    63,    64,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    75,    76,    77,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    85,    86,    -1,    88,    -1,    -1,    -1,
      -1,    -1,    94,    -1,    -1,    97,    -1,    -1,   100,    -1,
     102,    -1,    -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,
      -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   123,    -1,   125,   126,    -1,   128,    -1,   130,    -1,
      -1,   133,    -1,    -1,    -1,    -1,   138,   139,   140,   141,
     142,    -1,   144,   145,   146,   147,   148,   149,    -1,    -1,
      -1,   153,    -1,    -1,   156,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,   201,
     202,    -1,    -1,    -1,   206,   207,   208,    -1,   210,   211,
     212,    -1,    -1,    -1,   216,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     232,   233,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,
      -1,    -1,   244,     3,     4,    -1,   248,    -1,    -1,   251,
      10,    11,    12,    -1,    -1,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    55,    56,    57,    58,    -1,
      -1,    -1,    -1,    63,    64,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,    -1,
      -1,    -1,    -1,    -1,    94,    -1,    -1,    97,    -1,    -1,
     100,    -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,   109,
      -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   123,    -1,   125,   126,    -1,   128,    -1,
     130,    -1,    -1,   133,    -1,    -1,    -1,    -1,   138,   139,
     140,   141,   142,    -1,   144,   145,   146,   147,   148,   149,
      -1,    -1,    -1,   153,    -1,    -1,   156,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,
      -1,   201,   202,    -1,    -1,    -1,   206,   207,   208,    -1,
     210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   232,   233,    -1,    -1,    -1,    -1,    -1,   239,
      -1,    -1,    -1,    -1,   244,     3,     4,    -1,   248,    -1,
       8,   251,    10,    11,    12,    -1,    -1,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    56,    57,
      58,    -1,    -1,    -1,    -1,    63,    64,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    76,    77,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,    -1,
      88,    -1,    -1,    -1,    -1,    -1,    94,    -1,    -1,    97,
      -1,    -1,   100,    -1,   102,    -1,    -1,    -1,    -1,    -1,
      -1,   109,    -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   123,    -1,   125,   126,    -1,
     128,    -1,   130,    -1,    -1,   133,    -1,    -1,    -1,    -1,
     138,   139,   140,   141,   142,    -1,   144,   145,   146,   147,
     148,   149,    -1,    -1,    -1,   153,    -1,    -1,   156,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   253,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,
      -1,    -1,    -1,   201,   202,    -1,    -1,    -1,   206,   207,
     208,    -1,   210,   211,   212,    -1,    -1,    -1,   216,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   232,   233,    -1,    -1,    -1,    -1,
      -1,   239,     3,     4,    -1,    -1,   244,     8,    -1,    10,
      11,    12,    -1,   251,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    55,    56,    57,    58,    -1,    -1,
      -1,    -1,    63,    64,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    75,    76,    77,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    85,    86,    -1,    88,    -1,    -1,
      -1,    -1,    -1,    94,    -1,    -1,    97,    -1,    -1,   100,
      -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,   109,    -1,
      -1,    -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   123,    -1,   125,   126,    -1,   128,    -1,   130,
      -1,    -1,   133,    -1,    -1,    -1,    -1,   138,   139,   140,
     141,   142,    -1,   144,   145,   146,   147,   148,   149,    -1,
      -1,    -1,   153,    -1,    -1,   156,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,
     201,   202,    -1,    -1,    -1,   206,   207,   208,    -1,   210,
     211,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   232,   233,    -1,    -1,    -1,    -1,    -1,   239,     3,
       4,    -1,    -1,   244,    -1,    -1,    10,    11,    12,    -1,
     251,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    55,    56,    57,    58,    -1,    -1,    -1,    -1,    63,
      64,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    85,    86,    -1,    88,    -1,    -1,    -1,    -1,    -1,
      94,    -1,    -1,    97,    -1,    -1,   100,    -1,   102,    -1,
      -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,    -1,   113,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   123,
      -1,   125,   126,    -1,   128,    -1,   130,    -1,    -1,   133,
      -1,    -1,    -1,    -1,   138,   139,   140,   141,   142,    -1,
     144,   145,   146,   147,   148,   149,    -1,    -1,    -1,   153,
      -1,    -1,   156,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   253,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   197,    -1,    -1,    -1,   201,   202,    -1,
      -1,    -1,   206,   207,   208,    -1,   210,   211,   212,    -1,
      -1,    -1,   216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,   233,
      -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,    -1,
     244,     3,     4,    -1,   248,    -1,    -1,   251,    10,    11,
      12,    -1,    -1,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    55,    56,    57,    58,    -1,    -1,    -1,
      -1,    63,    64,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    75,    76,    77,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    85,    86,    -1,    88,    -1,    -1,    -1,
      -1,    -1,    94,    -1,    -1,    97,    -1,    -1,   100,    -1,
     102,    -1,    -1,    -1,    -1,    -1,    -1,   109,    -1,    -1,
      -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   123,    -1,   125,   126,    -1,   128,    -1,   130,    -1,
      -1,   133,    -1,    -1,    -1,    -1,   138,   139,   140,   141,
     142,    -1,   144,   145,   146,   147,   148,   149,    -1,    -1,
      -1,   153,    -1,    -1,   156,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,   201,
     202,    -1,    -1,    -1,   206,   207,   208,    -1,   210,   211,
     212,    -1,    -1,    -1,   216,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     232,   233,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,
      -1,    -1,   244,     3,     4,    -1,   248,    -1,    -1,   251,
      10,    11,    12,    -1,    -1,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    55,    56,    57,    58,    -1,
      -1,    -1,    -1,    63,    64,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,    -1,
      -1,    -1,    -1,    -1,    94,    -1,    -1,    97,    -1,    -1,
     100,    -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,   109,
      -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   123,    -1,   125,   126,    -1,   128,    -1,
     130,    -1,    -1,   133,    -1,    -1,    -1,    -1,   138,   139,
     140,   141,   142,    -1,   144,   145,   146,   147,   148,   149,
      -1,    -1,    -1,   153,    -1,    -1,   156,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,
      -1,   201,   202,    -1,    -1,    -1,   206,   207,   208,    -1,
     210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   232,   233,    -1,    -1,    -1,    -1,    -1,   239,
       3,     4,     5,    -1,   244,    -1,    -1,    10,    11,    12,
      -1,   251,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    -1,    -1,    -1,    67,    68,    69,    -1,    -1,    72,
      -1,    74,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    85,    86,    87,    88,    -1,    90,     3,     4,
       5,    -1,    -1,    -1,    -1,    10,    11,    12,    -1,    -1,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    -1,
      -1,    -1,    67,    68,    69,    -1,    -1,    72,    -1,    74,
      75,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      85,    86,    87,    88,    -1,    90,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   197,   198,   199,   200,   201,    -1,
      -1,    -1,    -1,   206,   207,   208,   209,   210,   211,   212,
     213,    -1,   215,   216,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,
     233,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,
      -1,   244,    -1,    -1,    -1,   248,    -1,    -1,   251,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   253,   197,   198,   199,   200,   201,    -1,    -1,    -1,
      -1,   206,   207,   208,   209,   210,   211,   212,   213,    -1,
     215,   216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,   233,    -1,
      -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,    -1,   244,
       3,     4,     5,   248,    -1,    -1,   251,    10,    11,    12,
      -1,    -1,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    -1,    -1,    -1,    67,    68,    69,    -1,    -1,    72,
      -1,    74,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    85,    86,    87,    88,    -1,    90,     3,     4,
       5,    -1,    -1,    -1,    -1,    10,    11,    12,    -1,    -1,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    -1,
      -1,    -1,    67,    68,    69,    -1,    -1,    72,    -1,    74,
      75,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      85,    86,    87,    88,    -1,    90,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   197,   198,   199,   200,   201,    -1,
      -1,    -1,    -1,   206,   207,   208,   209,   210,   211,   212,
     213,    -1,   215,   216,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,
     233,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,
      -1,   244,    -1,    -1,    -1,   248,    -1,    -1,   251,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   253,   197,   198,   199,   200,   201,    -1,    -1,    -1,
      -1,   206,   207,   208,   209,   210,   211,   212,   213,    -1,
     215,   216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,   233,    -1,
      -1,    -1,    -1,    -1,   239,     3,     4,    -1,    -1,   244,
      -1,    -1,    10,    11,    12,    -1,   251,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    56,    57,
      58,    -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    76,    77,
      -1,    -1,    -1,    -1,     3,     4,    -1,    85,    86,     8,
      88,    10,    11,    12,    -1,    -1,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    55,    56,    57,    58,
      -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,   197,
      -1,    -1,   253,   201,    -1,    -1,    -1,    -1,   206,   207,
     208,    -1,   210,   211,   212,    -1,    -1,    -1,   216,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   232,   233,    -1,    -1,    -1,    -1,
      -1,   239,    -1,    -1,    -1,    -1,   244,    -1,    -1,    -1,
     248,    -1,    -1,   251,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   197,    -1,
      -1,    -1,   201,    -1,    -1,    -1,   253,   206,   207,   208,
      -1,   210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   232,   233,    -1,    -1,    -1,    -1,    -1,
     239,     3,     4,    -1,    -1,   244,    -1,    -1,    10,    11,
      12,    -1,   251,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    55,    56,    57,    58,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    75,    76,    77,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    85,    86,    -1,    88,     3,     4,    -1,
       6,    -1,    -1,    -1,    10,    11,    12,    -1,    -1,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    13,    14,    -1,    55,
      56,    57,    58,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,
      76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,
      86,    -1,    88,    -1,    -1,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,   197,    -1,   243,    -1,   201,
      -1,    -1,    -1,    -1,   206,   207,   208,   253,   210,   211,
     212,    -1,    -1,    -1,   216,    -1,    -1,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     232,   233,    -1,   110,   111,   112,   113,   239,    -1,   116,
      -1,    -1,   244,   245,    -1,    -1,   123,   124,    -1,   251,
     127,   128,    -1,   130,   131,    -1,   133,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   150,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   197,    -1,    -1,    -1,   201,    -1,    -1,    -1,    -1,
     206,   207,   208,    -1,   210,   211,   212,   174,   175,   176,
     216,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   232,   233,    -1,    -1,
      -1,    -1,    -1,   239,     3,     4,    -1,     6,   244,    -1,
      -1,    10,    11,    12,    -1,   251,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   249,    -1,    -1,    -1,    55,    56,    57,    58,
      -1,    -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,
       3,     4,    -1,     6,    -1,    -1,    -1,    10,    11,    12,
      -1,    -1,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    13,
      14,    -1,    55,    56,    57,    58,    -1,    -1,    -1,    -1,
      63,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    75,    76,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    85,    86,    -1,    88,    -1,    -1,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,   197,    -1,
     243,    -1,   201,    -1,    -1,    -1,    -1,   206,   207,   208,
     253,   210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,
      94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   232,   233,    -1,   110,   111,   112,   113,
     239,    -1,   116,    -1,    -1,   244,    -1,    -1,    -1,   123,
     124,    -1,   251,   127,   128,    -1,   130,   131,    -1,   133,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   150,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,   201,    -1,
      -1,    -1,    -1,   206,   207,   208,    -1,   210,   211,   212,
     174,   175,   176,   216,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,
     233,    -1,    -1,    -1,    -1,    -1,   239,     3,     4,    -1,
      -1,   244,    -1,    -1,    10,    11,    12,    -1,   251,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   249,    -1,    -1,    -1,    55,
      56,    57,    58,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,
      76,    77,    -1,    -1,    -1,    -1,     3,     4,    -1,    85,
      86,    -1,    88,    10,    11,    12,    -1,    -1,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    56,
      57,    58,    13,    14,    -1,    -1,    63,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    75,    76,
      77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,
      -1,    88,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,
      -1,   197,    -1,    -1,   253,   201,    -1,    -1,    -1,    -1,
     206,   207,   208,    -1,   210,   211,   212,    -1,    -1,    -1,
     216,    -1,    -1,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   232,   233,    -1,   110,
     111,   112,   113,   239,    -1,   116,    -1,    -1,   244,   245,
      -1,    -1,   123,   124,    -1,   251,   127,   128,    -1,   130,
     131,    -1,   133,    -1,    -1,    -1,   137,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   150,
     197,    -1,    -1,   154,   201,    -1,    -1,    -1,    -1,   206,
     207,   208,    -1,   210,   211,   212,    -1,    -1,    -1,   216,
      -1,    -1,    -1,   174,   175,   176,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   232,   233,    -1,    -1,    -1,
      -1,    -1,   239,     3,     4,    -1,    -1,   244,   245,    -1,
      10,    11,    12,    -1,   251,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,     4,     5,    -1,
      -1,    -1,    -1,    -1,    -1,    55,    56,    57,    58,    -1,
      -1,    -1,    -1,    63,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    75,    76,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    88,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    -1,    -1,
      -1,    -1,    59,    60,    61,    62,    -1,    -1,    -1,    -1,
      67,    68,    69,    -1,    -1,    72,    -1,    74,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      87,    -1,    -1,    90,    -1,    -1,    -1,    94,    -1,    -1,
      97,    -1,    -1,   100,    -1,   102,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   113,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   130,    -1,    -1,    -1,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,   197,    -1,   243,
      -1,   201,    -1,    -1,    -1,    -1,   206,   207,   208,   253,
     210,   211,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    -1,     4,     5,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   232,   233,    -1,    -1,    -1,    -1,    -1,   239,
      -1,   198,   199,   200,   244,    -1,    -1,    -1,    -1,    -1,
      -1,   251,   209,    -1,    -1,    -1,   213,    -1,   215,   216,
      -1,    -1,    -1,    -1,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    -1,    -1,    -1,    -1,    59,    60,    61,
      62,    -1,    -1,   240,   241,    67,    68,    69,    -1,    -1,
      72,   248,    74,   250,     4,     5,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    87,    -1,    -1,    90,    -1,
      -1,    -1,    94,    -1,    -1,    97,    -1,    -1,   100,    -1,
     102,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   113,    -1,    -1,    -1,    -1,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    -1,    -1,    -1,   130,    59,
      60,    61,    62,    -1,    -1,    -1,    -1,    67,    68,    69,
      -1,    -1,    72,    -1,    74,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    87,    -1,    -1,
      90,    -1,    -1,    -1,    94,    -1,    -1,    97,    -1,    -1,
     100,    -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   113,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   198,   199,   200,    -1,
     130,    -1,    -1,    -1,    -1,    -1,    -1,   209,    -1,    13,
      14,   213,    -1,   215,   216,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,   245,    -1,   247,   248,   253,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    13,    14,   198,   199,
     200,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   209,
      -1,    -1,    -1,   213,    -1,   215,   216,    -1,    -1,    -1,
      94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,    -1,    -1,    -1,   110,   111,   112,   113,
      -1,    -1,   116,    -1,    -1,    -1,    -1,    -1,   248,   123,
     124,    -1,    -1,   127,   128,    -1,   130,   131,    -1,   133,
      -1,    -1,    -1,    13,    14,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   150,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
      -1,    -1,    -1,   110,   111,   112,   113,    -1,    -1,   116,
     174,   175,   176,    -1,    -1,    -1,   123,   124,    -1,    -1,
     127,   128,    -1,   130,   131,    -1,   133,    -1,    -1,    -1,
      13,    14,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   150,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   174,   175,   176,
     110,   111,   112,   113,    -1,    -1,   116,    -1,    -1,    -1,
      -1,    -1,    -1,   123,   124,   249,    -1,   127,   128,    -1,
     130,   131,    -1,   133,    -1,    -1,    -1,    13,    14,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     150,    94,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,    -1,    -1,    -1,   110,   111,   112,
     113,    -1,    -1,   116,   174,   175,   176,    -1,    -1,    -1,
     123,   124,   249,    -1,   127,   128,    -1,   130,   131,    -1,
     133,    -1,    -1,    -1,    13,    14,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   150,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    94,    95,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   174,   175,   176,   110,   111,   112,   113,    -1,    -1,
     116,    -1,    -1,    -1,    -1,    -1,    -1,   123,   124,   249,
      -1,   127,   128,    -1,   130,   131,    -1,   133,    -1,    -1,
      -1,    13,    14,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   150,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,    -1,    -1,
      -1,   110,   111,   112,   113,    -1,    -1,   116,   174,   175,
     176,    -1,    -1,    -1,   123,   124,   249,    -1,   127,   128,
      -1,   130,   131,    -1,   133,    -1,    -1,    -1,    13,    14,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   150,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    94,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   174,   175,   176,   110,   111,
     112,   113,    -1,    -1,   116,    -1,    -1,    -1,    -1,    -1,
      -1,   123,   124,   249,    -1,   127,   128,    -1,   130,   131,
      -1,   133,    -1,    -1,    -1,    13,    14,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   150,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,    -1,    -1,    -1,   110,   111,   112,   113,    -1,
      -1,   116,   174,   175,   176,    -1,    -1,    -1,   123,   124,
     249,    -1,   127,   128,    -1,   130,   131,    -1,   133,    -1,
      -1,    -1,    13,    14,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   150,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    94,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   174,
     175,   176,   110,   111,   112,   113,    -1,    -1,   116,    -1,
      -1,    -1,    -1,    -1,    -1,   123,   124,   249,    -1,   127,
     128,    -1,   130,   131,    -1,   133,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   150,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,    -1,    -1,    -1,   110,
     111,   112,   113,    -1,    -1,   116,   174,   175,   176,    -1,
      -1,    -1,   123,   124,   249,    -1,   127,   128,    -1,   130,
     131,    -1,   133,     0,     1,    -1,    -1,     4,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    13,    14,    -1,   150,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,   174,   175,   176,    -1,    44,    45,    -1,
      -1,    -1,   253,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   249,    -1,    -1,    -1,    -1,    -1,    64,    65,    66,
      -1,    -1,    -1,    70,    71,    -1,    73,    -1,    -1,    -1,
      -1,    78,    79,    80,    81,    -1,    83,    84,    -1,    86,
      -1,    -1,    -1,    -1,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
      -1,    -1,   109,   110,   111,   112,   113,   114,   249,   116,
      -1,   118,   119,   120,   121,   122,   123,   124,   125,   126,
     127,   128,   129,   130,   131,   132,    -1,    -1,   135,   136,
     137,   138,   139,   140,   141,   142,   143,   144,   145,   146,
     147,   148,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,    -1,     4,    -1,   162,   163,   164,    -1,    -1,
      -1,   168,    13,    14,    -1,    -1,   173,   174,   175,   176,
      -1,    -1,   179,    -1,   181,    -1,   183,   184,   185,   186,
     187,   188,   189,   190,   191,   192,   193,   194,   195,   196,
      -1,    -1,    -1,    44,    45,    -1,   203,   204,   205,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   214,    -1,   216,
      -1,    -1,    -1,    64,    65,    66,    -1,    -1,    -1,    70,
      71,    -1,    73,    -1,    -1,    -1,    -1,    78,    79,    80,
      81,    -1,    83,    84,    -1,    86,    -1,    -1,    -1,    -1,
      91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,    -1,    -1,   109,   110,
     111,   112,   113,   114,    -1,   116,    -1,   118,   119,   120,
     121,   122,   123,   124,   125,   126,   127,   128,   129,   130,
     131,   132,    -1,    -1,   135,   136,   137,   138,   139,   140,
     141,   142,   143,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,    -1,    -1,
      -1,   162,   163,   164,    -1,    -1,    -1,   168,    -1,    -1,
      -1,    -1,   173,   174,   175,   176,     4,     5,   179,    -1,
     181,    -1,   183,   184,   185,   186,   187,   188,   189,   190,
     191,   192,   193,   194,   195,   196,    -1,    -1,    -1,    -1,
      -1,    -1,   203,   204,   205,    -1,    13,    14,    -1,    -1,
      -1,    -1,    -1,   214,    -1,   216,    -1,    -1,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    -1,    -1,    -1,
      -1,    59,    60,    61,    62,    -1,    -1,    -1,    -1,    67,
      68,    69,    -1,    -1,    72,    -1,    74,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    87,
      -1,    -1,    90,    -1,    -1,    -1,    94,    -1,    -1,    97,
      -1,    -1,   100,    -1,   102,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   113,    -1,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
      -1,    -1,   130,   110,   111,   112,   113,    -1,    -1,   116,
      -1,    -1,    -1,    -1,    -1,    -1,   123,   124,    -1,    -1,
     127,   128,    -1,   130,   131,    -1,   133,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   150,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,   174,   175,   176,
     198,   199,   200,    -1,    -1,    -1,   253,    -1,    -1,    -1,
      -1,   209,    -1,    -1,    -1,   213,    -1,   215,   216,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   253,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   253,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     253,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   253,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   253,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   253,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   253,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   253,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   253,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     253,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   253,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   253,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   253,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,   249,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,   249,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,   249,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,   249,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,    -1,    -1,   249,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,   249,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
      -1,    -1,   249,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
      -1,   249,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,
     249,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,   249,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,   249,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,   249,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,   249,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,   249,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,    -1,    -1,   249,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,   249,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
      -1,    -1,   249,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
      -1,   249,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,
     249,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,   249,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,    -1,    -1,    -1,   249,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,    -1,   249,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,    -1,    -1,   249,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,   249,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,    -1,    -1,   249,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,    -1,    -1,   249,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
      -1,    -1,   249,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,
     248,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,    -1,   248,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,    -1,   248,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,   245,    -1,   247,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
     247,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,   247,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,   247,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,   247,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,   247,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,    -1,    -1,   247,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,
      -1,    -1,   247,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,   247,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,   235,   236,   237,   238,    -1,    -1,
      -1,    -1,   243,    -1,    -1,    -1,   247,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,   247,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,
     247,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,   235,   236,   237,   238,    -1,
      -1,    -1,    -1,   243,    -1,    -1,    -1,   247,   221,   222,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   234,   235,   236,   237,   238,    -1,    -1,    -1,    -1,
     243,    -1,    -1,    -1,   247,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,   235,
     236,   237,   238,    -1,    -1,    -1,    -1,   243,    -1,    -1,
      -1,   247,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
      -1,    -1,    -1,    -1,   243,    -1,    -1,    -1,   247,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243,    -1,   245,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,   231,
     232,   233,   234,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,   243,    -1,   245,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,    -1,    -1,    -1,    -1,   243
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     1,   255,   256,     6,     0,     4,    13,    14,    44,
      45,    64,    65,    66,    70,    71,    73,    78,    79,    80,
      81,    83,    84,    86,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     109,   110,   111,   112,   113,   114,   116,   118,   119,   120,
     121,   122,   123,   124,   125,   126,   127,   128,   129,   130,
     131,   132,   135,   136,   137,   138,   139,   140,   141,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   162,   163,   164,   168,
     173,   174,   175,   176,   179,   181,   183,   184,   185,   186,
     187,   188,   189,   190,   191,   192,   193,   194,   195,   196,
     203,   204,   205,   214,   216,   257,   259,   260,   280,   298,
     300,   304,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   322,   324,   325,   331,   332,   333,   334,   340,   365,
     366,   248,   252,    14,   100,   244,   244,     6,   248,     6,
       6,     6,     6,   244,     6,     6,     6,     6,   246,   246,
       4,   342,   366,   244,   246,   278,    94,    97,   100,   102,
     278,   244,   244,   244,     4,   244,   244,   244,     4,   244,
     244,   244,   244,   244,   244,   244,   244,   244,   244,   244,
     248,   115,   100,     6,   248,    94,    97,   100,   113,   303,
     102,   244,     3,    10,    11,    12,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    55,    56,    57,    58,    63,
      64,    75,    76,    77,    85,    88,    94,    97,   100,   102,
     113,   123,   128,   130,   133,   197,   201,   202,   206,   207,
     208,   210,   211,   212,   232,   233,   239,   244,   248,   251,
     300,   301,   304,   315,   322,   324,   335,   336,   340,   342,
     349,   351,   366,   244,   248,   248,   100,   100,   123,    97,
     100,   102,    94,    97,   100,   102,   300,    97,   100,   102,
     113,   301,    97,   100,   244,    97,   154,   179,   195,   196,
     248,   232,   233,   244,   248,   346,   347,   346,   248,   248,
     346,     4,    94,    98,   104,   105,   107,   108,   127,   248,
     244,   100,   102,   100,    97,     4,    86,   190,   248,   366,
       4,     6,    94,    97,   100,    97,   100,   113,   302,     4,
       4,     4,     5,   244,   248,   349,   350,     4,   244,   244,
     244,     4,   248,   353,   366,     4,   244,   244,   244,     6,
       6,   246,     5,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    59,    60,    61,    62,    67,    68,    69,    72,
      74,    87,    90,   198,   199,   200,   209,   213,   215,   357,
     366,   244,     4,   357,     5,   248,     5,   248,    32,   233,
     335,   366,   246,   248,   244,   248,     6,   244,   248,     6,
     252,     7,   130,   190,   217,   218,   219,   220,   240,   241,
     244,   246,   250,   276,   277,   278,   300,   335,   356,   357,
     366,     4,   304,   305,   306,   248,     6,   335,   356,   357,
     366,   356,   335,   356,   363,   364,   366,   282,   286,   244,
     345,     9,   357,   244,   244,   244,   244,   366,   335,   335,
     335,   244,   335,   335,   335,   244,   335,   335,   335,   335,
     335,   335,   335,   356,   335,   335,   335,   335,   350,   244,
     233,   335,   351,   352,   248,   350,   349,   356,   278,   278,
     278,   278,   278,   278,   278,   278,   278,   278,   278,   278,
     278,   278,   278,   278,   278,   278,   278,   278,   278,   278,
     278,   244,   246,   278,   278,   278,   278,   278,   278,   244,
     278,   278,   244,   300,   278,   278,     5,   248,   248,   123,
     300,   300,   244,   278,   278,   244,   244,   244,   335,   248,
     335,   351,   335,   335,   249,   352,   342,   366,   182,     5,
     248,     8,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,   235,   236,   237,   238,
     243,     9,   244,   246,   250,   277,   278,   335,   352,   352,
     244,   244,   244,   349,   350,   350,   350,   299,   248,   244,
     349,   248,   248,   335,     4,   349,   248,   353,   248,   248,
     346,   346,   346,   335,   335,   232,   233,   248,   248,   346,
     232,   233,   244,   306,   346,   248,   244,   248,   244,   244,
     244,   244,   244,   244,   244,   352,   335,   350,   350,   350,
     244,     4,   246,     6,   246,   306,     6,     6,   248,   248,
     248,   248,   350,   246,   246,   246,   335,     8,     6,     6,
     335,   335,   335,   250,   335,   248,   182,   335,   335,   335,
     335,   278,   278,   278,   244,   244,   244,   278,   278,   278,
     278,   278,   278,   278,   278,   278,   278,   244,   244,   278,
     244,   246,     6,     6,   248,     6,     8,   306,     6,     8,
     306,   278,   335,   234,   248,     9,   244,   246,   250,   356,
     352,   335,   306,   349,   349,   248,   357,   300,     7,   335,
     335,     4,   179,   180,   349,     6,   245,   247,   248,   279,
     248,     6,   248,     6,     9,   244,   246,   250,   366,   249,
     123,   128,   130,   131,   133,   298,   300,   335,     6,   245,
     253,     9,   244,   246,   250,   245,   253,   253,   245,   253,
       9,   244,   250,   247,   253,   281,   247,   281,    89,   344,
     341,   366,   253,   335,   335,   335,   335,   253,   245,   245,
     245,   335,   245,   245,   245,   335,   245,   245,   245,   245,
     245,   245,   245,   245,   245,   245,   245,   249,     7,   335,
     234,   249,   253,   335,     6,     6,   245,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   351,   335,   335,   335,   335,   335,
     335,   335,   351,   351,   366,   248,   335,   335,   356,   335,
     356,   349,   356,   356,   363,   248,   335,   279,   366,     8,
     335,   335,   350,   349,   356,   356,   351,   342,   357,   342,
     352,   245,   249,   250,   278,    64,     8,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   248,   335,   351,   335,   335,   335,   335,   335,   366,
     335,   335,     4,   343,   248,   279,   245,   249,   249,   335,
     335,   335,     7,     7,   328,   328,   244,   335,   335,     6,
     352,   352,   248,   245,     6,   306,   248,   306,   306,   253,
     253,   253,   346,   346,   305,   305,   253,   335,   249,   319,
     253,   306,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   249,   245,     7,   329,     6,     7,   335,     6,   335,
     335,   249,   352,   352,   352,   335,     6,   335,   335,   335,
     245,   249,   245,   245,   245,   179,   253,   306,   248,     8,
     245,   245,   247,   363,   356,   363,   356,   356,   356,   356,
     356,   356,   335,   356,   356,   356,   356,   251,   359,   366,
     357,   356,   356,   356,   342,   366,   352,   249,   249,   249,
     249,   335,   335,   306,   366,   343,   247,   249,   245,   137,
     154,   323,   245,   249,   253,   335,     6,   248,   349,   245,
     247,     7,   276,   277,   250,     7,     6,   352,     7,   217,
     276,   335,   261,   366,   335,   335,   343,   246,   244,   123,
     300,   301,   300,   248,   249,     6,   228,   231,   258,   352,
     366,   335,   335,     4,   343,     6,   352,   335,     6,   356,
     364,   366,   245,   343,     6,   366,     6,   356,   335,   245,
     246,   335,   253,   253,   253,   253,   357,     7,     7,     7,
     245,     7,     7,     7,   245,     7,     7,     7,     7,     7,
       7,     7,     7,     7,     7,   335,   245,   248,   335,   351,
     249,     6,   279,   279,   279,   279,   279,   279,   279,   279,
     279,   279,   279,   253,   279,   279,   279,   279,   279,   279,
     279,   279,   279,   253,   253,   253,   245,   247,   247,   352,
     253,   253,   279,   253,   279,   253,   253,   253,   245,   352,
     337,   279,   249,   249,   249,   253,   253,   279,   279,   245,
     250,   245,   250,   253,   278,   338,   249,     7,   343,   279,
     248,   249,     8,     8,   352,   250,   245,   247,   277,   244,
     246,   278,   352,     7,   248,   248,   245,   245,   245,   335,
     349,     4,   327,     6,   295,   335,   357,   249,   245,   249,
     249,   352,   250,   249,   306,   249,   249,   346,   335,   335,
     249,   249,   335,   346,   134,   134,   151,   159,   160,   161,
     165,   166,   320,   321,   346,   249,   316,   245,   249,   245,
     245,   245,   245,   245,   245,   245,   248,     7,   335,     6,
     335,   245,   247,   247,   249,   249,   249,   249,   247,   247,
     253,     7,     7,     7,   250,   335,   249,   335,   335,     7,
     250,   279,   253,   279,   279,   245,   245,   253,   279,   279,
     253,   253,   279,   279,   279,   279,   335,   279,     9,   358,
     253,   245,   253,   279,   250,   253,   339,   247,   249,   249,
     250,   244,   246,   252,   182,     7,   154,     6,   335,   249,
     248,     6,   349,   249,   335,     6,     7,   276,   277,   250,
     276,   277,   357,   335,     6,     4,   248,   354,   366,   249,
      46,    46,   349,   249,     4,   169,   170,   171,   172,   249,
     264,   268,   271,   273,   274,   250,   245,   247,   244,   335,
     335,   244,   248,   244,   248,     8,   352,   356,   245,   250,
     245,   247,   244,   245,   253,   250,   244,     7,   278,     4,
     289,   290,   291,   279,   335,   335,   335,   335,   279,   346,
     349,   349,     7,   349,   349,   349,     7,   349,   349,   349,
     349,   349,   349,   349,   349,   349,   349,     6,     7,   352,
     335,   335,   335,   335,   249,   335,   335,   335,   349,   356,
     356,   249,   253,   288,   335,   335,   343,   343,   335,   335,
     245,   349,   278,   335,   335,   335,   249,   343,   277,   250,
     277,   335,   335,   279,   249,   349,   352,   352,     7,     7,
       7,   134,   326,     6,   245,   253,     7,     7,     7,   249,
       4,   249,   253,   253,   253,   249,   249,   112,     4,     6,
     335,   248,     6,   244,     6,   167,     6,   167,   249,   321,
     253,   320,     7,     6,     7,     7,     7,     7,     7,     7,
       7,   305,   349,     6,   248,     6,     6,     6,   100,     7,
       6,     6,   335,   349,   349,   349,     4,   253,     8,     8,
     245,     4,   103,   104,     4,   352,   356,   335,   356,   251,
     253,   292,   356,   356,   343,   356,   245,   253,   343,   248,
     300,   248,     6,   335,     6,   248,   349,   249,   249,   335,
       6,     4,   179,   180,   335,     6,     6,     6,     7,   353,
     355,     6,   246,   279,   278,   278,     6,   265,   244,   244,
     248,   275,     6,   343,   250,   356,   335,   247,   245,   335,
       8,   352,   335,   352,   249,   249,     6,     6,   258,   343,
     250,   335,     6,   335,   343,   245,   248,   335,   357,   279,
      46,   248,   349,   357,   360,   247,   253,     6,   245,   245,
     245,   245,     6,     6,   127,   297,   297,   349,     6,     6,
       6,   349,   134,   182,   296,     6,     6,     6,     6,     6,
       6,     6,     6,     6,     5,   249,   279,   279,   279,   279,
     279,   253,   253,   253,   245,   279,   279,   290,   279,   245,
     279,   245,   278,   338,   279,     6,   279,   253,   244,   246,
     278,     4,   245,   247,   279,     6,   249,   249,   349,   349,
     349,     4,     6,   276,   335,   349,   248,   248,     7,     6,
       7,   335,   335,   335,   248,   248,   248,   246,     6,   335,
     349,   335,     6,     6,   335,   346,   249,     5,   349,   248,
     248,   248,   248,   248,   248,   248,   349,   249,     6,   352,
     248,   335,   247,     6,     6,   178,   335,   335,   335,     6,
       6,     6,     6,     7,   279,   253,   253,   279,   253,   335,
       4,   194,   293,   294,   279,   245,   279,   339,   357,   244,
     246,   335,   248,   306,     6,   306,   253,     6,     6,     7,
     276,   277,   250,     7,     6,   353,   249,   253,   335,   276,
     248,   279,   361,   362,   363,   361,   244,   335,   335,   348,
     349,   248,   244,     4,     6,   245,     6,   245,   249,   249,
     245,   249,     6,     6,   356,   244,     4,   245,   253,   244,
     349,   357,     7,   278,   287,   335,   351,   291,     6,     6,
       6,     6,   346,     6,     6,     6,     6,    94,   113,    98,
       6,     5,   248,   335,   335,   335,   335,   245,   338,   335,
     335,   335,   279,   277,   248,   248,     6,   296,     6,   335,
     349,     4,     6,   352,   352,   335,   335,   357,   249,   245,
     249,   253,   305,   305,   335,   335,   249,   253,   245,   249,
     253,     6,     6,   348,   346,   346,   346,   346,   346,   233,
     346,     6,   249,   335,     6,     6,   349,   249,   253,     8,
     249,   245,   248,   335,   357,   356,   335,   356,   335,   357,
     360,   362,   357,   253,   245,   253,   249,   335,   323,   323,
     349,   357,   335,     6,     4,   354,     6,   353,   247,   349,
     363,     6,   279,   279,   262,   335,   253,   253,   249,   253,
     263,   335,   335,     6,     6,     6,     6,   335,   335,   245,
     283,   285,   248,   362,   249,   253,     7,     7,   248,   248,
     248,     5,   348,   279,   279,   253,   279,   245,   253,   245,
     247,   352,   352,     6,     6,   248,   249,   249,   248,     6,
       6,   248,   335,   249,   249,   249,   247,     6,   349,     7,
     248,   335,   249,   253,   253,   253,   253,   253,   253,     6,
     249,   177,   335,   335,   352,     6,     6,   245,   279,   279,
     294,   357,   249,   249,   249,   249,     6,     6,     7,     6,
     250,     6,   249,     6,     6,   245,   253,   335,   335,   248,
     349,   249,   253,   245,   245,   253,   288,   292,   349,   279,
     335,   357,   366,   352,   352,   335,     6,   249,   335,   338,
     335,   249,   249,   348,   138,   139,   144,   330,   138,   139,
     330,   352,   305,   249,   253,     6,   249,   349,   306,   249,
       6,   352,   346,   346,   346,   346,   346,   335,   249,   249,
     249,   245,     6,   248,     6,   353,   180,   266,   335,   253,
     253,   348,     6,   335,   335,   249,   249,   284,     7,   244,
     249,   249,   249,   248,   253,   245,   253,   248,   249,   248,
     346,   349,     6,   248,   346,     6,   249,   249,   335,     6,
     134,   249,   317,   248,   249,   253,   253,   253,   253,   253,
       6,     6,     6,   306,     6,   248,   335,   335,   249,   253,
     288,   357,   245,   335,   335,   335,   352,     6,   346,     6,
     346,     6,     6,   249,   335,   320,   306,     6,   352,   352,
     352,   352,   346,   352,   323,   263,   245,   253,     6,   248,
     335,   249,   253,   253,   253,   249,   253,   253,     6,   249,
     249,   318,   249,   249,   249,   249,   253,   249,   249,   249,
     269,   335,   348,   249,   335,   335,   335,   346,   346,   320,
       6,     6,     6,     6,   352,     6,     6,     6,   248,   245,
     249,     6,   249,   279,   253,   253,   253,   249,   249,   267,
     356,   272,   248,     6,   335,   335,   335,     6,   249,   253,
     248,   348,   249,   249,   249,     6,   356,   270,   356,   249,
       6,     6,   249,   253,     6,     6,   356
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   254,   255,   255,   256,   256,   257,   257,   257,   257,
     257,   257,   257,   257,   257,   257,   257,   257,   257,   257,
     257,   257,   257,   257,   257,   257,   258,   258,   259,   259,
     259,   259,   259,   259,   260,   260,   260,   260,   261,   261,
     261,   261,   261,   261,   262,   262,   263,   263,   265,   266,
     264,   267,   267,   269,   268,   270,   270,   272,   271,   273,
     273,   275,   274,   276,   276,   276,   276,   276,   277,   277,
     278,   278,   279,   279,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   281,   281,   282,   282,   282,   283,   282,   284,   282,
     282,   285,   282,   286,   286,   287,   287,   287,   288,   288,
     289,   289,   290,   290,   291,   291,   291,   291,   291,   292,
     292,   293,   293,   294,   294,   294,   294,   294,   295,   295,
     295,   296,   296,   296,   296,   297,   297,   298,   298,   298,
     298,   298,   298,   298,   298,   298,   298,   298,   298,   298,
     298,   298,   298,   298,   298,   298,   298,   298,   298,   298,
     298,   298,   298,   298,   298,   298,   298,   298,   298,   299,
     298,   300,   300,   300,   300,   300,   301,   301,   301,   301,
     302,   302,   302,   303,   303,   303,   303,   304,   304,   304,
     304,   304,   304,   304,   304,   304,   305,   305,   306,   306,
     306,   306,   306,   306,   306,   307,   307,   307,   307,   307,
     307,   307,   307,   307,   307,   307,   307,   308,   308,   308,
     308,   308,   308,   308,   309,   309,   310,   311,   311,   311,
     311,   311,   311,   311,   311,   312,   312,   312,   312,   312,
     312,   312,   312,   312,   312,   312,   312,   312,   312,   312,
     312,   312,   312,   312,   312,   312,   313,   314,   314,   314,
     314,   314,   314,   314,   314,   314,   314,   314,   314,   314,
     314,   315,   315,   315,   316,   315,   317,   315,   318,   315,
     319,   315,   315,   315,   315,   315,   315,   320,   320,   321,
     321,   321,   321,   321,   321,   321,   321,   321,   321,   321,
     322,   322,   322,   322,   322,   323,   323,   323,   323,   323,
     324,   324,   325,   326,   326,   327,   327,   328,   328,   329,
     329,   330,   330,   331,   331,   331,   331,   331,   331,   331,
     331,   331,   331,   331,   331,   331,   331,   331,   331,   331,
     331,   331,   331,   331,   331,   331,   331,   331,   331,   332,
     332,   332,   333,   333,   333,   334,   334,   334,   334,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   336,   336,   336,   336,
     336,   336,   336,   336,   336,   336,   336,   337,   336,   336,
     336,   336,   336,   336,   336,   336,   336,   336,   336,   336,
     336,   336,   336,   336,   336,   336,   336,   336,   336,   336,
     336,   336,   336,   336,   336,   336,   336,   336,   336,   336,
     336,   336,   336,   336,   338,   338,   339,   339,   341,   340,
     342,   342,   343,   344,   344,   345,   345,   346,   346,   346,
     346,   346,   347,   347,   347,   347,   348,   348,   349,   349,
     349,   349,   349,   349,   350,   350,   350,   351,   351,   351,
     351,   351,   351,   351,   351,   351,   351,   351,   351,   351,
     351,   351,   351,   351,   351,   351,   351,   351,   351,   351,
     351,   351,   351,   351,   351,   351,   351,   352,   352,   352,
     352,   353,   353,   353,   353,   354,   354,   355,   355,   356,
     356,   356,   356,   356,   356,   356,   356,   356,   356,   356,
     357,   357,   357,   357,   357,   357,   357,   357,   357,   357,
     357,   357,   357,   357,   357,   357,   357,   357,   357,   357,
     357,   357,   357,   357,   357,   357,   357,   357,   357,   358,
     357,   357,   359,   359,   360,   361,   361,   362,   363,   363,
     363,   363,   364,   364,   364,   365,   365,   365,   366,   366,
     366
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     2,     0,     2,     1,     1,     1,     5,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     5,     5,
       7,     7,     7,     9,     6,     6,     6,     8,     0,     2,
       2,     2,     2,     2,     1,     3,     1,     3,     0,     0,
      10,     1,     3,     0,    13,     1,     3,     0,    15,     8,
      14,     0,     6,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     5,     5,     2,     7,     7,     4,
       3,     6,     9,     7,     6,     7,     6,     4,     8,     9,
       9,     6,     9,     6,     9,     5,     8,     8,    11,     6,
       9,     5,     7,     7,     7,     9,     9,    11,     7,     9,
       9,     0,     1,     0,     3,     5,     0,     9,     0,    11,
       5,     0,     9,     0,     3,     3,     5,     5,     0,     2,
       0,     1,     1,     3,     2,     1,     4,     2,     2,     0,
       2,     1,     3,     2,     2,     2,     2,     2,     1,     1,
       3,     0,     5,     5,     5,     0,     2,     7,     7,     7,
       8,     8,     7,     7,    11,     7,     8,     8,     8,     9,
       3,     4,    10,     7,     7,     7,     7,     7,     7,     7,
       7,     7,     7,     8,     7,     7,     8,     8,    12,     0,
       9,     1,     1,     1,     1,     4,     1,     1,     1,     4,
       1,     1,     4,     1,     1,     1,     4,     5,    11,     5,
       9,     9,     7,     4,     9,     9,     1,     1,     0,     2,
       6,     7,     7,     6,     7,     8,    10,    14,    16,    12,
      14,    14,    14,    14,     8,     8,     6,     4,     5,     6,
       6,     3,     4,     3,     5,     6,     5,     4,     3,     4,
       3,     4,     5,     4,     5,     3,     5,     7,     7,     3,
       7,     3,     2,     2,     2,     2,     2,    15,     2,     2,
       2,     2,     2,     2,     2,    16,    11,     6,     8,     8,
      10,     1,     2,     2,     1,     3,     3,     4,     4,     1,
       1,     5,    11,    13,     0,     7,     0,    13,     0,    15,
       0,     6,     9,     2,     3,    10,    13,     1,     2,     5,
       7,     2,     2,     3,     2,     3,     2,     3,     9,     6,
       1,     1,     1,     1,     1,     0,     2,     3,     3,     4,
       9,     4,    14,     0,     3,     0,     1,     0,     2,     0,
       2,     0,     2,     6,     7,     6,     5,     3,     8,     8,
       8,     8,     8,     5,     4,     6,    11,    11,    18,    18,
      12,    12,    12,    10,     4,     4,     4,     4,     4,     2,
       3,     6,     1,     1,     1,     2,     5,     7,    10,     1,
       3,     2,     2,     2,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     5,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     6,     4,     4,     4,     4,     4,     4,
       4,     4,     6,     6,     6,     4,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     6,     1,
       4,     6,     1,     4,     4,     4,     6,     5,     7,     8,
      10,     4,     4,     6,     4,     3,     2,     5,     5,     3,
       5,     6,     8,     6,     8,     6,     4,     7,     6,     6,
       6,     4,     6,     4,     0,     2,     0,     2,     0,     7,
       1,     3,     1,     1,     2,     0,     3,     1,     2,     2,
       3,     3,    11,     9,     7,     7,     1,     3,     1,     1,
       2,     3,     4,     5,     1,     3,     1,     2,     3,     3,
       5,     4,     4,     2,     4,     2,     3,     3,    16,     5,
       1,     1,     1,     3,     5,     7,     4,     4,     4,     6,
       6,     8,     8,     4,    14,     4,     4,     1,     1,     3,
       3,     9,     7,     1,     5,     3,     6,     1,     3,     1,
       1,     4,     4,     3,     5,     6,     8,     6,     4,     5,
       1,     4,     1,     1,     1,     1,     4,     6,     4,     6,
       5,     7,     4,     4,     4,     8,     4,     4,     4,     4,
       8,     8,     6,     4,     6,     4,     1,     4,     4,     0,
       6,     4,     2,     4,     4,     1,     1,     3,     1,     1,
       3,     3,     3,     5,     7,     5,     5,     8,     1,     1,
       4
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyo, yytype, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &yyvsp[(yyi + 1) - (yynrhs)]
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            else
              goto append;

          append:
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return (YYSIZE_T) (yystpcpy (yyres, yystr) - yyres);
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
                    yysize = yysize1;
                  else
                    return 2;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
    default: /* Avoid compiler warnings. */
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
      yysize = yysize1;
    else
      return 2;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yynewstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  *yyssp = (yytype_int16) yystate;

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = (YYSIZE_T) (yyssp - yyss + 1);

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
# undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 269 "Gmsh.y" /* yacc.c:1652  */
    { yyerrok; return 1; }
#line 5913 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 6:
#line 280 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5919 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 7:
#line 281 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5925 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 8:
#line 282 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5931 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 9:
#line 284 "Gmsh.y" /* yacc.c:1652  */
    {
      gmsh_yyfactory = (yyvsp[-2].c);
      if(gmsh_yyfactory == "OpenCASCADE"){
        if(!GModel::current()->getOCCInternals())
          GModel::current()->createOCCInternals();
        for(int dim = -2; dim <= 3; dim++)
          GModel::current()->getOCCInternals()->setMaxTag
            (dim, std::max(GModel::current()->getOCCInternals()->getMaxTag(dim),
                           GModel::current()->getGEOInternals()->getMaxTag(dim)));
      }
      else if(GModel::current()->getOCCInternals()){
        for(int dim = -2; dim <= 3; dim++)
          GModel::current()->getGEOInternals()->setMaxTag
            (dim, std::max(GModel::current()->getGEOInternals()->getMaxTag(dim),
                           GModel::current()->getOCCInternals()->getMaxTag(dim)));
      }
      Free((yyvsp[-2].c));
    }
#line 5954 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 10:
#line 302 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5960 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 11:
#line 303 "Gmsh.y" /* yacc.c:1652  */
    { List_Delete((yyvsp[0].l)); return 1; }
#line 5966 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 12:
#line 304 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5972 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 13:
#line 305 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5978 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 14:
#line 306 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5984 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 15:
#line 307 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 5990 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 16:
#line 308 "Gmsh.y" /* yacc.c:1652  */
    { List_Delete((yyvsp[0].l)); return 1; }
#line 5996 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 17:
#line 309 "Gmsh.y" /* yacc.c:1652  */
    { List_Delete((yyvsp[0].l)); return 1; }
#line 6002 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 18:
#line 310 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6008 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 19:
#line 311 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6014 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 20:
#line 312 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6020 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 21:
#line 313 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6026 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 22:
#line 314 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6032 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 23:
#line 315 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6038 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 24:
#line 316 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6044 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 25:
#line 317 "Gmsh.y" /* yacc.c:1652  */
    { return 1; }
#line 6050 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 26:
#line 322 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (char*)"w";
    }
#line 6058 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 27:
#line 326 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (char*)"a";
    }
#line 6066 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 28:
#line 333 "Gmsh.y" /* yacc.c:1652  */
    {
      Msg::Direct((yyvsp[-2].c));
      Free((yyvsp[-2].c));
    }
#line 6075 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 29:
#line 338 "Gmsh.y" /* yacc.c:1652  */
    {
      Msg::Error((yyvsp[-2].c));
      Free((yyvsp[-2].c));
    }
#line 6084 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 30:
#line 343 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
      FILE *fp = Fopen(tmp.c_str(), (yyvsp[-2].c));
      if(!fp){
	yymsg(0, "Unable to open file '%s'", tmp.c_str());
      }
      else{
	fprintf(fp, "%s\n", (yyvsp[-4].c));
	fclose(fp);
      }
      Free((yyvsp[-4].c));
      Free((yyvsp[-1].c));
    }
#line 6102 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 31:
#line 357 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstring[5000];
      int i = printListOfDouble((yyvsp[-4].c), (yyvsp[-2].l), tmpstring);
      if(i < 0)
	yymsg(0, "Too few arguments in Printf");
      else if(i > 0)
	yymsg(0, "%d extra argument%s in Printf", i, (i > 1) ? "s" : "");
      else
	Msg::Direct(tmpstring);
      Free((yyvsp[-4].c));
      List_Delete((yyvsp[-2].l));
    }
#line 6119 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 32:
#line 370 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstring[5000];
      int i = printListOfDouble((yyvsp[-4].c), (yyvsp[-2].l), tmpstring);
      if(i < 0)
	yymsg(0, "Too few arguments in Error");
      else if(i > 0)
	yymsg(0, "%d extra argument%s in Error", i, (i > 1) ? "s" : "");
      else
	Msg::Error(tmpstring);
      Free((yyvsp[-4].c));
      List_Delete((yyvsp[-2].l));
    }
#line 6136 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 33:
#line 383 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstring[5000];
      int i = printListOfDouble((yyvsp[-6].c), (yyvsp[-4].l), tmpstring);
      if(i < 0)
	yymsg(0, "Too few arguments in Printf");
      else if(i > 0)
	yymsg(0, "%d extra argument%s in Printf", i, (i > 1) ? "s" : "");
      else{
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
	FILE *fp = Fopen(tmp.c_str(), (yyvsp[-2].c));
	if(!fp){
	  yymsg(0, "Unable to open file '%s'", tmp.c_str());
	}
	else{
	  fprintf(fp, "%s\n", tmpstring);
	  fclose(fp);
	}
      }
      Free((yyvsp[-6].c));
      Free((yyvsp[-1].c));
      List_Delete((yyvsp[-4].l));
    }
#line 6163 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 34:
#line 411 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-5].c), "View") && ViewData->finalize()){
	ViewData->setName((yyvsp[-4].c));
	ViewData->setFileName(gmsh_yyname);
	ViewData->setFileIndex(gmsh_yyviewindex++);
	new PView(ViewData);
      }
      else
	delete ViewData;
#endif
      Free((yyvsp[-5].c)); Free((yyvsp[-4].c));
    }
#line 6181 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 35:
#line 425 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-4].c), "View")){
	int index = (int)(yyvsp[-2].d);
	if(index >= 0 && index < (int)PView::list.size())
	  new PView(PView::list[index], false);
        else
	  yymsg(0, "Unknown view %d", index);
      }
#endif
      Free((yyvsp[-4].c));
    }
#line 6198 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 36:
#line 438 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-4].c), "View")){
	int index = (int)(yyvsp[-2].d);
	if(index >= 0 && index < (int)PView::list.size())
	  new PView(PView::list[index], true);
        else
	  yymsg(0, "Unknown view %d", index);
      }
#endif
      Free((yyvsp[-4].c));
    }
#line 6215 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 37:
#line 451 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-6].c), "View")){
	int index = (int)(yyvsp[-4].d), index2 = (int)(yyvsp[-2].d);
	if(index >= 0 && index < (int)PView::list.size() &&
           index2 >= 0 && index2 < (int)PView::list.size()){
          PView::list[index2]->setOptions(PView::list[index]->getOptions());
        }
        else
	  yymsg(0, "Unknown view %d or %d", index, index2);
      }
#endif
      Free((yyvsp[-6].c));
    }
#line 6234 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 38:
#line 469 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      ViewData = new PViewDataList();
#endif
    }
#line 6244 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 44:
#line 483 "Gmsh.y" /* yacc.c:1652  */
    { ViewCoord.push_back((yyvsp[0].d)); }
#line 6250 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 45:
#line 485 "Gmsh.y" /* yacc.c:1652  */
    { ViewCoord.push_back((yyvsp[0].d)); }
#line 6256 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 46:
#line 490 "Gmsh.y" /* yacc.c:1652  */
    { if(ViewValueList) ViewValueList->push_back((yyvsp[0].d)); }
#line 6262 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 47:
#line 492 "Gmsh.y" /* yacc.c:1652  */
    { if(ViewValueList) ViewValueList->push_back((yyvsp[0].d)); }
#line 6268 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 48:
#line 497 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strncmp((yyvsp[0].c), "SP", 2)){
	ViewValueList = &ViewData->SP; ViewNumList = &ViewData->NbSP;
      }
      else if(!strncmp((yyvsp[0].c), "VP", 2)){
	ViewValueList = &ViewData->VP; ViewNumList = &ViewData->NbVP;
      }
      else if(!strncmp((yyvsp[0].c), "TP", 2)){
	ViewValueList = &ViewData->TP; ViewNumList = &ViewData->NbTP;
      }
      else if(!strncmp((yyvsp[0].c), "SL", 2)){
	ViewValueList = &ViewData->SL; ViewNumList = &ViewData->NbSL;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_LIN);
      }
      else if(!strncmp((yyvsp[0].c), "VL", 2)){
	ViewValueList = &ViewData->VL; ViewNumList = &ViewData->NbVL;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_LIN);
      }
      else if(!strncmp((yyvsp[0].c), "TL", 2)){
	ViewValueList = &ViewData->TL; ViewNumList = &ViewData->NbTL;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_LIN);
      }
      else if(!strncmp((yyvsp[0].c), "ST", 2)){
	ViewValueList = &ViewData->ST; ViewNumList = &ViewData->NbST;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_TRI);
      }
      else if(!strncmp((yyvsp[0].c), "VT", 2)){
	ViewValueList = &ViewData->VT; ViewNumList = &ViewData->NbVT;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_TRI);
      }
      else if(!strncmp((yyvsp[0].c), "TT", 2)){
	ViewValueList = &ViewData->TT; ViewNumList = &ViewData->NbTT;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_TRI);
      }
      else if(!strncmp((yyvsp[0].c), "SQ", 2)){
	ViewValueList = &ViewData->SQ; ViewNumList = &ViewData->NbSQ;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_QUA);
      }
      else if(!strncmp((yyvsp[0].c), "VQ", 2)){
	ViewValueList = &ViewData->VQ; ViewNumList = &ViewData->NbVQ;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_QUA);
      }
      else if(!strncmp((yyvsp[0].c), "TQ", 2)){
	ViewValueList = &ViewData->TQ; ViewNumList = &ViewData->NbTQ;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_QUA);
      }
      else if(!strncmp((yyvsp[0].c), "SS", 2)){
	ViewValueList = &ViewData->SS; ViewNumList = &ViewData->NbSS;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_TET);
      }
      else if(!strncmp((yyvsp[0].c), "VS", 2)){
	ViewValueList = &ViewData->VS; ViewNumList = &ViewData->NbVS;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_TET);
      }
      else if(!strncmp((yyvsp[0].c), "TS", 2)){
	ViewValueList = &ViewData->TS; ViewNumList = &ViewData->NbTS;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_TET);
      }
      else if(!strncmp((yyvsp[0].c), "SH", 2)){
	ViewValueList = &ViewData->SH; ViewNumList = &ViewData->NbSH;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_HEX);
      }
      else if(!strncmp((yyvsp[0].c), "VH", 2)){
	ViewValueList = &ViewData->VH; ViewNumList = &ViewData->NbVH;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_HEX);
      }
      else if(!strncmp((yyvsp[0].c), "TH", 2)){
	ViewValueList = &ViewData->TH; ViewNumList = &ViewData->NbTH;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_HEX);
      }
      else if(!strncmp((yyvsp[0].c), "SI", 2)){
	ViewValueList = &ViewData->SI; ViewNumList = &ViewData->NbSI;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_PRI);
      }
      else if(!strncmp((yyvsp[0].c), "VI", 2)){
	ViewValueList = &ViewData->VI; ViewNumList = &ViewData->NbVI;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_PRI);
      }
      else if(!strncmp((yyvsp[0].c), "TI", 2)){
	ViewValueList = &ViewData->TI; ViewNumList = &ViewData->NbTI;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_PRI);
      }
      else if(!strncmp((yyvsp[0].c), "SY", 2)){
	ViewValueList = &ViewData->SY; ViewNumList = &ViewData->NbSY;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_PYR);
      }
      else if(!strncmp((yyvsp[0].c), "VY", 2)){
	ViewValueList = &ViewData->VY; ViewNumList = &ViewData->NbVY;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_PYR);
      }
      else if(!strncmp((yyvsp[0].c), "TY", 2)){
	ViewValueList = &ViewData->TY; ViewNumList = &ViewData->NbTY;
        if(strlen((yyvsp[0].c)) > 2) ViewData->setOrder2(TYPE_PYR);
      }
      else{
	yymsg(0, "Unknown element type '%s'", (yyvsp[0].c));
	ViewValueList = 0; ViewNumList = 0;
      }
#endif
      ViewCoord.clear();
      Free((yyvsp[0].c));
    }
#line 6376 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 49:
#line 601 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(ViewValueList){
	for(int i = 0; i < 3; i++)
	  for(std::size_t j = 0; j < ViewCoord.size() / 3; j++)
	    ViewValueList->push_back(ViewCoord[3 * j + i]);
      }
#endif
    }
#line 6390 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 50:
#line 611 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(ViewValueList) (*ViewNumList)++;
#endif
    }
#line 6400 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 51:
#line 620 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      for(int i = 0; i < (int)strlen((yyvsp[0].c)) + 1; i++) ViewData->T2C.push_back((yyvsp[0].c)[i]);
#endif
      Free((yyvsp[0].c));
    }
#line 6411 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 52:
#line 627 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      for(int i = 0; i < (int)strlen((yyvsp[0].c)) + 1; i++) ViewData->T2C.push_back((yyvsp[0].c)[i]);
#endif
      Free((yyvsp[0].c));
    }
#line 6422 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 53:
#line 637 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      ViewData->T2D.push_back((yyvsp[-5].d));
      ViewData->T2D.push_back((yyvsp[-3].d));
      ViewData->T2D.push_back((yyvsp[-1].d));
      ViewData->T2D.push_back(ViewData->T2C.size());
#endif
    }
#line 6435 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 54:
#line 646 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      ViewData->NbT2++;
#endif
    }
#line 6445 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 55:
#line 655 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      for(int i = 0; i < (int)strlen((yyvsp[0].c)) + 1; i++) ViewData->T3C.push_back((yyvsp[0].c)[i]);
#endif
      Free((yyvsp[0].c));
    }
#line 6456 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 56:
#line 662 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      for(int i = 0; i < (int)strlen((yyvsp[0].c)) + 1; i++) ViewData->T3C.push_back((yyvsp[0].c)[i]);
#endif
      Free((yyvsp[0].c));
    }
#line 6467 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 57:
#line 672 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      ViewData->T3D.push_back((yyvsp[-7].d)); ViewData->T3D.push_back((yyvsp[-5].d));
      ViewData->T3D.push_back((yyvsp[-3].d)); ViewData->T3D.push_back((yyvsp[-1].d));
      ViewData->T3D.push_back(ViewData->T3C.size());
#endif
    }
#line 6479 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 58:
#line 680 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      ViewData->NbT3++;
#endif
    }
#line 6489 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 59:
#line 690 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      int type =
	(ViewData->NbSL || ViewData->NbVL) ? TYPE_LIN :
	(ViewData->NbST || ViewData->NbVT) ? TYPE_TRI :
	(ViewData->NbSQ || ViewData->NbVQ) ? TYPE_QUA :
	(ViewData->NbSS || ViewData->NbVS) ? TYPE_TET :
	(ViewData->NbSY || ViewData->NbVY) ? TYPE_PYR :
	(ViewData->NbSI || ViewData->NbVI) ? TYPE_PRI :
      	(ViewData->NbSH || ViewData->NbVH) ? TYPE_HEX :
	0;
      ViewData->setInterpolationMatrices(type, ListOfListOfDouble2Matrix((yyvsp[-5].l)),
                                         ListOfListOfDouble2Matrix((yyvsp[-2].l)));
#endif
    }
#line 6509 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 60:
#line 709 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      int type =
	(ViewData->NbSL || ViewData->NbVL) ? TYPE_LIN :
	(ViewData->NbST || ViewData->NbVT) ? TYPE_TRI :
	(ViewData->NbSQ || ViewData->NbVQ) ? TYPE_QUA :
	(ViewData->NbSS || ViewData->NbVS) ? TYPE_TET :
      	(ViewData->NbSH || ViewData->NbVH) ? TYPE_HEX :
	0;
      ViewData->setInterpolationMatrices(type, ListOfListOfDouble2Matrix((yyvsp[-11].l)),
                                         ListOfListOfDouble2Matrix((yyvsp[-8].l)),
                                         ListOfListOfDouble2Matrix((yyvsp[-5].l)),
                                         ListOfListOfDouble2Matrix((yyvsp[-2].l)));
#endif
    }
#line 6529 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 61:
#line 728 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      ViewValueList = &ViewData->Time;
#endif
    }
#line 6539 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 62:
#line 734 "Gmsh.y" /* yacc.c:1652  */
    {
    }
#line 6546 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 63:
#line 741 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 0; }
#line 6552 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 64:
#line 742 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 6558 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 65:
#line 743 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 2; }
#line 6564 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 66:
#line 744 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 3; }
#line 6570 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 67:
#line 745 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 4; }
#line 6576 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 68:
#line 749 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 6582 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 69:
#line 750 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = -1; }
#line 6588 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 70:
#line 756 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)"("; }
#line 6594 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 71:
#line 756 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)"["; }
#line 6600 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 72:
#line 757 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)")"; }
#line 6606 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 73:
#line 757 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)"]"; }
#line 6612 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 77:
#line 767 "Gmsh.y" /* yacc.c:1652  */
    {
      Msg::SetOnelabNumber((yyvsp[-4].c), (yyvsp[-2].d));
      Free((yyvsp[-4].c));
    }
#line 6621 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 78:
#line 772 "Gmsh.y" /* yacc.c:1652  */
    {
      Msg::SetOnelabString((yyvsp[-4].c), (yyvsp[-2].c));
      Free((yyvsp[-4].c));
      Free((yyvsp[-2].c));
    }
#line 6631 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 79:
#line 778 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!gmsh_yysymbols.count((yyvsp[-3].c)) && (yyvsp[-2].i) && List_Nbr((yyvsp[-1].l)) == 1){
        yymsg(0, "Unknown variable '%s'", (yyvsp[-3].c));
      }
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-3].c)]);
        if(!(yyvsp[-2].i)) s.list = (List_Nbr((yyvsp[-1].l)) != 1); // list if 0 or > 1 elements
        if(!s.list){ // single expression
          if(List_Nbr((yyvsp[-1].l)) != 1){
            yymsg(0, "Cannot assign list to variable '%s'", (yyvsp[-3].c));
          }
          else{
            double d;
            List_Read((yyvsp[-1].l), 0, &d);
            if(s.value.empty()){
              if((yyvsp[-2].i)) yymsg(1, "Uninitialized variable '%s'", (yyvsp[-3].c));
              s.value.resize(1, 0.);
            }
            switch((yyvsp[-2].i)){
            case 0 : s.value[0] = d; break;
            case 1 : s.value[0] += d; break;
            case 2 : s.value[0] -= d; break;
            case 3 : s.value[0] *= d; break;
            case 4 :
              if(d) s.value[0] /= d;
              else yymsg(0, "Division by zero in '%s /= %g'", (yyvsp[-3].c), d);
              break;
            }
          }
        }
        else{
          // list of expressions; this is not recommended (should use [] or ()
          // notation instead)
          switch((yyvsp[-2].i)){
          case 0: // affect
            s.value.clear(); // fall-through
          case 1: // append
            for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
              double d;
              List_Read((yyvsp[-1].l), i, &d);
              s.value.push_back(d);
            }
            break;
          case 2: // remove
            for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
              double d;
              List_Read((yyvsp[-1].l), i, &d);
              std::vector<double>::iterator it = std::find(s.value.begin(),
                                                           s.value.end(), d);
              if(it != s.value.end()) s.value.erase(it);
            }
            break;
          default:
            yymsg(0, "Operators *= and /= not available for lists");
            break;
          }
        }
      }
      Free((yyvsp[-3].c));
      List_Delete((yyvsp[-1].l));
    }
#line 6697 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 80:
#line 840 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!gmsh_yysymbols.count((yyvsp[-2].c)))
	yymsg(0, "Unknown variable '%s'", (yyvsp[-2].c));
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-2].c)]);
        if(!s.list && s.value.empty())
          yymsg(0, "Uninitialized variable '%s'", (yyvsp[-2].c));
        else if(!s.list)
          s.value[0] += (yyvsp[-1].i);
        else
          yymsg(0, "Variable '%s' is a list", (yyvsp[-2].c));
      }
      Free((yyvsp[-2].c));
    }
#line 6716 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 81:
#line 855 "Gmsh.y" /* yacc.c:1652  */
    {
      gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-5].c)]);
      s.list = true;
      double d;
      switch((yyvsp[-2].i)){
      case 0: // affect
        s.value.clear(); // fall-through
      case 1: // append
        for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
          List_Read((yyvsp[-1].l), i, &d);
          s.value.push_back(d);
        }
        break;
      case 2: // remove
        for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
          List_Read((yyvsp[-1].l), i, &d);
          std::vector<double>::iterator it = std::find(s.value.begin(),
                                                       s.value.end(), d);
          if(it != s.value.end()) s.value.erase(it);
        }
        break;
      default:
        yymsg(0, "Operators *= and /= not available for lists");
        break;
      }
      Free((yyvsp[-5].c));
      List_Delete((yyvsp[-1].l));
    }
#line 6749 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 82:
#line 884 "Gmsh.y" /* yacc.c:1652  */
    {
      assignVariables((yyvsp[-8].c), (yyvsp[-5].l), (yyvsp[-2].i), (yyvsp[-1].l));
      Free((yyvsp[-8].c));
      List_Delete((yyvsp[-5].l));
      List_Delete((yyvsp[-1].l));
    }
#line 6760 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 83:
#line 894 "Gmsh.y" /* yacc.c:1652  */
    {
      assignVariable((yyvsp[-6].c), (int)(yyvsp[-4].d), (yyvsp[-2].i), (yyvsp[-1].d));
      Free((yyvsp[-6].c));
    }
#line 6769 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 84:
#line 899 "Gmsh.y" /* yacc.c:1652  */
    {
      incrementVariable((yyvsp[-5].c), (int)(yyvsp[-3].d), (yyvsp[-1].i));
      Free((yyvsp[-5].c));
    }
#line 6778 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 85:
#line 907 "Gmsh.y" /* yacc.c:1652  */
    {
      assignVariable((yyvsp[-6].c), (int)(yyvsp[-4].d), (yyvsp[-2].i), (yyvsp[-1].d));
      Free((yyvsp[-6].c));
    }
#line 6787 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 86:
#line 912 "Gmsh.y" /* yacc.c:1652  */
    {
      incrementVariable((yyvsp[-5].c), (yyvsp[-3].d), (yyvsp[-1].i));
      Free((yyvsp[-5].c));
    }
#line 6796 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 87:
#line 920 "Gmsh.y" /* yacc.c:1652  */
    {
      gmsh_yystringsymbols[(yyvsp[-3].c)] = std::vector<std::string>(1, (yyvsp[-1].c));
      Free((yyvsp[-3].c));
      Free((yyvsp[-1].c));
    }
#line 6806 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 88:
#line 929 "Gmsh.y" /* yacc.c:1652  */
    {
      gmsh_yystringsymbols[(yyvsp[-7].c)] = std::vector<std::string>();
      Free((yyvsp[-7].c));
    }
#line 6815 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 89:
#line 934 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::string> s;
      for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
        char **c = (char**)List_Pointer((yyvsp[-2].l), i);
        s.push_back(*c);
        Free(*c);
      }
      gmsh_yystringsymbols[(yyvsp[-8].c)] = s;
      Free((yyvsp[-8].c));
      List_Delete((yyvsp[-2].l));
    }
#line 6831 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 90:
#line 946 "Gmsh.y" /* yacc.c:1652  */
    {
      if(gmsh_yystringsymbols.count((yyvsp[-8].c))){
        for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
          char **c = (char**)List_Pointer((yyvsp[-2].l), i);
          gmsh_yystringsymbols[(yyvsp[-8].c)].push_back(*c);
          Free(*c);
        }
      }
      else
        yymsg(0, "Uninitialized variable '%s'", (yyvsp[-8].c));
      Free((yyvsp[-8].c));
      List_Delete((yyvsp[-2].l));
    }
#line 6849 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 91:
#line 963 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp((yyvsp[-1].c));
      StringOption(GMSH_SET|GMSH_GUI, (yyvsp[-5].c), 0, (yyvsp[-3].c), tmp);
      Free((yyvsp[-5].c)); Free((yyvsp[-3].c)); Free((yyvsp[-1].c));
    }
#line 6859 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 92:
#line 969 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp((yyvsp[-1].c));
      StringOption(GMSH_SET|GMSH_GUI, (yyvsp[-8].c), (int)(yyvsp[-6].d), (yyvsp[-3].c), tmp);
      Free((yyvsp[-8].c)); Free((yyvsp[-3].c)); Free((yyvsp[-1].c));
    }
#line 6869 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 93:
#line 978 "Gmsh.y" /* yacc.c:1652  */
    {
      double d = 0.;
      if(NumberOption(GMSH_GET, (yyvsp[-5].c), 0, (yyvsp[-3].c), d)){
	switch((yyvsp[-2].i)){
	case 0 : d = (yyvsp[-1].d); break;
	case 1 : d += (yyvsp[-1].d); break;
	case 2 : d -= (yyvsp[-1].d); break;
	case 3 : d *= (yyvsp[-1].d); break;
	case 4 :
	  if((yyvsp[-1].d)) d /= (yyvsp[-1].d);
	  else yymsg(0, "Division by zero in '%s.%s /= %g'", (yyvsp[-5].c), (yyvsp[-3].c), (yyvsp[-1].d));
	  break;
	}
	NumberOption(GMSH_SET|GMSH_GUI, (yyvsp[-5].c), 0, (yyvsp[-3].c), d);
      }
      Free((yyvsp[-5].c)); Free((yyvsp[-3].c));
    }
#line 6891 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 94:
#line 996 "Gmsh.y" /* yacc.c:1652  */
    {
      double d = 0.;
      if(NumberOption(GMSH_GET, (yyvsp[-8].c), (int)(yyvsp[-6].d), (yyvsp[-3].c), d)){
	switch((yyvsp[-2].i)){
	case 0 : d = (yyvsp[-1].d); break;
	case 1 : d += (yyvsp[-1].d); break;
	case 2 : d -= (yyvsp[-1].d); break;
	case 3 : d *= (yyvsp[-1].d); break;
	case 4 :
	  if((yyvsp[-1].d)) d /= (yyvsp[-1].d);
	  else yymsg(0, "Division by zero in '%s[%d].%s /= %g'", (yyvsp[-8].c), (int)(yyvsp[-6].d), (yyvsp[-3].c), (yyvsp[-1].d));
	  break;
	}
	NumberOption(GMSH_SET|GMSH_GUI, (yyvsp[-8].c), (int)(yyvsp[-6].d), (yyvsp[-3].c), d);
      }
      Free((yyvsp[-8].c)); Free((yyvsp[-3].c));
    }
#line 6913 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 95:
#line 1014 "Gmsh.y" /* yacc.c:1652  */
    {
      double d = 0.;
      if(NumberOption(GMSH_GET, (yyvsp[-4].c), 0, (yyvsp[-2].c), d)){
	d += (yyvsp[-1].i);
	NumberOption(GMSH_SET|GMSH_GUI, (yyvsp[-4].c), 0, (yyvsp[-2].c), d);
      }
      Free((yyvsp[-4].c)); Free((yyvsp[-2].c));
    }
#line 6926 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 96:
#line 1023 "Gmsh.y" /* yacc.c:1652  */
    {
      double d = 0.;
      if(NumberOption(GMSH_GET, (yyvsp[-7].c), (int)(yyvsp[-5].d), (yyvsp[-2].c), d)){
	d += (yyvsp[-1].i);
	NumberOption(GMSH_SET|GMSH_GUI, (yyvsp[-7].c), (int)(yyvsp[-5].d), (yyvsp[-2].c), d);
      }
      Free((yyvsp[-7].c)); Free((yyvsp[-2].c));
    }
#line 6939 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 97:
#line 1035 "Gmsh.y" /* yacc.c:1652  */
    {
      ColorOption(GMSH_SET|GMSH_GUI, (yyvsp[-7].c), 0, (yyvsp[-3].c), (yyvsp[-1].u));
      Free((yyvsp[-7].c)); Free((yyvsp[-3].c));
    }
#line 6948 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 98:
#line 1040 "Gmsh.y" /* yacc.c:1652  */
    {
      ColorOption(GMSH_SET|GMSH_GUI, (yyvsp[-10].c), (int)(yyvsp[-8].d), (yyvsp[-3].c), (yyvsp[-1].u));
      Free((yyvsp[-10].c)); Free((yyvsp[-3].c));
    }
#line 6957 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 99:
#line 1048 "Gmsh.y" /* yacc.c:1652  */
    {
      GmshColorTable *ct = GetColorTable(0);
      if(!ct)
	yymsg(0, "View[%d] does not exist", 0);
      else{
	ct->size = List_Nbr((yyvsp[-1].l));
	if(ct->size > COLORTABLE_NBMAX_COLOR)
	  yymsg(0, "Too many (%d>%d) colors in View[%d].ColorTable",
		ct->size, COLORTABLE_NBMAX_COLOR, 0);
	else
	  for(int i = 0; i < ct->size; i++) List_Read((yyvsp[-1].l), i, &ct->table[i]);
	if(ct->size == 1){
	  ct->size = 2;
	  ct->table[1] = ct->table[0];
	}
      }
      Free((yyvsp[-5].c));
      List_Delete((yyvsp[-1].l));
    }
#line 6981 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 100:
#line 1068 "Gmsh.y" /* yacc.c:1652  */
    {
      GmshColorTable *ct = GetColorTable((int)(yyvsp[-6].d));
      if(!ct)
	yymsg(0, "View[%d] does not exist", (int)(yyvsp[-6].d));
      else{
	ct->size = List_Nbr((yyvsp[-1].l));
	if(ct->size > COLORTABLE_NBMAX_COLOR)
	  yymsg(0, "Too many (%d>%d) colors in View[%d].ColorTable",
		   ct->size, COLORTABLE_NBMAX_COLOR, (int)(yyvsp[-6].d));
	else
	  for(int i = 0; i < ct->size; i++) List_Read((yyvsp[-1].l), i, &ct->table[i]);
	if(ct->size == 1){
	  ct->size = 2;
	  ct->table[1] = ct->table[0];
	}
      }
      Free((yyvsp[-8].c));
      List_Delete((yyvsp[-1].l));
    }
#line 7005 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 101:
#line 1091 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      if(!strcmp((yyvsp[-4].c),"Background")) {
	if (tags.size() > 1)
	  yymsg(0, "Only 1 field can be set as a background field.");
	else if (tags.size() == 0)
	  yymsg(1, "No field given (Background Field).");
	else
	  GModel::current()->getFields()->setBackgroundFieldId((int)tags[0]);
	  }
      else if(!strcmp((yyvsp[-4].c),"BoundaryLayer"))
	GModel::current()->getFields()->addBoundaryLayerFieldId(tags);
      else
	yymsg(0, "Unknown command '%s Field'", (yyvsp[-4].c));
#endif
      Free((yyvsp[-4].c));
      List_Delete((yyvsp[-1].l));
    }
#line 7029 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 102:
#line 1111 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      if(!GModel::current()->getFields()->newField((int)(yyvsp[-4].d), (yyvsp[-1].c)))
	yymsg(0, "Cannot create field %i of type '%s'", (int)(yyvsp[-4].d), (yyvsp[-1].c));
#endif
      Free((yyvsp[-1].c));
    }
#line 7041 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 103:
#line 1119 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      if(!GModel::current()->getFields()->newField((int)(yyvsp[-4].d), "Box"))
	yymsg(0, "Cannot create field %i of type '%s'", (int)(yyvsp[-4].d), "Box");
#endif
    }
#line 7052 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 104:
#line 1126 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      if(!GModel::current()->getFields()->newField((int)(yyvsp[-4].d), "Cylinder"))
	yymsg(0, "Cannot create field %i of type '%s'", (int)(yyvsp[-4].d), "Cylinder");
#endif
    }
#line 7063 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 105:
#line 1133 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      Field *field = GModel::current()->getFields()->get((int)(yyvsp[-6].d));
      if(field){
	FieldOption *option = field->options[(yyvsp[-3].c)];
	if(option){
	  try { option->numericalValue((yyvsp[-1].d)); }
	  catch(...){
	    yymsg(0, "Cannot assign a numerical value to option '%s' "
		  "in field %i of type '%s'", (yyvsp[-3].c), (int)(yyvsp[-6].d), field->getName());
	  }
	}
	else
	  yymsg(0, "Unknown option '%s' in field %i of type '%s'",
		(yyvsp[-3].c), (int)(yyvsp[-6].d), field->getName());
      }
      else
	yymsg(0, "No field with id %i", (int)(yyvsp[-6].d));
#endif
      Free((yyvsp[-3].c));
    }
#line 7089 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 106:
#line 1155 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      Field *field = GModel::current()->getFields()->get((int)(yyvsp[-6].d));
      if(field){
	FieldOption *option = field->options[(yyvsp[-3].c)];
	if(option){
	  try { option->string((yyvsp[-1].c)); }
	  catch (...){
	    yymsg(0, "Cannot assign a string value to  option '%s' "
		  "in field %i of type '%s'", (yyvsp[-3].c), (int)(yyvsp[-6].d), field->getName());
	  }
	}
	else
	  yymsg(0, "Unknown option '%s' in field %i of type '%s'",
		(yyvsp[-3].c), (int)(yyvsp[-6].d), field->getName());
      }
      else
	yymsg(0, "No field with id %i", (int)(yyvsp[-6].d));
#endif
      Free((yyvsp[-3].c));
      Free((yyvsp[-1].c));
    }
#line 7116 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 107:
#line 1178 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      Field *field = GModel::current()->getFields()->get((int)(yyvsp[-8].d));
      if(field){
	FieldOption *option = field->options[(yyvsp[-5].c)];
	if(option){
	  if (option->getType() == FIELD_OPTION_LIST) {
	    std::list<int> vl = option->list();
	    vl.clear();
	    for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
	      double id;
	      List_Read((yyvsp[-2].l), i, &id);
	      vl.push_back((int)id);
	    }
	    option->list(vl);
	  }
	  else {
	    std::list<double> vl = option->listdouble();
	    vl.clear();
	    for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
	      double id;
	      List_Read((yyvsp[-2].l), i, &id);
	      vl.push_back(id);
	    }
	    option->listdouble(vl);
	  }
	}
	else
	  yymsg(0, "Unknown option '%s' in field %i of type '%s'",
		(yyvsp[-5].c), (int)(yyvsp[-8].d), field->getName());
      }
      else
	yymsg(0, "No field with id %i", (int)(yyvsp[-8].d));
#endif
      Free((yyvsp[-5].c));
      List_Delete((yyvsp[-2].l));
    }
#line 7158 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 108:
#line 1216 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      Field *field = GModel::current()->getFields()->get((int)(yyvsp[-4].d));
      if(field){
        FieldCallback *callback = field->callbacks[(yyvsp[-1].c)];
        if(callback) {
          callback->run();
        }
        else
          yymsg(0, "Unknown callback '%s' in field %i of type '%s'",
              (yyvsp[-1].c), (int)(yyvsp[-4].d), field->getName());
      }
      else
	yymsg(0, "No field with id %i", (int)(yyvsp[-4].d));
#endif
      Free((yyvsp[-1].c));
    }
#line 7180 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 109:
#line 1237 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_PLUGINS)
      try {
	PluginManager::instance()->setPluginOption((yyvsp[-6].c), (yyvsp[-3].c), (yyvsp[-1].d));
      }
      catch (...) {
	yymsg(0, "Unknown option '%s' or plugin '%s'", (yyvsp[-3].c), (yyvsp[-6].c));
      }
#endif
      Free((yyvsp[-6].c)); Free((yyvsp[-3].c));
    }
#line 7196 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 110:
#line 1249 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_PLUGINS)
      try {
	PluginManager::instance()->setPluginOption((yyvsp[-6].c), (yyvsp[-3].c), (yyvsp[-1].c));
      }
      catch (...) {
	yymsg(0, "Unknown option '%s' or plugin '%s'", (yyvsp[-3].c), (yyvsp[-6].c));
      }
#endif
      Free((yyvsp[-6].c)); Free((yyvsp[-3].c)); Free((yyvsp[-1].c));
    }
#line 7212 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 114:
#line 1267 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[0].c));
      std::vector<double> val(1, 0.);
      if(!gmsh_yysymbols.count(key)){
        gmsh_yysymbols[key].value = val;
      }
      Free((yyvsp[0].c));
    }
#line 7225 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 115:
#line 1276 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-2].c));
      std::vector<double> val(1, (yyvsp[0].d));
      if(!gmsh_yysymbols.count(key)){
        gmsh_yysymbols[key].value = val;
      }
      Free((yyvsp[-2].c));
    }
#line 7238 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 116:
#line 1285 "Gmsh.y" /* yacc.c:1652  */
    { init_options(); }
#line 7244 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 117:
#line 1287 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-3].l)) != 1)
	yymsg(1, "List notation should be used to define list '%s[]'", (yyvsp[-6].c));
      std::string key((yyvsp[-6].c));
      std::vector<double> val;
      for(int i = 0; i < List_Nbr((yyvsp[-3].l)); i++){
        double d;
        List_Read((yyvsp[-3].l), i, &d);
        val.push_back(d);
      }
      if(!gmsh_yysymbols.count(key)){
        Msg::ExchangeOnelabParameter(key, val, floatOptions, charOptions);
        gmsh_yysymbols[key].value = val;
      }
      Free((yyvsp[-6].c));
      List_Delete((yyvsp[-3].l));
    }
#line 7266 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 118:
#line 1305 "Gmsh.y" /* yacc.c:1652  */
    { init_options(); }
#line 7272 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 119:
#line 1307 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-8].c));
      std::vector<double> val;
      for(int i = 0; i < List_Nbr((yyvsp[-3].l)); i++){
        double d;
        List_Read((yyvsp[-3].l), i, &d);
        val.push_back(d);
      }
      if(!gmsh_yysymbols.count(key)){
        Msg::ExchangeOnelabParameter(key, val, floatOptions, charOptions);
        gmsh_yysymbols[key].value = val;
      }
      Free((yyvsp[-8].c));
      List_Delete((yyvsp[-3].l));
    }
#line 7292 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 120:
#line 1323 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-2].c)), val((yyvsp[0].c));
      if(!gmsh_yystringsymbols.count(key)){
        gmsh_yystringsymbols[key] = std::vector<std::string>(1, val);
      }
      Free((yyvsp[-2].c));
      Free((yyvsp[0].c));
    }
#line 7305 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 121:
#line 1332 "Gmsh.y" /* yacc.c:1652  */
    { init_options(); }
#line 7311 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 122:
#line 1334 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-6].c)), val((yyvsp[-3].c));
      if(!gmsh_yystringsymbols.count(key)){
        Msg::ExchangeOnelabParameter(key, val, floatOptions, charOptions);
        gmsh_yystringsymbols[key] = std::vector<std::string>(1, val);
      }
      Free((yyvsp[-6].c));
      Free((yyvsp[-3].c));
    }
#line 7325 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 124:
#line 1348 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string name((yyvsp[0].c));
      Msg::UndefineOnelabParameter(name);
      Free((yyvsp[0].c));
    }
#line 7335 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 125:
#line 1356 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20,20,sizeof(doubleXstring));
      doubleXstring v = {(yyvsp[-2].d), (yyvsp[0].c)};
      List_Add((yyval.l), &v);
    }
#line 7345 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 126:
#line 1362 "Gmsh.y" /* yacc.c:1652  */
    {
      doubleXstring v = {(yyvsp[-2].d), (yyvsp[0].c)};
      List_Add((yyval.l), &v);
    }
#line 7354 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 127:
#line 1367 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20,20,sizeof(doubleXstring));
      int n = List_Nbr((yyvsp[-4].l));
      if(!gmsh_yystringsymbols.count((yyvsp[-2].c))){
	yymsg(0, "Unknown string variable '%s'", (yyvsp[-2].c));
      }
      else{
        std::vector<std::string> &s(gmsh_yystringsymbols[(yyvsp[-2].c)]);
        int m = s.size();
        if(n == m){
          for(int i = 0; i < n; i++){
            double d;
            List_Read((yyvsp[-4].l), i, &d);
            doubleXstring v = {d, strsave((char*)s[i].c_str())};
            List_Add((yyval.l), &v);
          }
        }
        else{
          yymsg(0, "Size mismatch in enumeration: %d != %d", n, m);
        }
      }
      List_Delete((yyvsp[-4].l));
      Free((yyvsp[-2].c));
    }
#line 7383 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 134:
#line 1410 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
        double v;
        List_Read((yyvsp[0].l), i, &v);
        floatOptions[key].push_back(v);
        if (flag_Enum && !i) { member_ValMax = (int)v; }
      }
      Free((yyvsp[-1].c));
      List_Delete((yyvsp[0].l));
    }
#line 7399 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 135:
#line 1422 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[0].c));
      double v;
      if (!flag_Enum) {
        v = 1.;
        if (key == "Enum") flag_Enum = 1;
      }
      else
        v = (double)++member_ValMax;
      floatOptions[key].push_back(v);
      Free((yyvsp[0].c));
    }
#line 7416 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 136:
#line 1435 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-3].c));
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
        doubleXstring v;
        List_Read((yyvsp[-1].l), i, &v);
        floatOptions[key].push_back(v.d);
        charOptions[key].push_back(v.s);
      }
      Free((yyvsp[-3].c));
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++)
        Free(((doubleXstring*)List_Pointer((yyvsp[-1].l), i))->s);
      List_Delete((yyvsp[-1].l));
    }
#line 7434 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 137:
#line 1450 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      std::string val((yyvsp[0].c));
      charOptions[key].push_back(val);
      Free((yyvsp[-1].c));
      Free((yyvsp[0].c));
    }
#line 7446 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 138:
#line 1459 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
        char *v;
        List_Read((yyvsp[0].l), i, &v);
        charOptions[key].push_back(v);
      }
      Free((yyvsp[-1].c));
      List_Delete((yyvsp[0].l));
    }
#line 7461 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 143:
#line 1484 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      double val = (yyvsp[0].d);
      floatOptions[key].push_back(val);
      Free((yyvsp[-1].c));
    }
#line 7472 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 144:
#line 1492 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      std::string val((yyvsp[0].c));
      charOptions[key].push_back(val);
      Free((yyvsp[-1].c));
      Free((yyvsp[0].c));
    }
#line 7484 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 145:
#line 1501 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key("Macro");
      std::string val((yyvsp[0].c));
      charOptions[key].push_back(val);
      Free((yyvsp[0].c));
    }
#line 7495 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 146:
#line 1509 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
        char *s;
        List_Read((yyvsp[0].l), i, &s);
        std::string val(s);
        Free(s);
        charOptions[key].push_back(val);
      }
      Free((yyvsp[-1].c));
      List_Delete((yyvsp[0].l));
    }
#line 7512 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 147:
#line 1523 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string key((yyvsp[-1].c));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
        char *s;
        List_Read((yyvsp[0].l), i, &s);
        std::string val(s);
        Free(s);
        charOptions[key].push_back(val);
      }
      Free((yyvsp[-1].c));
      List_Delete((yyvsp[0].l));
    }
#line 7529 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 148:
#line 1541 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = (int)(yyvsp[0].d);
    }
#line 7537 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 149:
#line 1545 "Gmsh.y" /* yacc.c:1652  */
    {
      int t = GModel::current()->getGEOInternals()->getMaxPhysicalTag();
      GModel::current()->getGEOInternals()->setMaxPhysicalTag(t + 1);
      (yyval.i) = GModel::current()->setPhysicalName(std::string((yyvsp[0].c)), dim_entity, t + 1);
      Free((yyvsp[0].c));
    }
#line 7548 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 150:
#line 1552 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = GModel::current()->setPhysicalName(std::string((yyvsp[-2].c)), dim_entity, (yyvsp[0].d));
      Free((yyvsp[-2].c));
    }
#line 7557 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 151:
#line 1560 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = 0;
    }
#line 7565 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 152:
#line 1564 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(1, 1, sizeof(double));
      double p = (yyvsp[-1].d);
      List_Add((yyval.l), &p);
    }
#line 7575 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 153:
#line 1570 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
    }
#line 7583 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 154:
#line 1574 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      double flag = -1;
      List_Add((yyval.l), &flag);
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++)
        List_Add((yyval.l), List_Pointer((yyvsp[-1].l), i));
      List_Delete((yyvsp[-1].l));
    }
#line 7596 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 155:
#line 1585 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < 4; i++) (yyval.v)[i] = 0.;
    }
#line 7604 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 156:
#line 1589 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < 4; i++) (yyval.v)[i] = (yyvsp[0].v)[i];
    }
#line 7612 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 157:
#line 1596 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      double x = CTX::instance()->geom.scalingFactor * (yyvsp[-1].v)[0];
      double y = CTX::instance()->geom.scalingFactor * (yyvsp[-1].v)[1];
      double z = CTX::instance()->geom.scalingFactor * (yyvsp[-1].v)[2];
      double lc = CTX::instance()->geom.scalingFactor * (yyvsp[-1].v)[3];
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addVertex(num, x, y, z, lc);
      }
      else{
        if(!myGmshSurface)
          r = GModel::current()->getGEOInternals()->addVertex(num, x, y, z, lc);
        else
          r = GModel::current()->getGEOInternals()->addVertex(num, x, y,
                                                              myGmshSurface, lc);
      }
      if(!r) yymsg(0, "Could not add point");
      AddToTemporaryBoundingBox(x, y, z);
      (yyval.s).Type = MSH_POINT;
      (yyval.s).Num = num;
    }
#line 7639 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 158:
#line 1619 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addLine(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addLine(num, tags);
      }
      if(!r) yymsg(0, "Could not add line");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SEGM_LINE;
      (yyval.s).Num = num;
    }
#line 7659 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 159:
#line 1635 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addSpline(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addSpline(num, tags);
      }
      if(!r) yymsg(0, "Could not add spline");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SEGM_SPLN;
      (yyval.s).Num = num;
    }
#line 7679 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 160:
#line 1651 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-5].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-2].l), tags);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-2].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(tags.size() == 3){
          r = GModel::current()->getOCCInternals()->addCircleArc
            (num, tags[0], tags[1], tags[2]);
        }
        else if(param.size() >= 4 && param.size() <= 6){
          double r = param[3];
          double a1 = (param.size() == 6) ? param[4] : 0.;
          double a2 = (param.size() == 6) ? param[5] :
            (param.size() == 5) ? param[4] : 2.*M_PI;
          r = GModel::current()->getOCCInternals()->addCircle
            (num, param[0], param[1], param[2], r, a1, a2);
        }
        else{
          yymsg(0, "Circle requires 3 points or 4 to 6 parameters");
        }
      }
      else{
        if(tags.size() == 3){
          r = GModel::current()->getGEOInternals()->addCircleArc
            (num, tags[0], tags[1], tags[2], (yyvsp[-1].v)[0], (yyvsp[-1].v)[1], (yyvsp[-1].v)[2]);
        }
        else{
          yymsg(0, "Circle requires 3 points");
        }
      }
      if(!r) yymsg(0, "Could not add circle");
      List_Delete((yyvsp[-2].l));
      (yyval.s).Type = MSH_SEGM_CIRC;
      (yyval.s).Num = num;
    }
#line 7720 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 161:
#line 1688 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-5].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-2].l), tags);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-2].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(tags.size() == 3){
          r = GModel::current()->getOCCInternals()->addEllipseArc
            (num, tags[0], tags[1], tags[2]);
        }
        else if(tags.size() == 4){
          r = GModel::current()->getOCCInternals()->addEllipseArc
            (num, tags[0], tags[1], tags[3]);
        }
        else if(param.size() >= 5 && param.size() <= 7){
          double a1 = (param.size() == 7) ? param[5] : 0.;
          double a2 = (param.size() == 7) ? param[6] :
            (param.size() == 6) ? param[5] : 2.*M_PI;
          r = GModel::current()->getOCCInternals()->addEllipse
            (num, param[0], param[1], param[2], param[3], param[4], a1, a2);
        }
        else{
          yymsg(0, "Ellipse requires 3 or 4 points, or 5 to 7 parameters");
        }
      }
      else{
        if(tags.size() == 4){
          r = GModel::current()->getGEOInternals()->addEllipseArc
            (num, tags[0], tags[1], tags[2], tags[3], (yyvsp[-1].v)[0], (yyvsp[-1].v)[1], (yyvsp[-1].v)[2]);
        }
        else{
          yymsg(0, "Ellipse requires 4 points");
        }
      }
      if(!r) yymsg(0, "Could not add ellipse");
      List_Delete((yyvsp[-2].l));
      (yyval.s).Type = MSH_SEGM_ELLI;
      (yyval.s).Num = num;
    }
#line 7764 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 162:
#line 1728 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addBSpline(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addBSpline(num, tags);
      }
      if(!r) yymsg(0, "Could not add BSpline");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SEGM_BSPLN;
      (yyval.s).Num = num;
    }
#line 7784 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 163:
#line 1744 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addBezier(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addBezier(num, tags);
      }
      if(!r) yymsg(0, "Could not add Bezier");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SEGM_BEZIER;
      (yyval.s).Num = num;
    }
#line 7804 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 164:
#line 1761 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-8].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-5].l), tags);
      std::vector<double> seqknots; ListOfDouble2Vector((yyvsp[-3].l), seqknots);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        int degree = (int)(yyvsp[-1].d);
        std::vector<double> weights, knots;
        std::vector<int> mults;
        for(std::size_t i = 0; i < seqknots.size(); i++){
          if(!i || (i && fabs(seqknots[i] - seqknots[i - 1]) > 1e-12)){
            knots.push_back(seqknots[i]);
            mults.push_back(1);
          }
          else{
            mults.back() += 1;
          }
        }
        r = GModel::current()->getOCCInternals()->addBSpline
          (num, tags, degree, weights, knots, mults);
      }
      else{
        r = GModel::current()->getGEOInternals()->addBSpline(num, tags, seqknots);
      }
      if(!r) yymsg(0, "Could not add nurbs");
      List_Delete((yyvsp[-5].l));
      List_Delete((yyvsp[-3].l));
      (yyval.s).Type = MSH_SEGM_NURBS;
      (yyval.s).Num = num;
    }
#line 7839 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 165:
#line 1792 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addWire(num, tags, false);
      }
      else{
        yymsg(0, "Wire only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add wire");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SEGM_LOOP;
      (yyval.s).Num = num;
    }
#line 7859 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 166:
#line 1808 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addLineLoop(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addLineLoop(num, tags);
      }
      if(!r) yymsg(0, "Could not add line loop");
      List_Delete((yyvsp[-1].l));
      Free((yyvsp[-6].c));
      (yyval.s).Type = MSH_SEGM_LOOP;
      (yyval.s).Num = num;
    }
#line 7880 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 167:
#line 1825 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addPlaneSurface(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addPlaneSurface(num, tags);
      }
      if(!r) yymsg(0, "Could not add plane surface");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SURF_PLAN;
      (yyval.s).Num = num;
    }
#line 7900 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 168:
#line 1841 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-5].d);
      std::vector<int> wires; ListOfDouble2Vector((yyvsp[-2].l), wires);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(wires.size() != 1){
          yymsg(0, "OpenCASCADE surface filling requires a single line loop");
        }
        else{
          std::vector<int> constraints; ListOfDouble2Vector((yyvsp[-1].l), constraints);
          std::vector<int> points, surfaces, continuity;
          if(constraints.size() >= 3 && constraints[0] < 0){
            // {-1, type, ent, type, ent, ...}
            for(std::size_t i = 2; i < constraints.size(); i+=2){
              int type = constraints[i - 1];
              if(type == 0){
                points.push_back(constraints[i]);
              }
              else if(type == 1 || type == 2){
                surfaces.push_back(constraints[i]);
                continuity.push_back(type);
              }
              else
                yymsg(0, "Unknown type of constraint for surface filling");
            }
          }
          else if(constraints.size() > 0){
            // {point, point, ...}
            points = constraints;
          }
          r = GModel::current()->getOCCInternals()->addSurfaceFilling
            (num, wires[0], points, surfaces, continuity);
        }
      }
      else{
        int sphereCenter = -1;
        if(List_Nbr((yyvsp[-1].l)) == 1){
          double d; List_Read((yyvsp[-1].l), 0, &d);
          sphereCenter = (int)d;
        }
        r = GModel::current()->getGEOInternals()->addSurfaceFilling
          (num, wires, sphereCenter);
      }
      if(!r) yymsg(0, "Could not add surface");
      List_Delete((yyvsp[-2].l));
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SURF_REGL;
      (yyval.s).Num = num;
    }
#line 7954 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 169:
#line 1891 "Gmsh.y" /* yacc.c:1652  */
    {
      yymsg(2, "'Ruled Surface' command is deprecated: use 'Surface' instead");
      int num = (int)(yyvsp[-5].d);
      std::vector<int> wires; ListOfDouble2Vector((yyvsp[-2].l), wires);
      int sphereCenter = -1;
      if(List_Nbr((yyvsp[-1].l)) == 1){
        double d; List_Read((yyvsp[-1].l), 0, &d);
        sphereCenter = (int)d;
      }
      bool r = GModel::current()->getGEOInternals()->addSurfaceFilling
        (num, wires, sphereCenter);
      if(!r) yymsg(0, "Could not add surface");
      List_Delete((yyvsp[-2].l));
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type =  MSH_SURF_REGL;
      (yyval.s).Num = num;
    }
#line 7976 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 170:
#line 1909 "Gmsh.y" /* yacc.c:1652  */
    {
      myGmshSurface = 0;
      (yyval.s).Type = 0;
      (yyval.s).Num = 0;
    }
#line 7986 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 171:
#line 1915 "Gmsh.y" /* yacc.c:1652  */
    {
      myGmshSurface = gmshSurface::getSurface((int)(yyvsp[-1].d));
      (yyval.s).Type = 0;
      (yyval.s).Num = 0;
    }
#line 7996 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 172:
#line 1921 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-6].d);
      myGmshSurface = gmshParametricSurface::NewParametricSurface(num, (yyvsp[-3].c), (yyvsp[-2].c), (yyvsp[-1].c));
      (yyval.s).Type = 0;
      (yyval.s).Num = num;
    }
#line 8007 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 173:
#line 1928 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      (yyval.s).Type = 0;
      bool r = true;
      if(param.size() >= 4 && param.size() <= 7){
        if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
          double a1 = (param.size() >= 5) ? param[4] : -M_PI/2;
          double a2 = (param.size() >= 6) ? param[5] : M_PI/2;
          double a3 = (param.size() >= 7) ? param[6] : 2.*M_PI;
          r = GModel::current()->getOCCInternals()->addSphere
            (num, param[0], param[1], param[2], param[3], a1, a2, a3);
        }
        else{
          yymsg(0, "Sphere only available with OpenCASCADE geometry kernel");
        }
        (yyval.s).Type = MSH_VOLUME;
      }
      else if(tags.size() == 2){
        myGmshSurface = GModel::current()->getGEOInternals()->newGeometrySphere
          (num, tags[0], tags[1]);
      }
      else{
        yymsg(0, "Sphere requires 2 points or from 4 to 7 parameters");
      }
      if(!r) yymsg(0, "Could not add sphere");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Num = num;
    }
#line 8042 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 174:
#line 1959 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      if(tags.size() == 2){
        myGmshSurface = GModel::current()->getGEOInternals()->newGeometryPolarSphere
          (num, tags[0], tags[1]);
      }
      else{
        yymsg(0, "PolarSphere requires 2 points");
      }
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = 0;
      (yyval.s).Num = num;
    }
#line 8061 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 175:
#line 1974 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 6){
          r = GModel::current()->getOCCInternals()->addBox
            (num, param[0], param[1], param[2], param[3], param[4], param[5]);
        }
        else{
          yymsg(0, "Box requires 6 parameters");
        }
      }
      else{
        yymsg(0, "Box only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add block");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8087 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 176:
#line 1996 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 5 || param.size() == 6){
          double alpha = (param.size() == 6) ? param[5] : 2*M_PI;
          r = GModel::current()->getOCCInternals()->addTorus
            (num, param[0], param[1], param[2], param[3], param[4], alpha);
        }
        else{
          yymsg(0, "Torus requires 5 ou 6 parameters");
        }
      }
      else{
        yymsg(0, "Torus only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add torus");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8114 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 177:
#line 2019 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 5 || param.size() == 6){
          double r = (param.size() == 6) ? param[5] : 0.;
          r = GModel::current()->getOCCInternals()->addRectangle
            (num, param[0], param[1], param[2], param[3], param[4], r);
        }
        else{
          yymsg(0, "Rectangle requires 5 ou 6 parameters");
        }
      }
      else{
        yymsg(0, "Rectangle only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add rectangle");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SURF_PLAN;
      (yyval.s).Num = num;
    }
#line 8141 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 178:
#line 2042 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 4 || param.size() == 5){
          double ry = (param.size() == 5) ? param[4] : param[3];
          r = GModel::current()->getOCCInternals()->addDisk
            (num, param[0], param[1], param[2], param[3], ry);
        }
        else{
          yymsg(0, "Disk requires 4 or 5 parameters");
        }
      }
      else{
        yymsg(0, "Disk only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add disk");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_SURF_PLAN;
      (yyval.s).Num = num;
    }
#line 8168 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 179:
#line 2065 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 7 || param.size() == 8){
          double angle = (param.size() == 8) ? param[7] : 2*M_PI;
          r = GModel::current()->getOCCInternals()->addCylinder
            (num, param[0], param[1], param[2], param[3], param[4], param[5],
             param[6], angle);
        }
        else{
          yymsg(0, "Cylinder requires 7 or 8 parameters");
        }
      }
      else{
        yymsg(0, "Cylinder only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add cylinder");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8196 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 180:
#line 2089 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 8 || param.size() == 9){
          double alpha = (param.size() == 9) ? param[8] : 2*M_PI;
          r = GModel::current()->getOCCInternals()->addCone
            (num, param[0], param[1], param[2], param[3], param[4], param[5],
             param[6], param[7], alpha);
        }
        else{
          yymsg(0, "Cone requires 8 or 9 parameters");
        }
      }
      else{
        yymsg(0, "Cone only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add cone");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8224 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 181:
#line 2113 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() == 6 || param.size() == 7){
          double ltx = (param.size() == 7) ? param[6] : 0.;
          r = GModel::current()->getOCCInternals()->addWedge
            (num, param[0], param[1], param[2], param[3], param[4], param[5],
             ltx);
        }
        else{
          yymsg(0, "Wedge requires 7 parameters");
        }
      }
      else{
        yymsg(0, "Wedge only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add wedge");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8252 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 182:
#line 2137 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<double> param; ListOfDouble2Vector((yyvsp[-1].l), param);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        if(param.size() >= 2){
          int in = (int)param[0];
          double offset = param[1];
          std::vector<int> exclude;
          for(std::size_t i = 2; i < param.size(); i++)
            exclude.push_back(param[i]);
          std::vector<std::pair<int, int> > outDimTags;
          r = GModel::current()->getOCCInternals()->addThickSolid
            (num, in, exclude, offset, outDimTags);
        }
        else{
          yymsg(0, "ThickSolid requires at least 2 parameters");
        }
      }
      else{
        yymsg(0, "ThickSolid only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add thick solid");
      List_Delete((yyvsp[-1].l));
    }
#line 8282 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 183:
#line 2163 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addSurfaceLoop(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addSurfaceLoop(num, tags);
      }
      if(!r) yymsg(0, "Could not add surface loop");
      List_Delete((yyvsp[-1].l));
      Free((yyvsp[-6].c));
      (yyval.s).Type = MSH_SURF_LOOP;
      (yyval.s).Num = num;
    }
#line 8303 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 184:
#line 2180 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addVolume(num, tags);
      }
      else{
        r = GModel::current()->getGEOInternals()->addVolume(num, tags);
      }
      if(!r) yymsg(0, "Could not add volume");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8323 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 185:
#line 2196 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> wires; ListOfDouble2Vector((yyvsp[-1].l), wires);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<std::pair<int, int> > outDimTags;
        r = GModel::current()->getOCCInternals()->addThruSections
          (num, wires, true, false, outDimTags);
      }
      else{
        yymsg(0, "ThruSections only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add thrusections");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8345 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 186:
#line 2214 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].d);
      std::vector<int> wires; ListOfDouble2Vector((yyvsp[-1].l), wires);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<std::pair<int, int> > outDimTags;
        r = GModel::current()->getOCCInternals()->addThruSections
          (num, wires, true, true, outDimTags);
      }
      else{
        yymsg(0, "ThruSections only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add ruled thrusections");
      List_Delete((yyvsp[-1].l));
      (yyval.s).Type = MSH_VOLUME;
      (yyval.s).Num = num;
    }
#line 8367 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 187:
#line 2232 "Gmsh.y" /* yacc.c:1652  */
    {
      yymsg(0, "Compounds entities are deprecated: use Compound meshing constraints "
            "instead, i.e. Compound %s { ... };", ((yyvsp[-6].i) == 2) ? "Surface" : "Curve");
      (yyval.s).Type = 0;
      (yyval.s).Num = 0;
    }
#line 8378 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 188:
#line 2240 "Gmsh.y" /* yacc.c:1652  */
    {
      yymsg(0, "Compounds entities are deprecated: use Compound meshing constraints "
            "instead, i.e. Compound %s { ... };", ((yyvsp[-10].i) == 2) ? "Surface" : "Curve");
      (yyval.s).Type = 0;
      (yyval.s).Num = 0;
    }
#line 8389 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 189:
#line 2247 "Gmsh.y" /* yacc.c:1652  */
    {
      dim_entity = (yyvsp[0].i);
    }
#line 8397 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 190:
#line 2251 "Gmsh.y" /* yacc.c:1652  */
    {
      int num = (int)(yyvsp[-4].i);
      int op = (yyvsp[-2].i);
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      bool r = GModel::current()->getGEOInternals()->modifyPhysicalGroup
        ((yyvsp[-7].i), num, op, tags);
      if(!r)
        switch ((yyvsp[-7].i)) {
        case 0: yymsg(0, "Could not modify physical point"); break;
        case 1: yymsg(0, "Could not modify physical line"); break;
        case 2: yymsg(0, "Could not modify physical surface"); break;
        case 3: yymsg(0, "Could not modify physical volume"); break;
        }
      List_Delete((yyvsp[-1].l));
      switch ((yyvsp[-7].i)) {
      case 0: (yyval.s).Type = MSH_PHYSICAL_POINT  ; break;
      case 1: (yyval.s).Type = MSH_PHYSICAL_LINE   ; break;
      case 2: (yyval.s).Type = MSH_PHYSICAL_SURFACE; break;
      case 3: (yyval.s).Type = MSH_PHYSICAL_VOLUME ; break;
      }
      (yyval.s).Num = num;
    }
#line 8424 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 191:
#line 2277 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 0; }
#line 8430 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 192:
#line 2279 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 8436 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 193:
#line 2281 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 2; }
#line 8442 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 194:
#line 2283 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 3; }
#line 8448 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 195:
#line 2285 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = (int)(yyvsp[-1].d);
      if ((yyval.i)<0 || (yyval.i)>3) yymsg(0, "GeoEntity dim out of range [0,3]");
    }
#line 8457 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 196:
#line 2293 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 8463 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 197:
#line 2295 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 2; }
#line 8469 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 198:
#line 2297 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 3; }
#line 8475 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 199:
#line 2299 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = (int)(yyvsp[-1].d);
      if ((yyval.i)<1 || (yyval.i)>3) yymsg(0, "GeoEntity dim out of range [1,3]");
    }
#line 8484 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 200:
#line 2307 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 8490 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 201:
#line 2309 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 2; }
#line 8496 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 202:
#line 2311 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = (int)(yyvsp[-1].d);
      if ((yyval.i)<1 || (yyval.i)>2) yymsg(0, "GeoEntity dim out of range [1,2]");
    }
#line 8505 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 203:
#line 2319 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 0; }
#line 8511 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 204:
#line 2321 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 8517 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 205:
#line 2323 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 2; }
#line 8523 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 206:
#line 2325 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = (int)(yyvsp[-1].d);
      if ((yyval.i)<0 || (yyval.i)>2) yymsg(0, "GeoEntity dim out of range [0,2]");
    }
#line 8532 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 207:
#line 2335 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->translate
          (dimTags, (yyvsp[-3].v)[0], (yyvsp[-3].v)[1], (yyvsp[-3].v)[2]);
      }
      else{
        r = GModel::current()->getGEOInternals()->translate
          (dimTags, (yyvsp[-3].v)[0], (yyvsp[-3].v)[1], (yyvsp[-3].v)[2]);
      }
      if(!r) yymsg(0, "Could not translate shapes");
      (yyval.l) = (yyvsp[-1].l);
    }
#line 8552 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 208:
#line 2351 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->rotate
          (dimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-4].d));
      }
      else{
        r = GModel::current()->getGEOInternals()->rotate
          (dimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-4].d));
      }
      if(!r) yymsg(0, "Could not rotate shapes");
      (yyval.l) = (yyvsp[-1].l);
    }
#line 8572 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 209:
#line 2367 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->symmetry
          (dimTags, (yyvsp[-3].v)[0], (yyvsp[-3].v)[1], (yyvsp[-3].v)[2], (yyvsp[-3].v)[3]);
      }
      else{
        r = GModel::current()->getGEOInternals()->symmetry
          (dimTags, (yyvsp[-3].v)[0], (yyvsp[-3].v)[1], (yyvsp[-3].v)[2], (yyvsp[-3].v)[3]);
      }
      if(!r) yymsg(0, "Could not apply symmetry transform");
      (yyval.l) = (yyvsp[-1].l);
    }
#line 8592 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 210:
#line 2383 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->dilate
          (dimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-4].d), (yyvsp[-4].d), (yyvsp[-4].d));
      }
      else{
        r = GModel::current()->getGEOInternals()->dilate
          (dimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-4].d), (yyvsp[-4].d), (yyvsp[-4].d));
      }
      if(!r) yymsg(0, "Could not dilate shapes");
      (yyval.l) = (yyvsp[-1].l);
    }
#line 8612 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 211:
#line 2399 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->dilate
          (dimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]);
      }
      else{
        r = GModel::current()->getGEOInternals()->dilate
          (dimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]);
      }
      if(!r) yymsg(0, "Could not dilate shapes");
      (yyval.l) = (yyvsp[-1].l);
    }
#line 8632 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 212:
#line 2415 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<double> mat;
        ListOfDouble2Vector((yyvsp[-4].l), mat);
        r = GModel::current()->getOCCInternals()->affine(dimTags, mat);
      }
      else{
        yymsg(0, "Affine transform only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not transform shapes");
      List_Delete((yyvsp[-4].l));
      (yyval.l) = (yyvsp[-1].l);
    }
#line 8653 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 213:
#line 2432 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), inDimTags);
      (yyval.l) = (yyvsp[-1].l);
      List_Reset((yyval.l));
      std::string action((yyvsp[-3].c));
      bool r = true;
      if(action == "Duplicata"){
        if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
          r = GModel::current()->getOCCInternals()->copy(inDimTags, outDimTags);
        }
        else{
          r = GModel::current()->getGEOInternals()->copy(inDimTags, outDimTags);
        }
      }
      else if(action == "Boundary" || action == "CombinedBoundary" ||
              action == "PointsOf"){
        // boundary operations are performed directly on GModel, which enables
        // to compute the boundary of hybrid CAD models; this also automatically
        // binds all boundary entities for OCC models
        if(GModel::current()->getOCCInternals() &&
           GModel::current()->getOCCInternals()->getChanged())
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
        if(GModel::current()->getGEOInternals()->getChanged())
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
        r = GModel::current()->getBoundaryTags
          (inDimTags, outDimTags, action == "CombinedBoundary", true,
           action == "PointsOf");
      }
      else{
        yymsg(0, "Unknown action on multiple shapes '%s'", (yyvsp[-3].c));
      }
      if(!r) yymsg(0, "Could not apply operation on shapes");
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
      Free((yyvsp[-3].c));
    }
#line 8694 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 214:
#line 2469 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        yymsg(0, "Intersect line not available with OpenCASCADE geometry kernel");
      }
      else{
        std::vector<int> in, out; ListOfDouble2Vector((yyvsp[-5].l), in);
        r = GModel::current()->getGEOInternals()->intersectCurvesWithSurface
          (in, (int)(yyvsp[-1].d), out);
        for(std::size_t i = 0; i < out.size(); i++){
          Shape s;
          s.Type = MSH_POINT;
          s.Num = out[i];
          List_Add((yyval.l), &s);
        }
      }
      if(!r) yymsg(0, "Could not intersect line");
      List_Delete((yyvsp[-5].l));
    }
#line 8719 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 215:
#line 2491 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        yymsg(0, "Split Line not available with OpenCASCADE geometry kernel");
      }
      else{
        std::vector<int> vertices, curves; ListOfDouble2Vector((yyvsp[-2].l), vertices);
        r = GModel::current()->getGEOInternals()->splitCurve
          ((int)(yyvsp[-5].d), vertices, curves);
        for(std::size_t i = 0; i < curves.size(); i++){
          Shape s;
          s.Type = MSH_SEGM_LINE;
          s.Num = curves[i];
          List_Add((yyval.l), &s);
        }
      }
      if(!r) yymsg(0, "Could not split line");
      List_Delete((yyvsp[-2].l));
    }
#line 8744 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 216:
#line 2514 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[0].l); }
#line 8750 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 217:
#line 2515 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[0].l); }
#line 8756 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 218:
#line 2520 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(3, 3, sizeof(Shape));
    }
#line 8764 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 219:
#line 2524 "Gmsh.y" /* yacc.c:1652  */
    {
      List_Add((yyval.l), &(yyvsp[0].s));
    }
#line 8772 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 220:
#line 2528 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
	double d;
	List_Read((yyvsp[-2].l), i, &d);
	Shape s;
	s.Num = (int)d;
        switch ((yyvsp[-4].i)) {
        case 0: s.Type = MSH_POINT    ; break;
        case 1: s.Type = MSH_SEGM_LINE; break;
        case 2: s.Type = MSH_SURF_PLAN; break; // we don't care about the actual type
        case 3: s.Type = MSH_VOLUME   ; break;
        }
        List_Add((yyval.l), &s);
      }
      List_Delete((yyvsp[-2].l));
    }
#line 8793 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 221:
#line 2545 "Gmsh.y" /* yacc.c:1652  */
    {
      List_T *tmp = List_Create(10, 10, sizeof(double));
      getElementaryTagsForPhysicalGroups((yyvsp[-4].i), (yyvsp[-2].l), tmp);
      for(int i = 0; i < List_Nbr(tmp); i++){
	double d;
	List_Read(tmp, i, &d);
 	Shape s;
	s.Num = (int)d; // FIXME
        switch ((yyvsp[-4].i)) {
        case 0: s.Type = MSH_POINT    ; break;
        case 1: s.Type = MSH_SEGM_LINE; break;
        case 2: s.Type = MSH_SURF_PLAN; break; // we don't care about the actual type
        case 3: s.Type = MSH_VOLUME   ; break;
        }
        List_Add((yyval.l), &s);
      }
      List_Delete(tmp);
      List_Delete((yyvsp[-2].l));
    }
#line 8817 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 222:
#line 2565 "Gmsh.y" /* yacc.c:1652  */
    {
      List_T *tmp = List_Create(10, 10, sizeof(double));
      getParentTags((yyvsp[-4].i), (yyvsp[-2].l), tmp);
      for(int i = 0; i < List_Nbr(tmp); i++){
	double d;
	List_Read(tmp, i, &d);
 	Shape s;
	s.Num = (int)d; // FIXME
        switch ((yyvsp[-4].i)) {
        case 0: s.Type = MSH_POINT    ; break;
        case 1: s.Type = MSH_SEGM_LINE; break;
        case 2: s.Type = MSH_SURF_PLAN; break; // we don't care about the actual type
        case 3: s.Type = MSH_VOLUME   ; break;
        }
        List_Add((yyval.l), &s);
      }
      List_Delete(tmp);
      List_Delete((yyvsp[-2].l));
    }
#line 8841 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 223:
#line 2585 "Gmsh.y" /* yacc.c:1652  */
    {
      List_T *tmp = List_Create(10, 10, sizeof(double));
      getAllElementaryTags((yyvsp[-4].i), tmp);
      for(int i = 0; i < List_Nbr(tmp); i++){
	double d;
	List_Read(tmp, i, &d);
	Shape s;
	s.Num = (int)d;
        switch ((yyvsp[-4].i)) {
        case 0: s.Type = MSH_POINT    ; break;
        case 1: s.Type = MSH_SEGM_LINE; break;
        case 2: s.Type = MSH_SURF_PLAN; break; // we don't care about the actual type
        case 3: s.Type = MSH_VOLUME   ; break;
        }
        List_Add((yyval.l), &s);
      }
      List_Delete(tmp);
    }
#line 8864 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 224:
#line 2604 "Gmsh.y" /* yacc.c:1652  */
    {
      List_T *tmp = List_Create(10, 10, sizeof(double));
      List_T *tmp2 = List_Create(10, 10, sizeof(double));
      getAllPhysicalTags((yyvsp[-4].i), tmp2);
      getElementaryTagsForPhysicalGroups((yyvsp[-4].i), tmp2, tmp);
      for(int i = 0; i < List_Nbr(tmp); i++){
	double d;
	List_Read(tmp, i, &d);
 	Shape s;
	s.Num = (int)d; // FIXME
        switch ((yyvsp[-4].i)) {
        case 0: s.Type = MSH_POINT    ; break;
        case 1: s.Type = MSH_SEGM_LINE; break;
        case 2: s.Type = MSH_SURF_PLAN; break; // we don't care about the actual type
        case 3: s.Type = MSH_VOLUME   ; break;
        }
        List_Add((yyval.l), &s);
      }
      List_Delete(tmp);
      List_Delete(tmp2);
    }
#line 8890 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 225:
#line 2631 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-1].l)) == 4){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d[4];
          for(int i = 0; i < 4; i++)
            List_Read((yyvsp[-1].l), i, &d[i]);
          gLevelset *ls = new gLevelsetPlane(d[0], d[1], d[2], d[3], t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      List_Delete((yyvsp[-1].l));
    }
#line 8913 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 226:
#line 2650 "Gmsh.y" /* yacc.c:1652  */
    {
      int t = (int)(yyvsp[-6].d);
      if(gLevelset::find(t)){
	yymsg(0, "Levelset %d already exists", t);
      }
      else {
	fullMatrix<double> centers(List_Nbr((yyvsp[-2].l)),3);
	for (int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
	  List_T *l = *(List_T**)List_Pointer((yyvsp[-2].l), i);
	  for (int j = 0; j < List_Nbr(l); j++){
	    centers(i,j) = (double)(*(double*)List_Pointer(l, j));
	  }
	}
        gLevelset *ls = new gLevelsetPoints(centers, t);
        gLevelset::add(ls);
      }
      for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++)
        List_Delete(*(List_T**)List_Pointer((yyvsp[-2].l), i));
      List_Delete((yyvsp[-2].l));
    }
#line 8938 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 227:
#line 2672 "Gmsh.y" /* yacc.c:1652  */
    {
      int t = (int)(yyvsp[-10].d);
      if(gLevelset::find(t)){
        yymsg(0, "Levelset %d already exists", t);
      }
      else {
        double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
        double n[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
        gLevelset *ls = new gLevelsetPlane(pt, n, t);
        gLevelset::add(ls);
      }
      List_Delete((yyvsp[-2].l));
    }
#line 8956 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 228:
#line 2687 "Gmsh.y" /* yacc.c:1652  */
    {
      int t = (int)(yyvsp[-12].d);
      if(gLevelset::find(t)){
        yymsg(0, "Levelset %d already exists", t);
      }
      else {
        double pt1[3] = {(yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2]};
        double pt2[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
        double pt3[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
        gLevelset *ls = new gLevelsetPlane(pt1, pt2, pt3, t);
        gLevelset::add(ls);
      }
      List_Delete((yyvsp[-2].l));
    }
#line 8975 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 229:
#line 2702 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-2].l)) == 1){
        int t = (int)(yyvsp[-8].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d;
          List_Read((yyvsp[-2].l), 0, &d);
          gLevelset *ls = new gLevelsetSphere((yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2], d, t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      List_Delete((yyvsp[-2].l));
    }
#line 8997 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 230:
#line 2721 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-2].l)) == 1){
        int t = (int)(yyvsp[-10].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d;
          List_Read((yyvsp[-2].l), 0, &d);
          double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
          double dir[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
          gLevelset *ls = new gLevelsetGenCylinder(pt, dir, d, t);
          gLevelset::add(ls);
        }
      }
      else if(List_Nbr((yyvsp[-2].l)) == 2){
        int t = (int)(yyvsp[-10].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d[2];
          for(int i = 0; i < 2; i++)
            List_Read((yyvsp[-2].l), i, &d[i]);
          double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
          double dir[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
          gLevelset *ls = new gLevelsetCylinder(pt, dir, d[0], d[1], t);
          gLevelset::add(ls);
        }
      }
      else if(List_Nbr((yyvsp[-2].l)) == 3){
        int t = (int)(yyvsp[-10].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d[3];
          for(int i = 0; i < 3; i++)
            List_Read((yyvsp[-2].l), i, &d[i]);
          double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
          double dir[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
          gLevelset *ls = new gLevelsetCylinder(pt, dir, d[0], d[1], d[2], t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      List_Delete((yyvsp[-2].l));
    }
#line 9051 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 231:
#line 2772 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-2].l)) == 1){
        int t = (int)(yyvsp[-10].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d;
          List_Read((yyvsp[-2].l), 0, &d);
          double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
          double dir[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
          gLevelset *ls = new gLevelsetCone(pt, dir, d, t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      List_Delete((yyvsp[-2].l));
    }
#line 9075 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 232:
#line 2793 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-2].l)) == 3){
        int t = (int)(yyvsp[-10].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d[3];
          for(int i = 0; i < 3; i++)
            List_Read((yyvsp[-2].l), i, &d[i]);
          double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
          double dir[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
          gLevelset *ls = new gLevelsetEllipsoid(pt, dir, d[0], d[1], d[2], t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      List_Delete((yyvsp[-2].l));
    }
#line 9100 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 233:
#line 2815 "Gmsh.y" /* yacc.c:1652  */
    {
      if(List_Nbr((yyvsp[-2].l)) == 5){
        int t = (int)(yyvsp[-10].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d[5];
          for(int i = 0; i < 5; i++)
            List_Read((yyvsp[-2].l), i, &d[i]);
          double pt[3] = {(yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2]};
          double dir[3] = {(yyvsp[-4].v)[0], (yyvsp[-4].v)[1], (yyvsp[-4].v)[2]};
          gLevelset *ls = new gLevelsetGeneralQuadric(pt, dir, d[0], d[1],
                                                      d[2], d[3], d[4], t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      List_Delete((yyvsp[-2].l));
    }
#line 9126 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 234:
#line 2837 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-6].c), "Union")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          std::vector<gLevelset *> vl;
          for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++) {
            double d; List_Read((yyvsp[-1].l), i, &d);
            gLevelset *pl = gLevelset::find((int)d);
	    if(!pl) yymsg(0, "Unknown levelset %d", (int)d);
            else vl.push_back(pl);
          }
          gLevelset *ls = new gLevelsetUnion(vl, true, t);
          gLevelset::add(ls);
        }
      }
      else if(!strcmp((yyvsp[-6].c), "Intersection")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          std::vector<gLevelset *> vl;
          for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++) {
            double d; List_Read((yyvsp[-1].l), i, &d);
            gLevelset *pl = gLevelset::find((int)d);
	    if(!pl) yymsg(0, "Unknown levelset %d", (int)d);
            else vl.push_back(pl);
          }
          gLevelset *ls = new gLevelsetIntersection(vl, true, t);
          gLevelset::add(ls);
        }
      }
      else if(!strcmp((yyvsp[-6].c), "Cut")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          std::vector<gLevelset *> vl;
          for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++) {
            double d; List_Read((yyvsp[-1].l), i, &d);
            gLevelset *pl = gLevelset::find((int)d);
	    if(!pl) yymsg(0, "Unknown levelset %d", (int)d);
            else vl.push_back(pl);
          }
          gLevelset *ls = new gLevelsetCut(vl, true, t);
          gLevelset::add(ls);
        }
      }
      else if(!strcmp((yyvsp[-6].c), "Crack")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          std::vector<gLevelset *> vl;
          for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++) {
            double d; List_Read((yyvsp[-1].l), i, &d);
            gLevelset *pl = gLevelset::find((int)d);
	    if(!pl) yymsg(0, "Unknown levelset %d", (int)d);
            else vl.push_back(pl);
          }
          gLevelset *ls = new gLevelsetCrack(vl, false, t);
          gLevelset::add(ls);
        }
      }
      else if(!strcmp((yyvsp[-6].c), "Reverse")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          double d;
          List_Read((yyvsp[-1].l), 0, &d);
          gLevelset *pl = gLevelset::find((int)d);
          gLevelset *ls = nullptr;
          if(!pl) yymsg(0, "Unknown levelset %d", (int)d);
          else ls = new gLevelsetReverse(pl, t);
          if(ls) gLevelset::add(ls);
        }
      }
#if defined(HAVE_POST)
      else if(!strcmp((yyvsp[-6].c), "PostView")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          if(List_Nbr((yyvsp[-1].l)) > 0){
            double d; List_Read((yyvsp[-1].l), 0, &d);
            gLevelset *ls = new gLevelsetPostView((int)d, t);
            gLevelset::add(ls);
          }
        }
      }
#endif
      else
        yymsg(0, "Wrong number of arguments for levelset definition");
      Free((yyvsp[-6].c));
      List_Delete((yyvsp[-1].l));
    }
#line 9235 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 235:
#line 2942 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-6].c), "MathEval")){
        int t = (int)(yyvsp[-4].d);
        if(gLevelset::find(t)){
	  yymsg(0, "Levelset %d already exists", t);
        }
        else {
          gLevelset *ls = new gLevelsetMathEval((yyvsp[-1].c), t);
          gLevelset::add(ls);
        }
      }
      else
        yymsg(0, "Unknown levelset '%s'", (yyvsp[-6].c));
      Free((yyvsp[-6].c)); Free((yyvsp[-1].c));
    }
#line 9255 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 236:
#line 2958 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-4].c), "CutMesh")){
        int t = (int)(yyvsp[-2].d);
        if(gLevelset::find(t)){
          GModel::current()->buildCutGModel(gLevelset::find(t), true, false);
        }
        else
          yymsg(0, "Unknown levelset %d", t);
      }
      else if(!strcmp((yyvsp[-4].c), "CutMeshTri")){
        int t = (int)(yyvsp[-2].d);
        if(gLevelset::find(t)){
          GModel::current()->buildCutGModel(gLevelset::find(t), true, true);
        }
        else
          yymsg(0, "Unknown levelset %d", t);
      }
      else if(!strcmp((yyvsp[-4].c), "SplitMesh")){
        int t = (int)(yyvsp[-2].d);
        if(gLevelset::find(t)){
          GModel::current()->buildCutGModel(gLevelset::find(t), false, true);
        }
        else
          yymsg(0, "Unknown levelset %d", t);
      }
      else
        yymsg(0, "Unknown levelset '%s'", (yyvsp[-4].c));
      Free((yyvsp[-4].c));
    }
#line 9289 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 237:
#line 2993 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool changed = false;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        GModel::current()->getOCCInternals()->remove(dimTags);
        changed = GModel::current()->getOCCInternals()->getChanged();
        if(changed)
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
      }
      else{
        GModel::current()->getGEOInternals()->remove(dimTags);
        changed = GModel::current()->getGEOInternals()->getChanged();
        if(changed)
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
      }
      if(!changed){
        GModel::current()->remove(dimTags);
      }
      List_Delete((yyvsp[-1].l));
    }
#line 9315 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 238:
#line 3015 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      bool changed = false;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        GModel::current()->getOCCInternals()->remove(dimTags, true);
        changed = GModel::current()->getOCCInternals()->getChanged();
        if(changed)
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
      }
      else{
        GModel::current()->getGEOInternals()->remove(dimTags, true);
        changed = GModel::current()->getGEOInternals()->getChanged();
        if(changed)
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
      }
      if(!changed){
        GModel::current()->remove(dimTags, true);
      }
      List_Delete((yyvsp[-1].l));
    }
#line 9341 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 239:
#line 3037 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_MESH)
      GModel::current()->getFields()->deleteField((int)(yyvsp[-2].d));
#endif
    }
#line 9351 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 240:
#line 3043 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-4].c), "View")){
	int index = (int)(yyvsp[-2].d);
	if(index >= 0 && index < (int)PView::list.size())
	  delete PView::list[index];
	else
	  yymsg(0, "Unknown view %d", index);
      }
      else
	yymsg(0, "Unknown command 'Delete %s'", (yyvsp[-4].c));
#endif
      Free((yyvsp[-4].c));
    }
#line 9370 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 241:
#line 3058 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-1].c), "Meshes") || !strcmp((yyvsp[-1].c), "All")){
        ClearProject();
      }
      else if(!strcmp((yyvsp[-1].c), "Model")){
	GModel::current()->destroy(true); // destroy, but keep name/filename
	GModel::current()->getGEOInternals()->destroy();
      }
      else if(!strcmp((yyvsp[-1].c), "Physicals")){
	GModel::current()->getGEOInternals()->resetPhysicalGroups();
	GModel::current()->removePhysicalGroups();
      }
      else if(!strcmp((yyvsp[-1].c), "Variables")){
	gmsh_yysymbols.clear();
      }
      else if(!strcmp((yyvsp[-1].c), "Options")){
        ReInitOptions(0);
        InitOptionsGUI(0);
      }
      else{
	if(gmsh_yysymbols.count((yyvsp[-1].c)))
	  gmsh_yysymbols.erase((yyvsp[-1].c));
	else
	  yymsg(0, "Unknown object or expression to delete '%s'", (yyvsp[-1].c));
      }
      Free((yyvsp[-1].c));
    }
#line 9402 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 242:
#line 3086 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-2].c), "Empty") && !strcmp((yyvsp[-1].c), "Views")){
	for(int i = PView::list.size() - 1; i >= 0; i--)
	  if(PView::list[i]->getData()->empty()) delete PView::list[i];
      }
      else
	yymsg(0, "Unknown command 'Delete %s %s'", (yyvsp[-2].c), (yyvsp[-1].c));
#endif
      Free((yyvsp[-2].c)); Free((yyvsp[-1].c));
    }
#line 9418 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 243:
#line 3098 "Gmsh.y" /* yacc.c:1652  */
    {
      gmsh_yynamespaces.clear();
    }
#line 9426 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 244:
#line 3107 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      setColor(dimTags, (yyvsp[-3].u), false);
      List_Delete((yyvsp[-1].l));
    }
#line 9437 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 245:
#line 3114 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      setColor(dimTags, (yyvsp[-3].u), true);
      List_Delete((yyvsp[-1].l));
    }
#line 9448 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 246:
#line 3126 "Gmsh.y" /* yacc.c:1652  */
    {
      yymsg(2, "'SetPartition' command is deprecated");
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      for(std::size_t i = 0; i < dimTags.size(); i++){
        GEntity *ge = GModel::current()->getEntityByTag
          (dimTags[i].first, dimTags[i].second);
        if(ge){
          for(std::size_t j = 0; j < ge->getNumMeshElements(); j++)
            ge->getMeshElement(j)->setPartition((int)(yyvsp[-3].d));
        }
      }
      List_Delete((yyvsp[-1].l));
    }
#line 9467 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 247:
#line 3146 "Gmsh.y" /* yacc.c:1652  */
    {
      setVisibility(-1, 1, false);
    }
#line 9475 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 248:
#line 3150 "Gmsh.y" /* yacc.c:1652  */
    {
      setVisibility(-1, 1, false);
      Free((yyvsp[-1].c));
    }
#line 9484 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 249:
#line 3155 "Gmsh.y" /* yacc.c:1652  */
    {
      setVisibility(-1, 0, false);
    }
#line 9492 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 250:
#line 3159 "Gmsh.y" /* yacc.c:1652  */
    {
      setVisibility(-1, 0, false);
      Free((yyvsp[-1].c));
    }
#line 9501 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 251:
#line 3164 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      setVisibility(dimTags, 1, false);
      List_Delete((yyvsp[-1].l));
    }
#line 9512 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 252:
#line 3171 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      setVisibility(dimTags, 1, true);
      List_Delete((yyvsp[-1].l));
    }
#line 9523 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 253:
#line 3178 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      setVisibility(dimTags, 0, false);
      List_Delete((yyvsp[-1].l));
    }
#line 9534 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 254:
#line 3185 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > dimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), dimTags);
      setVisibility(dimTags, 0, true);
      List_Delete((yyvsp[-1].l));
    }
#line 9545 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 255:
#line 3197 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-2].c), "Include")){
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
	Msg::StatusBar(true, "Reading '%s'...", tmp.c_str());
	// Warning: we explicitly ask ParseFile not to fclose() the included
        // file, in order to allow user functions definitions in these files.
        // The files will be closed in the next time OpenFile terminates. If
        // you need to include many many files and don't have functions in
        // the files, use "Merge" instead of "Include", as some OSes limit
        // the number of files a process can open simultaneously. (A better
        // solution would be to modify FunctionManager to reopen the files
        // instead of using the FILE pointer...)
	ParseFile(tmp, false, true);
	SetBoundingBox();
	Msg::StatusBar(true, "Done reading '%s'", tmp.c_str());
      }
      else if(!strcmp((yyvsp[-2].c), "Print")){
	// make sure we have the latest data from CAD internals in GModel (fixes
	// bug where we would have no geometry in the picture if the print
	// command is in the same file as the geometry)
        if(GModel::current()->getOCCInternals() &&
           GModel::current()->getOCCInternals()->getChanged())
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
        if(GModel::current()->getGEOInternals()->getChanged())
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
	CreateOutputFile(tmp, CTX::instance()->print.fileFormat);
      }
      else if(!strcmp((yyvsp[-2].c), "Save")){
        if(GModel::current()->getOCCInternals() &&
           GModel::current()->getOCCInternals()->getChanged())
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
        if(GModel::current()->getGEOInternals()->getChanged())
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
	CreateOutputFile(tmp, CTX::instance()->mesh.fileFormat);
      }
      else if(!strcmp((yyvsp[-2].c), "Merge") || !strcmp((yyvsp[-2].c), "MergeWithBoundingBox")){
	// sync CAD internals here, so that if we e.g. import a STEP file, we
        // have the correct entity tags and the numberings don't clash
        if(GModel::current()->getOCCInternals() &&
           GModel::current()->getOCCInternals()->getChanged())
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
        if(GModel::current()->getGEOInternals()->getChanged())
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
        MergeFile(tmp, true);
      }
      else if(!strcmp((yyvsp[-2].c), "NonBlockingSystemCall")){
	SystemCall((yyvsp[-1].c));
      }
      else if(!strcmp((yyvsp[-2].c), "System") || !strcmp((yyvsp[-2].c), "SystemCall")){
	SystemCall((yyvsp[-1].c), true);
      }
      else if(!strcmp((yyvsp[-2].c), "SetName")){
	GModel::current()->setName((yyvsp[-1].c));
      }
      else if(!strcmp((yyvsp[-2].c), "CreateDir")){
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
	CreateSingleDir(tmp);
      }
      else if(!strcmp((yyvsp[-2].c), "OnelabRun")){
        Msg::RunOnelabClient((yyvsp[-1].c));
      }
      else if(!strcmp((yyvsp[-2].c), "OptimizeMesh")){
        GModel::current()->optimizeMesh((yyvsp[-1].c));
      }
      else{
	yymsg(0, "Unknown command '%s'", (yyvsp[-2].c));
      }
      Free((yyvsp[-2].c)); Free((yyvsp[-1].c));
    }
#line 9622 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 256:
#line 3270 "Gmsh.y" /* yacc.c:1652  */
    {
      int n = List_Nbr((yyvsp[-2].l));
      if(n == 1){
        char *s; List_Read((yyvsp[-2].l), 0, &s);
        Msg::RunOnelabClient(s);
        Free(s);
      }
      else if(n == 2){
        char *s, *t; List_Read((yyvsp[-2].l), 0, &s); List_Read((yyvsp[-2].l), 1, &t);
        Msg::RunOnelabClient(s, t);
        Free(s); Free(t);
      }
      else{
        yymsg(0, "OnelabRun takes one or two arguments");
      }
      List_Delete((yyvsp[-2].l));
    }
#line 9644 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 257:
#line 3288 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-5].c), "View")){
	int index = (int)(yyvsp[-3].d);
	if(index >= 0 && index < (int)PView::list.size()){
          if(!strcmp((yyvsp[-6].c), "Save")){
            std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
            PView::list[index]->write(tmp, CTX::instance()->post.fileFormat);
          }
          else if(!strcmp((yyvsp[-6].c), "SendToServer")){
            PView::list[index]->sendToServer((yyvsp[-1].c));
          }
          else{
            yymsg(0, "Unknown operation '%s' on view %d", (yyvsp[-6].c), index);
          }
	}
	else
	  yymsg(0, "Unknown view %d", index);
      }
      else
	yymsg(0, "Unknown command '%s %s'", (yyvsp[-6].c), (yyvsp[-5].c));
#endif
      Free((yyvsp[-6].c)); Free((yyvsp[-5].c)); Free((yyvsp[-1].c));
    }
#line 9673 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 258:
#line 3313 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST) && defined(HAVE_MESH)
      if(!strcmp((yyvsp[-6].c), "Background") && !strcmp((yyvsp[-5].c), "Mesh")  && !strcmp((yyvsp[-4].c), "View")){
	int index = (int)(yyvsp[-2].d);
	if(index >= 0 && index < (int)PView::list.size())
	  GModel::current()->getFields()->setBackgroundMesh(index);
	else
	  yymsg(0, "Unknown view %d", index);
      }
      else
	yymsg(0, "Unknown command '%s'", (yyvsp[-6].c));
#endif
      Free((yyvsp[-6].c)); Free((yyvsp[-5].c)); Free((yyvsp[-4].c));
    }
#line 9692 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 259:
#line 3328 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-2].c), "Sleep")){
	SleepInSeconds((yyvsp[-1].d));
      }
      else if(!strcmp((yyvsp[-2].c), "Remesh")){
	yymsg(0, "Surface remeshing must be reinterfaced");
      }
      else if(!strcmp((yyvsp[-2].c), "Mesh")){
	int lock = CTX::instance()->lock;
	CTX::instance()->lock = 0;
        if(GModel::current()->getOCCInternals() &&
           GModel::current()->getOCCInternals()->getChanged())
          GModel::current()->getOCCInternals()->synchronize(GModel::current());
        if(GModel::current()->getGEOInternals()->getChanged())
          GModel::current()->getGEOInternals()->synchronize(GModel::current());
	GModel::current()->mesh((int)(yyvsp[-1].d));
	CTX::instance()->lock = lock;
      }
      else if(!strcmp((yyvsp[-2].c), "SetOrder")){
#if defined(HAVE_MESH)
        SetOrderN(GModel::current(), (yyvsp[-1].d), CTX::instance()->mesh.secondOrderLinear,
                  CTX::instance()->mesh.secondOrderIncomplete,
                  CTX::instance()->mesh.meshOnlyVisible);
#endif
      }
      else if(!strcmp((yyvsp[-2].c), "PartitionMesh")){
        GModel::current()->partitionMesh((yyvsp[-1].d));
      }
      else
	yymsg(0, "Unknown command '%s'", (yyvsp[-2].c));
      Free((yyvsp[-2].c));
    }
#line 9729 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 260:
#line 3361 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_PLUGINS)
       try {
	 PluginManager::instance()->action((yyvsp[-4].c), (yyvsp[-1].c), 0);
       }
       catch(...) {
	 yymsg(0, "Unknown action '%s' or plugin '%s'", (yyvsp[-1].c), (yyvsp[-4].c));
       }
#endif
       Free((yyvsp[-4].c)); Free((yyvsp[-1].c));
     }
#line 9745 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 261:
#line 3373 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POST)
      if(!strcmp((yyvsp[-1].c), "ElementsFromAllViews"))
	PView::combine(false, 1, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "ElementsFromVisibleViews"))
	PView::combine(false, 0, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "ElementsByViewName"))
	PView::combine(false, 2, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "TimeStepsFromAllViews"))
	PView::combine(true, 1, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "TimeStepsFromVisibleViews"))
	PView::combine(true, 0, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "TimeStepsByViewName"))
	PView::combine(true, 2, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "Views"))
	PView::combine(false, 1, CTX::instance()->post.combineRemoveOrig);
      else if(!strcmp((yyvsp[-1].c), "TimeSteps"))
	PView::combine(true, 2, CTX::instance()->post.combineRemoveOrig);
      else
	yymsg(0, "Unknown 'Combine' command");
#endif
      Free((yyvsp[-1].c));
    }
#line 9773 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 262:
#line 3397 "Gmsh.y" /* yacc.c:1652  */
    {
      Msg::Exit(0);
    }
#line 9781 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 263:
#line 3401 "Gmsh.y" /* yacc.c:1652  */
    {
      gmsh_yyerrorstate = 999; // this will be checked when yyparse returns
      YYABORT;
    }
#line 9790 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 264:
#line 3406 "Gmsh.y" /* yacc.c:1652  */
    {
      // force sync
      if(GModel::current()->getOCCInternals())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      GModel::current()->getGEOInternals()->synchronize(GModel::current());
    }
#line 9801 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 265:
#line 3413 "Gmsh.y" /* yacc.c:1652  */
    {
      new GModel();
      GModel::current(GModel::list.size() - 1);
    }
#line 9810 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 266:
#line 3418 "Gmsh.y" /* yacc.c:1652  */
    {
      CTX::instance()->forcedBBox = 0;
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(GModel::current()->getGEOInternals()->getChanged())
        GModel::current()->getGEOInternals()->synchronize(GModel::current());
      SetBoundingBox();
    }
#line 9824 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 267:
#line 3428 "Gmsh.y" /* yacc.c:1652  */
    {
      CTX::instance()->forcedBBox = 1;
      SetBoundingBox((yyvsp[-12].d), (yyvsp[-10].d), (yyvsp[-8].d), (yyvsp[-6].d), (yyvsp[-4].d), (yyvsp[-2].d));
    }
#line 9833 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 268:
#line 3433 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_OPENGL)
      drawContext::global()->draw();
#endif
    }
#line 9843 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 269:
#line 3439 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_OPENGL)
     CTX::instance()->mesh.changed = ENT_ALL;
     for(std::size_t index = 0; index < PView::list.size(); index++)
       PView::list[index]->setChanged(true);
#endif
    }
#line 9855 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 270:
#line 3447 "Gmsh.y" /* yacc.c:1652  */
    {
      GModel::current()->createTopologyFromMesh();
    }
#line 9863 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 271:
#line 3451 "Gmsh.y" /* yacc.c:1652  */
    {
      GModel::current()->createGeometryOfDiscreteEntities();
    }
#line 9871 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 272:
#line 3455 "Gmsh.y" /* yacc.c:1652  */
    {
      GModel::current()->renumberMeshVertices();
    }
#line 9879 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 273:
#line 3459 "Gmsh.y" /* yacc.c:1652  */
    {
      GModel::current()->renumberMeshElements();
    }
#line 9887 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 274:
#line 3463 "Gmsh.y" /* yacc.c:1652  */
    {
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(GModel::current()->getGEOInternals()->getChanged())
        GModel::current()->getGEOInternals()->synchronize(GModel::current());
      GModel::current()->refineMesh(CTX::instance()->mesh.secondOrderLinear);
    }
#line 9900 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 275:
#line 3473 "Gmsh.y" /* yacc.c:1652  */
    {
      int lock = CTX::instance()->lock;
      CTX::instance()->lock = 0;
      std::vector<int> technique;
      for(int i = 0; i < List_Nbr((yyvsp[-13].l)); i++){
        double d;
        List_Read((yyvsp[-13].l), i, &d);
        technique.push_back((int)d);
      }
      if(technique.empty()){
        yymsg(0, "Need at least one adaptation technique");
      }
      else{
        std::vector<simpleFunction<double>*> f;
        for(int i = 0; i < List_Nbr((yyvsp[-10].l)); i++){
          double d;
          List_Read((yyvsp[-10].l), i, &d);
          gLevelset *l = gLevelset::find((int)d);
          if(l) f.push_back(l);
          else yymsg(0, "Unknown levelset %d", (int)d);
        }
        if(technique.size() != f.size()){
          yymsg(0, "Number of techniques != number of levelsets");
        }
        else{
          if(List_Nbr((yyvsp[-7].l)) != (int)f.size()){
            yymsg(0, "Number of parameters != number of levelsets");
          }
          else{
            std::vector<std::vector<double> > parameters;
            parameters.resize(List_Nbr((yyvsp[-7].l)));
            for(int i = 0; i < List_Nbr((yyvsp[-7].l)); i++){
              List_T *l = *(List_T**)List_Pointer((yyvsp[-7].l), i);
              for(int j = 0; j < List_Nbr(l); j++){
                double d;
                List_Read(l, j, &d);
                parameters[i].push_back(d);
              }
            }
            int niter = (int)(yyvsp[-4].d);
            bool meshAll = ((yyvsp[-2].d) == 0) ? false : true;
            if(GModel::current()->getOCCInternals() &&
               GModel::current()->getOCCInternals()->getChanged())
              GModel::current()->getOCCInternals()->synchronize(GModel::current());
            if(GModel::current()->getGEOInternals()->getChanged())
              GModel::current()->getGEOInternals()->synchronize(GModel::current());
            GModel::current()->adaptMesh(technique, f, parameters, niter, meshAll);
          }
        }
      }
      List_Delete((yyvsp[-13].l));
      List_Delete((yyvsp[-10].l));
      for(int i = 0; i < List_Nbr((yyvsp[-7].l)); i++)
        List_Delete(*(List_T**)List_Pointer((yyvsp[-7].l), i));
      List_Delete((yyvsp[-7].l));
      CTX::instance()->lock = lock;
    }
#line 9962 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 276:
#line 3536 "Gmsh.y" /* yacc.c:1652  */
    {
#if defined(HAVE_POPPLER)
       std::vector<int> is;
       for(int i = 0; i < List_Nbr((yyvsp[-7].l)); i++){
	 double d;
	 List_Read((yyvsp[-7].l), i, &d);
	 is.push_back ((int) d);
       }
       gmshPopplerWrapper::instance()->setMacroForPages(is, (yyvsp[-4].c), (yyvsp[-2].c) );
#endif
     }
#line 9978 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 277:
#line 3552 "Gmsh.y" /* yacc.c:1652  */
    {
      LoopControlVariablesTab[ImbricatedLoop][0] = (yyvsp[-3].d);
      LoopControlVariablesTab[ImbricatedLoop][1] = (yyvsp[-1].d);
      LoopControlVariablesTab[ImbricatedLoop][2] = 1.0;
      LoopControlVariablesNameTab[ImbricatedLoop] = "";
      fgetpos(gmsh_yyin, &yyposImbricatedLoopsTab[ImbricatedLoop]);
      yylinenoImbricatedLoopsTab[ImbricatedLoop] = gmsh_yylineno;
      if((yyvsp[-3].d) > (yyvsp[-1].d))
	skip("For", "EndFor");
      else
	ImbricatedLoop++;
      if(ImbricatedLoop > MAX_RECUR_LOOPS - 1){
	yymsg(0, "Reached maximum number of imbricated loops");
	ImbricatedLoop = MAX_RECUR_LOOPS - 1;
      }
    }
#line 9999 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 278:
#line 3569 "Gmsh.y" /* yacc.c:1652  */
    {
      LoopControlVariablesTab[ImbricatedLoop][0] = (yyvsp[-5].d);
      LoopControlVariablesTab[ImbricatedLoop][1] = (yyvsp[-3].d);
      LoopControlVariablesTab[ImbricatedLoop][2] = (yyvsp[-1].d);
      LoopControlVariablesNameTab[ImbricatedLoop] = "";
      fgetpos(gmsh_yyin, &yyposImbricatedLoopsTab[ImbricatedLoop]);
      yylinenoImbricatedLoopsTab[ImbricatedLoop] = gmsh_yylineno;
      if(((yyvsp[-1].d) > 0. && (yyvsp[-5].d) > (yyvsp[-3].d)) || ((yyvsp[-1].d) < 0. && (yyvsp[-5].d) < (yyvsp[-3].d)))
	skip("For", "EndFor");
      else
	ImbricatedLoop++;
      if(ImbricatedLoop > MAX_RECUR_LOOPS - 1){
	yymsg(0, "Reached maximum number of imbricated loops");
	ImbricatedLoop = MAX_RECUR_LOOPS - 1;
      }
    }
#line 10020 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 279:
#line 3586 "Gmsh.y" /* yacc.c:1652  */
    {
      LoopControlVariablesTab[ImbricatedLoop][0] = (yyvsp[-3].d);
      LoopControlVariablesTab[ImbricatedLoop][1] = (yyvsp[-1].d);
      LoopControlVariablesTab[ImbricatedLoop][2] = 1.0;
      LoopControlVariablesNameTab[ImbricatedLoop] = (yyvsp[-6].c);
      gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-6].c)]);
      s.list = false;
      s.value.resize(1);
      s.value[0] = (yyvsp[-3].d);
      fgetpos(gmsh_yyin, &yyposImbricatedLoopsTab[ImbricatedLoop]);
      yylinenoImbricatedLoopsTab[ImbricatedLoop] = gmsh_yylineno;
      if((yyvsp[-3].d) > (yyvsp[-1].d))
	skip("For", "EndFor");
      else
	ImbricatedLoop++;
      if(ImbricatedLoop > MAX_RECUR_LOOPS - 1){
	yymsg(0, "Reached maximum number of imbricated loops");
	ImbricatedLoop = MAX_RECUR_LOOPS - 1;
      }
      Free((yyvsp[-6].c));
    }
#line 10046 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 280:
#line 3608 "Gmsh.y" /* yacc.c:1652  */
    {
      LoopControlVariablesTab[ImbricatedLoop][0] = (yyvsp[-5].d);
      LoopControlVariablesTab[ImbricatedLoop][1] = (yyvsp[-3].d);
      LoopControlVariablesTab[ImbricatedLoop][2] = (yyvsp[-1].d);
      LoopControlVariablesNameTab[ImbricatedLoop] = (yyvsp[-8].c);
      gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-8].c)]);
      s.list = false;
      s.value.resize(1);
      s.value[0] = (yyvsp[-5].d);
      fgetpos(gmsh_yyin, &yyposImbricatedLoopsTab[ImbricatedLoop]);
      yylinenoImbricatedLoopsTab[ImbricatedLoop] = gmsh_yylineno;
      if(((yyvsp[-1].d) > 0. && (yyvsp[-5].d) > (yyvsp[-3].d)) || ((yyvsp[-1].d) < 0. && (yyvsp[-5].d) < (yyvsp[-3].d)))
	skip("For", "EndFor");
      else
	ImbricatedLoop++;
      if(ImbricatedLoop > MAX_RECUR_LOOPS - 1){
	yymsg(0, "Reached maximum number of imbricated loops");
	ImbricatedLoop = MAX_RECUR_LOOPS - 1;
      }
      Free((yyvsp[-8].c));
    }
#line 10072 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 281:
#line 3630 "Gmsh.y" /* yacc.c:1652  */
    {
      if(ImbricatedLoop <= 0){
	yymsg(0, "Invalid For/EndFor loop");
	ImbricatedLoop = 0;
      }
      else{
	double step = LoopControlVariablesTab[ImbricatedLoop - 1][2];
        std::string name = LoopControlVariablesNameTab[ImbricatedLoop - 1];
        if(name.size()){
          if(!gmsh_yysymbols.count(name))
            yymsg(0, "Unknown loop variable '%s'", name.c_str());
          else{
            gmsh_yysymbol &s(gmsh_yysymbols[name]);
            if(!s.list && s.value.size()){
              s.value[0] += step;
              LoopControlVariablesTab[ImbricatedLoop - 1][0] = s.value[0];
            }
            else
              yymsg(0, "Bad loop variable %s", name.c_str());
          }
        }
        else{
          LoopControlVariablesTab[ImbricatedLoop - 1][0] += step;
        }
	double x0 = LoopControlVariablesTab[ImbricatedLoop - 1][0];
	double x1 = LoopControlVariablesTab[ImbricatedLoop - 1][1];
        if((step > 0. && x0 <= x1) || (step < 0. && x0 >= x1)){
	  fsetpos(gmsh_yyin, &yyposImbricatedLoopsTab[ImbricatedLoop - 1]);
	  gmsh_yylineno = yylinenoImbricatedLoopsTab[ImbricatedLoop - 1];
	}
	else
	  ImbricatedLoop--;
      }
    }
#line 10111 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 282:
#line 3665 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!FunctionManager::Instance()->createFunction
         (std::string((yyvsp[0].c)), gmsh_yyin, gmsh_yyname, gmsh_yylineno))
	yymsg(0, "Redefinition of function %s", (yyvsp[0].c));
      skip(nullptr, "Return");
      Free((yyvsp[0].c));
    }
#line 10123 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 283:
#line 3673 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!FunctionManager::Instance()->createFunction
         (std::string((yyvsp[0].c)), gmsh_yyin, gmsh_yyname, gmsh_yylineno))
	yymsg(0, "Redefinition of function %s", (yyvsp[0].c));
      skip(nullptr, "Return");
      Free((yyvsp[0].c));
    }
#line 10135 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 284:
#line 3681 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!FunctionManager::Instance()->leaveFunction
         (&gmsh_yyin, gmsh_yyname, gmsh_yylineno))
	yymsg(0, "Error while exiting function");
    }
#line 10145 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 285:
#line 3687 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!FunctionManager::Instance()->enterFunction
         (std::string((yyvsp[-1].c)), &gmsh_yyin, gmsh_yyname, gmsh_yylineno))
	yymsg(0, "Unknown function '%s'", (yyvsp[-1].c));
      Free((yyvsp[-1].c));
    }
#line 10156 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 286:
#line 3694 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!FunctionManager::Instance()->enterFunction
         (std::string((yyvsp[-1].c)), &gmsh_yyin, gmsh_yyname, gmsh_yylineno))
	yymsg(0, "Unknown function '%s'", (yyvsp[-1].c));
      Free((yyvsp[-1].c));
    }
#line 10167 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 287:
#line 3701 "Gmsh.y" /* yacc.c:1652  */
    {
      ImbricatedTest++;
      if(ImbricatedTest > MAX_RECUR_TESTS-1){
	yymsg(0, "Reached maximum number of imbricated tests");
        ImbricatedTest = MAX_RECUR_TESTS-1;
      }

      if((yyvsp[-1].d)){
        // Current test is true
        statusImbricatedTests[ImbricatedTest] = 1;
      }
      else{
        statusImbricatedTests[ImbricatedTest] = 0;
        // Go after the next ElseIf or Else or EndIf
        int type_until2 = 0;
        skipTest("If", "EndIf", "ElseIf", 4, &type_until2);
        if(!type_until2) ImbricatedTest--; // EndIf reached
      }
    }
#line 10191 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 288:
#line 3721 "Gmsh.y" /* yacc.c:1652  */
    {
      if(ImbricatedTest > 0){
        if (statusImbricatedTests[ImbricatedTest]){
          // Last test (If or ElseIf) was true, thus go after EndIf (out of If EndIf)
          skip("If", "EndIf");
          ImbricatedTest--;
        }
        else{
          // Previous test(s) (If and ElseIf) not yet true
          if((yyvsp[-1].d)){
            statusImbricatedTests[ImbricatedTest] = 1;
          }
          else{
            // Current test still not true: statusImbricatedTests[ImbricatedTest] = 0;
            // Go after the next ElseIf or Else or EndIf
            int type_until2 = 0;
            skipTest("If", "EndIf", "ElseIf", 4, &type_until2);
            if(!type_until2) ImbricatedTest--;
          }
        }
      }
      else{
	yymsg(0, "Orphan ElseIf");
      }
    }
#line 10221 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 289:
#line 3747 "Gmsh.y" /* yacc.c:1652  */
    {
      if(ImbricatedTest > 0){
        if(statusImbricatedTests[ImbricatedTest]){
          skip("If", "EndIf");
          ImbricatedTest--;
        }
      }
      else{
	yymsg(0, "Orphan Else");
      }
    }
#line 10237 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 290:
#line 3759 "Gmsh.y" /* yacc.c:1652  */
    {
      ImbricatedTest--;
      if(ImbricatedTest < 0)
        yymsg(1, "Orphan EndIf");
    }
#line 10247 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 291:
#line 3770 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->extrude
          (inDimTags, (yyvsp[-3].v)[0], (yyvsp[-3].v)[1], (yyvsp[-3].v)[2], outDimTags);
      }
      else{
        r = GModel::current()->getGEOInternals()->extrude
          (inDimTags, (yyvsp[-3].v)[0], (yyvsp[-3].v)[1], (yyvsp[-3].v)[2], outDimTags);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-1].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10269 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 292:
#line 3788 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->revolve
          (inDimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-4].d), outDimTags);
      }
      else{
        r = GModel::current()->getGEOInternals()->revolve
          (inDimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-4].d), outDimTags);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-1].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10291 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 293:
#line 3806 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-1].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        yymsg(0, "Twisting extrude not available with OpenCASCADE geometry kernel");
      }
      else{
        r = GModel::current()->getGEOInternals()->twist
          (inDimTags, (yyvsp[-6].v)[0], (yyvsp[-6].v)[1], (yyvsp[-6].v)[2], (yyvsp[-10].v)[0], (yyvsp[-10].v)[1], (yyvsp[-10].v)[2], (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2],
           (yyvsp[-4].d), outDimTags);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-1].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10313 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 294:
#line 3824 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.ExtrudeMesh = extr.mesh.Recombine = false;
      extr.mesh.QuadToTri = NO_QUADTRI;
      extr.mesh.ScaleLast = false;
    }
#line 10323 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 295:
#line 3830 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-3].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->extrude
          (inDimTags, (yyvsp[-5].v)[0], (yyvsp[-5].v)[1], (yyvsp[-5].v)[2], outDimTags, &extr);
      }
      else{
        r = GModel::current()->getGEOInternals()->extrude
          (inDimTags, (yyvsp[-5].v)[0], (yyvsp[-5].v)[1], (yyvsp[-5].v)[2], outDimTags, &extr);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-3].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10345 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 296:
#line 3848 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.ExtrudeMesh = extr.mesh.Recombine = false;
      extr.mesh.QuadToTri = NO_QUADTRI;
      extr.mesh.ScaleLast = false;
    }
#line 10355 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 297:
#line 3854 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-3].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->revolve
          (inDimTags, (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-10].v)[0], (yyvsp[-10].v)[1], (yyvsp[-10].v)[2], (yyvsp[-6].d), outDimTags,
           &extr);
      }
      else{
        r = GModel::current()->getGEOInternals()->revolve
          (inDimTags, (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-10].v)[0], (yyvsp[-10].v)[1], (yyvsp[-10].v)[2], (yyvsp[-6].d), outDimTags,
           &extr);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-3].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10379 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 298:
#line 3874 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.ExtrudeMesh = extr.mesh.Recombine = false;
      extr.mesh.QuadToTri = NO_QUADTRI;
      extr.mesh.ScaleLast = false;
    }
#line 10389 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 299:
#line 3880 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-3].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        yymsg(0, "Twisting extrude not available with OpenCASCADE geometry kernel");
      }
      else{
        r = GModel::current()->getGEOInternals()->twist
          (inDimTags, (yyvsp[-8].v)[0], (yyvsp[-8].v)[1], (yyvsp[-8].v)[2], (yyvsp[-12].v)[0], (yyvsp[-12].v)[1], (yyvsp[-12].v)[2], (yyvsp[-10].v)[0], (yyvsp[-10].v)[1], (yyvsp[-10].v)[2],
           (yyvsp[-6].d),  outDimTags, &extr);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-3].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10411 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 300:
#line 3898 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.ExtrudeMesh = extr.mesh.Recombine = false;
      extr.mesh.QuadToTri = NO_QUADTRI;
      extr.mesh.ScaleLast = false;
    }
#line 10421 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 301:
#line 3904 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-3].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        yymsg(0, "Boundary layer extrusion not available with OpenCASCADE geometry kernel");
      }
      else{
        r = GModel::current()->getGEOInternals()->boundaryLayer
          (inDimTags, outDimTags, &extr);
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-3].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10442 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 302:
#line 3921 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<std::pair<int, int> > inDimTags, outDimTags;
      ListOfShapes2VectorOfPairs((yyvsp[-6].l), inDimTags);
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        r = GModel::current()->getOCCInternals()->addPipe(inDimTags, (int)(yyvsp[-1].d), outDimTags);
      }
      else{
        yymsg(0, "Pipe only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not extrude shapes");
      (yyval.l) = (yyvsp[-6].l);
      List_Reset((yyval.l));
      VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
    }
#line 10462 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 303:
#line 3937 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<int> wires; ListOfDouble2Vector((yyvsp[0].l), wires);
        std::vector<std::pair<int, int> > outDimTags;
        r = GModel::current()->getOCCInternals()->addThruSections
          (-1, wires, false, false, outDimTags);
        VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
      }
      else{
        yymsg(0, "ThruSections only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add thrusections");
      List_Delete((yyvsp[0].l));
    }
#line 10483 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 304:
#line 3954 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<int> wires; ListOfDouble2Vector((yyvsp[0].l), wires);
        std::vector<std::pair<int, int> > outDimTags;
        r = GModel::current()->getOCCInternals()->addThruSections
          (-1, wires, false, true, outDimTags);
        VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
      }
      else{
        yymsg(0, "ThruSections only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not add ruled thrusections");
      List_Delete((yyvsp[0].l));
    }
#line 10504 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 305:
#line 3972 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<int> regions, edges;
        ListOfDouble2Vector((yyvsp[-7].l), regions); ListOfDouble2Vector((yyvsp[-4].l), edges);
        std::vector<double> radii;
        ListOfDouble2Vector((yyvsp[-1].l), radii);
        std::vector<std::pair<int, int> > outDimTags;
        r = GModel::current()->getOCCInternals()->fillet
          (regions, edges, radii, outDimTags, true);
        VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
      }
      else{
        yymsg(0, "Fillet only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not fillet shapes");
      List_Delete((yyvsp[-7].l));
      List_Delete((yyvsp[-4].l));
      List_Delete((yyvsp[-1].l));
    }
#line 10530 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 306:
#line 3995 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<int> regions, edges, surfaces;
        ListOfDouble2Vector((yyvsp[-10].l), regions); ListOfDouble2Vector((yyvsp[-7].l), edges);
        ListOfDouble2Vector((yyvsp[-4].l), surfaces);
        std::vector<double> distances;
        ListOfDouble2Vector((yyvsp[-1].l), distances);
        std::vector<std::pair<int, int> > outDimTags;
        r = GModel::current()->getOCCInternals()->chamfer
          (regions, edges, surfaces, distances, outDimTags, true);
        VectorOfPairs2ListOfShapes(outDimTags, (yyval.l));
      }
      else{
        yymsg(0, "Chamfer only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not chamfer shapes");
      List_Delete((yyvsp[-10].l));
      List_Delete((yyvsp[-7].l));
      List_Delete((yyvsp[-4].l));
      List_Delete((yyvsp[-1].l));
    }
#line 10558 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 307:
#line 4022 "Gmsh.y" /* yacc.c:1652  */
    {
    }
#line 10565 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 308:
#line 4025 "Gmsh.y" /* yacc.c:1652  */
    {
    }
#line 10572 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 309:
#line 4031 "Gmsh.y" /* yacc.c:1652  */
    {
      int n = (int)fabs((yyvsp[-2].d));
      if(n){ // we accept n==0 to easily disable layers
        extr.mesh.ExtrudeMesh = true;
        extr.mesh.NbLayer = 1;
        extr.mesh.NbElmLayer.clear();
        extr.mesh.hLayer.clear();
        extr.mesh.NbElmLayer.push_back((int)fabs((yyvsp[-2].d)));
        extr.mesh.hLayer.push_back(1.);
      }
    }
#line 10588 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 310:
#line 4043 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.ExtrudeMesh = true;
      extr.mesh.NbLayer = List_Nbr((yyvsp[-4].l));
      if(List_Nbr((yyvsp[-4].l)) == List_Nbr((yyvsp[-2].l))){
	extr.mesh.NbElmLayer.clear();
	extr.mesh.hLayer.clear();
	for(int i = 0; i < List_Nbr((yyvsp[-4].l)); i++){
	  double d;
	  List_Read((yyvsp[-4].l), i, &d);
	  extr.mesh.NbElmLayer.push_back((d > 0) ? (int)d : 1);
	  List_Read((yyvsp[-2].l), i, &d);
	  extr.mesh.hLayer.push_back(d);
	}
      }
      else
	yymsg(0, "Wrong layer definition {%d, %d}", List_Nbr((yyvsp[-4].l)), List_Nbr((yyvsp[-2].l)));
      List_Delete((yyvsp[-4].l));
      List_Delete((yyvsp[-2].l));
    }
#line 10612 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 311:
#line 4063 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.ScaleLast = true;
    }
#line 10620 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 312:
#line 4067 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.Recombine = true;
    }
#line 10628 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 313:
#line 4071 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.Recombine = (yyvsp[-1].d) ? true : false;
    }
#line 10636 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 314:
#line 4075 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.QuadToTri = QUADTRI_ADDVERTS_1;
    }
#line 10644 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 315:
#line 4079 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.QuadToTri = QUADTRI_ADDVERTS_1_RECOMB;
    }
#line 10652 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 316:
#line 4083 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.QuadToTri = QUADTRI_NOVERTS_1;
    }
#line 10660 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 317:
#line 4087 "Gmsh.y" /* yacc.c:1652  */
    {
      extr.mesh.QuadToTri = QUADTRI_NOVERTS_1_RECOMB;
    }
#line 10668 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 318:
#line 4091 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-3].l), tags);
      int num = (int)(yyvsp[-6].d);
      GModel::current()->getGEOInternals()->addDiscreteSurface(num);
      extr.mesh.Holes[num].first = (yyvsp[-1].d);
      extr.mesh.Holes[num].second = tags;
      List_Delete((yyvsp[-3].l));
    }
#line 10681 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 319:
#line 4100 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-4].c), "Index"))
        extr.mesh.BoundaryLayerIndex = (yyvsp[-2].d);
      else if(!strcmp((yyvsp[-4].c), "View"))
        extr.mesh.ViewIndex = (yyvsp[-2].d);
      Free((yyvsp[-4].c));
    }
#line 10693 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 320:
#line 4112 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = OCC_Internals::Union; }
#line 10699 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 321:
#line 4113 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = OCC_Internals::Intersection; }
#line 10705 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 322:
#line 4114 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = OCC_Internals::Difference; }
#line 10711 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 323:
#line 4115 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = OCC_Internals::Section; }
#line 10717 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 324:
#line 4116 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = OCC_Internals::Fragments; }
#line 10723 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 325:
#line 4120 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 0; }
#line 10729 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 326:
#line 4121 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 1; }
#line 10735 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 327:
#line 4122 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 2; }
#line 10741 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 328:
#line 4123 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = (yyvsp[-1].d) ? 1 : 0; }
#line 10747 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 329:
#line 4124 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = (yyvsp[-1].d) ? 2 : 0; }
#line 10753 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 330:
#line 4129 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<std::pair<int, int > > object, tool, out;
        std::vector<std::vector<std::pair<int, int > > > outMap;
        ListOfShapes2VectorOfPairs((yyvsp[-6].l), object);
        ListOfShapes2VectorOfPairs((yyvsp[-2].l), tool);
        // currently we don't distinguish between Delete and Recursive Delete:
        // we always delete recursively. Let us know if you have examples where
        // having the choice would be interesting
        r = GModel::current()->getOCCInternals()->booleanOperator
          (-1, (OCC_Internals::BooleanOperator)(yyvsp[-8].i), object, tool, out, outMap, (yyvsp[-5].i), (yyvsp[-1].i));
        VectorOfPairs2ListOfShapes(out, (yyval.l));
      }
      else{
        yymsg(0, "Boolean operators only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could not apply boolean operator");
      List_Delete((yyvsp[-6].l));
      List_Delete((yyvsp[-2].l));
    }
#line 10780 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 331:
#line 4152 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(Shape));
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<std::pair<int, int> > out;
        std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
        GModel::current()->getOCCInternals()->importShapes(tmp, true, out);
        VectorOfPairs2ListOfShapes(out, (yyval.l));
      }
      else{
        yymsg(0, "ShapeFromFile only available with OpenCASCADE geometry kernel");
      }
      if(!r) yymsg(0, "Could import shape");
      Free((yyvsp[-1].c));
    }
#line 10800 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 332:
#line 4172 "Gmsh.y" /* yacc.c:1652  */
    {
      bool r = true;
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals()){
        std::vector<std::pair<int, int> > object, tool, out;
        std::vector<std::vector<std::pair<int, int > > > outMap;
        ListOfShapes2VectorOfPairs((yyvsp[-7].l), object);
        ListOfShapes2VectorOfPairs((yyvsp[-3].l), tool);
        // currently we don't distinguish between Delete and Recursive Delete:
        // we always delete recursively. Let us know if you have examples where
        // having the choice would be interesting
        r = GModel::current()->getOCCInternals()->booleanOperator
          ((int)(yyvsp[-11].d), (OCC_Internals::BooleanOperator)(yyvsp[-13].i), object, tool, out, outMap, (yyvsp[-6].i), (yyvsp[-2].i));
      }
      if(!r) yymsg(0, "Could not apply boolean operator");
      List_Delete((yyvsp[-7].l));
      List_Delete((yyvsp[-3].l));
    }
#line 10822 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 333:
#line 4193 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.v)[0] = (yyval.v)[1] = 1.;
    }
#line 10830 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 334:
#line 4197 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-1].c), "Progression") || !strcmp((yyvsp[-1].c), "Power"))
        (yyval.v)[0] = 1.;
      else if(!strcmp((yyvsp[-1].c), "Bump"))
        (yyval.v)[0] = 2.;
      else{
        yymsg(0, "Unknown transfinite mesh type");
        (yyval.v)[0] = 1.;
      }
      (yyval.v)[1] = (yyvsp[0].d);
      Free((yyvsp[-1].c));
    }
#line 10847 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 335:
#line 4212 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = -1; // left
    }
#line 10855 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 336:
#line 4216 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[0].c), "Right"))
        (yyval.i) = 1;
      else if(!strcmp((yyvsp[0].c), "Left"))
	(yyval.i) = -1;
      else if(!strcmp((yyvsp[0].c), "AlternateRight"))
	(yyval.i) = 2;
      else if(!strcmp((yyvsp[0].c), "AlternateLeft"))
	(yyval.i) = -2;
      else // "Alternate" -> "Alternate Right"
	(yyval.i) = 2;
      Free((yyvsp[0].c));
    }
#line 10873 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 337:
#line 4232 "Gmsh.y" /* yacc.c:1652  */
    {
     (yyval.l) = List_Create(1, 1, sizeof(double));
   }
#line 10881 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 338:
#line 4236 "Gmsh.y" /* yacc.c:1652  */
    {
     (yyval.l) = (yyvsp[0].l);
   }
#line 10889 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 339:
#line 4241 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = 45;
    }
#line 10897 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 340:
#line 4245 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.i) = (int)(yyvsp[0].d);
    }
#line 10905 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 341:
#line 4251 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(1, 1, sizeof(double));
    }
#line 10913 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 342:
#line 4255 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[0].l);
    }
#line 10921 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 343:
#line 4262 "Gmsh.y" /* yacc.c:1652  */
    {
      // mesh sizes at vertices are stored in internal CAD data, as they can be
      // specified during vertex creation and copied around during CAD
      // operations
      List_T *tmp = (yyvsp[-3].l);
      if(!(yyvsp[-3].l)){
        tmp = List_Create(100, 100, sizeof(double));
        getAllElementaryTags(0, tmp);
      }
      for(int i = 0; i < List_Nbr(tmp); i++){
        double d;
        List_Read(tmp, i, &d);
        int tag = (int)d;
        if(GModel::current()->getOCCInternals())
          GModel::current()->getOCCInternals()->setMeshSize(0, tag, (yyvsp[-1].d));
        GModel::current()->getGEOInternals()->setMeshSize(0, tag, (yyvsp[-1].d));
        GVertex *gv = GModel::current()->getVertexByTag(tag);
        if(gv) gv->setPrescribedMeshSizeAtVertex((yyvsp[-1].d));
      }
      List_Delete(tmp);
    }
#line 10947 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 344:
#line 4284 "Gmsh.y" /* yacc.c:1652  */
    {
      // transfinite constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      int type = (int)(yyvsp[-1].v)[0];
      double coef = fabs((yyvsp[-1].v)[1]);
      int npoints = ((int)(yyvsp[-2].d) < 2) ? 2 : (int)(yyvsp[-2].d);
      if(!(yyvsp[-4].l)){
        GModel::current()->getGEOInternals()->setTransfiniteLine
          (0, npoints, type, coef);
        for(GModel::eiter it = GModel::current()->firstEdge();
            it != GModel::current()->lastEdge(); it++){
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          (*it)->meshAttributes.nbPointsTransfinite = npoints;
          (*it)->meshAttributes.typeTransfinite = type;
          (*it)->meshAttributes.coeffTransfinite = coef;
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-4].l)); i++){
          double d;
          List_Read((yyvsp[-4].l), i, &d);
          int j = (int)fabs(d);
          for(int sig = -1; sig <= 1; sig += 2){
            GModel::current()->getGEOInternals()->setTransfiniteLine
              (sig * j, npoints, type * gmsh_sign(d), coef);
            GEdge *ge = GModel::current()->getEdgeByTag(sig * j);
            if(ge){
              ge->meshAttributes.method = MESH_TRANSFINITE;
              ge->meshAttributes.nbPointsTransfinite = npoints;
              ge->meshAttributes.typeTransfinite = type * gmsh_sign(d);
              ge->meshAttributes.coeffTransfinite = coef;
            }
          }
        }
        List_Delete((yyvsp[-4].l));
      }
    }
#line 10992 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 345:
#line 4325 "Gmsh.y" /* yacc.c:1652  */
    {
      // transfinite constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      std::vector<int> corners; ListOfDouble2Vector((yyvsp[-2].l), corners);
      if(!(yyvsp[-3].l)){
        GModel::current()->getGEOInternals()->setTransfiniteSurface(0, (yyvsp[-1].i), corners);
        for(GModel::fiter it = GModel::current()->firstFace();
            it != GModel::current()->lastFace(); it++){
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          (*it)->meshAttributes.transfiniteArrangement = (yyvsp[-1].i);
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-3].l)); i++){
          double d;
          List_Read((yyvsp[-3].l), i, &d);
          int tag = (int)d;
          GModel::current()->getGEOInternals()->setTransfiniteSurface(tag, (yyvsp[-1].i), corners);
          GFace *gf = GModel::current()->getFaceByTag(tag);
          if(gf){
            gf->meshAttributes.method = MESH_TRANSFINITE;
            gf->meshAttributes.transfiniteArrangement = (yyvsp[-1].i);
            if(corners.empty() || corners.size() == 3 || corners.size() == 4){
              for(std::size_t j = 0; j < corners.size(); j++){
                GVertex *gv = GModel::current()->getVertexByTag(corners[j]);
                if(gv)
                  gf->meshAttributes.corners.push_back(gv);
                else
                  yymsg(0, "Unknown model point with tag %d", corners[j]);
              }
            }
            else{
              yymsg(0, "Transfinite surface requires 3 or 4 corners vertices");
            }
          }
        }
        List_Delete((yyvsp[-3].l));
      }
      List_Delete((yyvsp[-2].l));
    }
#line 11040 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 346:
#line 4369 "Gmsh.y" /* yacc.c:1652  */
    {
      // transfinite constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      std::vector<int> corners; ListOfDouble2Vector((yyvsp[-1].l), corners);
      if(!(yyvsp[-2].l)){
        GModel::current()->getGEOInternals()->setTransfiniteVolume(0, corners);
        for(GModel::riter it = GModel::current()->firstRegion();
            it != GModel::current()->lastRegion(); it++){
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
          double d;
          List_Read((yyvsp[-2].l), i, &d);
          int tag = (int)d;
          GModel::current()->getGEOInternals()->setTransfiniteVolume(tag, corners);
          GRegion *gr = GModel::current()->getRegionByTag(tag);
          if(gr){
            gr->meshAttributes.method = MESH_TRANSFINITE;
            if(corners.empty() || corners.size() == 6 || corners.size() == 8){
              for(std::size_t i = 0; i < corners.size(); i++){
                GVertex *gv = GModel::current()->getVertexByTag(corners[i]);
                if(gv)
                  gr->meshAttributes.corners.push_back(gv);
                else
                  yymsg(0, "Unknown model point with tag %d", corners[i]);
              }
            }
          }
        }
        List_Delete((yyvsp[-2].l));
      }
      List_Delete((yyvsp[-1].l));
    }
#line 11083 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 347:
#line 4408 "Gmsh.y" /* yacc.c:1652  */
    {
      // transfinite constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(!(yyvsp[-1].l)){
        GModel::current()->getGEOInternals()->setTransfiniteVolumeQuadTri(0);
        for(GModel::riter it = GModel::current()->firstRegion();
            it != GModel::current()->lastRegion(); it++)
          (*it)->meshAttributes.QuadTri = TRANSFINITE_QUADTRI_1;
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
          double d;
          List_Read((yyvsp[-1].l), i, &d);
          int tag = (int)d;
          GModel::current()->getGEOInternals()->setTransfiniteVolumeQuadTri(tag);
          GRegion *gr = GModel::current()->getRegionByTag(tag);
          if(gr) gr->meshAttributes.QuadTri = TRANSFINITE_QUADTRI_1;
        }
        List_Delete((yyvsp[-1].l));
      }
    }
#line 11112 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 348:
#line 4433 "Gmsh.y" /* yacc.c:1652  */
    {
      int tag = (int)(yyvsp[-4].d);
      GVertex *gf = GModel::current()->getVertexByTag(tag);
      if(gf){
	int new_tag = (int)(yyvsp[-2].d);
	gf->setTag(new_tag);
      }
      else{
	yymsg(0, "Unknown Model Vertex %d",tag);
      }
    }
#line 11128 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 349:
#line 4445 "Gmsh.y" /* yacc.c:1652  */
    {
      int tag = (int)(yyvsp[-4].d);
      GEdge *gf = GModel::current()->getEdgeByTag(tag);
      if(gf){
	int new_tag = (int)(yyvsp[-2].d);
	gf->setTag(new_tag);
      }
      else{
	yymsg(0, "Unknown Model Edge %d",tag);
      }
    }
#line 11144 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 350:
#line 4457 "Gmsh.y" /* yacc.c:1652  */
    {
      int tag = (int)(yyvsp[-4].d);
      GFace *gf = GModel::current()->getFaceByTag(tag);
      if(gf){
	int new_tag = (int)(yyvsp[-2].d);
	gf->setTag(new_tag);
      }
      else{
	yymsg(0, "Unknown Model Face %d",tag);
      }
    }
#line 11160 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 351:
#line 4469 "Gmsh.y" /* yacc.c:1652  */
    {
      int tag = (int)(yyvsp[-4].d);
      GRegion *gf = GModel::current()->getRegionByTag(tag);
      if(gf){
	int new_tag = (int)(yyvsp[-2].d);
	gf->setTag(new_tag);
      }
      else{
	yymsg(0, "Unknown Model Region %d",tag);
      }
    }
#line 11176 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 352:
#line 4481 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < List_Nbr((yyvsp[-4].l)); i++){
	double d;
	List_Read((yyvsp[-4].l), i, &d);
	CTX::instance()->mesh.algo2dPerFace[(int)d] = (int)(yyvsp[-1].d);
      }
      List_Delete((yyvsp[-4].l));
    }
#line 11189 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 353:
#line 4490 "Gmsh.y" /* yacc.c:1652  */
    {
      // recombine constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(!(yyvsp[-2].l)){
        GModel::current()->getGEOInternals()->setRecombine(2, 0, (yyvsp[-1].i));
        for(GModel::fiter it = GModel::current()->firstFace();
            it != GModel::current()->lastFace(); it++){
          (*it)->meshAttributes.recombine = 1;
          (*it)->meshAttributes.recombineAngle = (yyvsp[-1].i);
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
          double d;
          List_Read((yyvsp[-2].l), i, &d);
          int tag = (int)d;
          GModel::current()->getGEOInternals()->setRecombine(2, tag, (yyvsp[-1].i));
          GFace *gf = GModel::current()->getFaceByTag(tag);
          if(gf){
            gf->meshAttributes.recombine = 1;
            gf->meshAttributes.recombineAngle = (yyvsp[-1].i);
          }
        }
        List_Delete((yyvsp[-2].l));
      }
    }
#line 11223 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 354:
#line 4520 "Gmsh.y" /* yacc.c:1652  */
    {
      // recombine constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(!(yyvsp[-1].l)){
        GModel::current()->getGEOInternals()->setRecombine(3, 0, 0.);
        for(GModel::riter it = GModel::current()->firstRegion();
            it != GModel::current()->lastRegion(); it++){
          (*it)->meshAttributes.recombine3D = 1;
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
          double d;
          List_Read((yyvsp[-1].l), i, &d);
          int tag = (int)d;
          GModel::current()->getGEOInternals()->setRecombine(3, tag, 0.);
          GRegion *gr = GModel::current()->getRegionByTag(tag);
          if(gr) gr->meshAttributes.recombine3D = 1;
        }
        List_Delete((yyvsp[-1].l));
      }
    }
#line 11253 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 355:
#line 4546 "Gmsh.y" /* yacc.c:1652  */
    {
      // smoothing constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(!(yyvsp[-3].l)){
        GModel::current()->getGEOInternals()->setSmoothing(0, (int)(yyvsp[-1].d));
        for(GModel::fiter it = GModel::current()->firstFace();
            it != GModel::current()->lastFace(); it++){
          (*it)->meshAttributes.transfiniteSmoothing = (int)(yyvsp[-1].d);
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-3].l)); i++){
          double d;
          List_Read((yyvsp[-3].l), i, &d);
          int tag = (int)d;
          GModel::current()->getGEOInternals()->setSmoothing(tag, (int)(yyvsp[-1].d));
          GFace *gf = GModel::current()->getFaceByTag(tag);
          if(gf) gf->meshAttributes.transfiniteSmoothing = (int)(yyvsp[-1].d);
        }
        List_Delete((yyvsp[-3].l));
      }
    }
#line 11283 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 356:
#line 4573 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-7].l)) != List_Nbr((yyvsp[-3].l))){
        yymsg(0, "Number of master lines (%d) different from number of "
              "slaves (%d) ", List_Nbr((yyvsp[-3].l)), List_Nbr((yyvsp[-7].l)));
      }
      else{
        std::vector<double> transfo;
        if(List_Nbr((yyvsp[-1].l)) != 0) {
          if (List_Nbr((yyvsp[-1].l)) < 12){
            yymsg(0, "Affine transformation requires at least 12 entries (we have %d)",
                  List_Nbr((yyvsp[-1].l)));
          }
          else {
            transfo.resize(List_Nbr((yyvsp[-1].l)));
            for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++)
              List_Read((yyvsp[-1].l), i, &transfo[i]);
          }
        }
        for(int i = 0; i < List_Nbr((yyvsp[-7].l)); i++){
          double d_master, d_slave;
          List_Read((yyvsp[-3].l), i, &d_master);
          List_Read((yyvsp[-7].l), i, &d_slave);
          int j_master = (int)d_master;
          int j_slave  = (int)d_slave;
          addPeriodicEdge(j_slave, j_master, transfo);
        }
      }
      List_Delete((yyvsp[-7].l));
      List_Delete((yyvsp[-3].l));
    }
#line 11318 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 357:
#line 4605 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-7].l)) != List_Nbr((yyvsp[-3].l))){
        yymsg(0, "Number of master surfaces (%d) different from number of "
              "slaves (%d) ", List_Nbr((yyvsp[-3].l)), List_Nbr((yyvsp[-7].l)));
      }
      else{
        if (List_Nbr((yyvsp[-1].l)) < 12){
          // FIXME full automatic case here if List_Nbr($10) == 0)
          yymsg(0, "Affine transformation requires at least 12 entries");
        }
        else {
          std::vector<double> transfo(16,0);
          for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++)
            List_Read((yyvsp[-1].l), i, &transfo[i]);
          for(int i = 0; i < List_Nbr((yyvsp[-7].l)); i++){
            double d_master, d_slave;
            List_Read((yyvsp[-3].l), i, &d_master);
            List_Read((yyvsp[-7].l), i, &d_slave);
            addPeriodicFace(d_slave, d_master, transfo);
          }
        }
      }
      List_Delete((yyvsp[-7].l));
      List_Delete((yyvsp[-3].l));
    }
#line 11348 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 358:
#line 4632 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-14].l)) != List_Nbr((yyvsp[-10].l))){
        yymsg(0, "Number of master curves (%d) different from number of "
              "slaves (%d) ", List_Nbr((yyvsp[-10].l)), List_Nbr((yyvsp[-14].l)));
      }
      else{
        SPoint3 axis((yyvsp[-6].v)[0],(yyvsp[-6].v)[1],(yyvsp[-6].v)[2]);
        SPoint3 origin((yyvsp[-4].v)[0],(yyvsp[-4].v)[1],(yyvsp[-4].v)[2]);
        double  angle((yyvsp[-2].d));
        SPoint3 translation(0,0,0);

        std::vector<double> transfo;
        computeAffineTransformation(origin,axis,angle,translation,transfo);

        for(int i = 0; i < List_Nbr((yyvsp[-14].l)); i++){
          double d_master, d_slave;
          List_Read((yyvsp[-10].l), i, &d_master);
          List_Read((yyvsp[-14].l), i, &d_slave);
          addPeriodicEdge(d_slave,d_master,transfo);
        }
      }
      List_Delete((yyvsp[-14].l));
      List_Delete((yyvsp[-10].l));
    }
#line 11377 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 359:
#line 4658 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-14].l)) != List_Nbr((yyvsp[-10].l))){
        yymsg(0, "Number of master surfaces (%d) different from number of "
              "slaves (%d) ", List_Nbr((yyvsp[-10].l)), List_Nbr((yyvsp[-14].l)));
      }
      else{
        SPoint3 origin((yyvsp[-4].v)[0],(yyvsp[-4].v)[1],(yyvsp[-4].v)[2]);
        SPoint3 axis((yyvsp[-6].v)[0],(yyvsp[-6].v)[1],(yyvsp[-6].v)[2]);
        double  angle((yyvsp[-2].d));
        SPoint3 translation(0,0,0);

        std::vector<double> transfo;
        computeAffineTransformation(origin,axis,angle,translation,transfo);

        for(int i = 0; i < List_Nbr((yyvsp[-14].l)); i++){
          double d_master, d_slave;
          List_Read((yyvsp[-10].l), i, &d_master);
          List_Read((yyvsp[-14].l), i, &d_slave);
          addPeriodicFace(d_slave, d_master, transfo);
        }
      }
      List_Delete((yyvsp[-14].l));
      List_Delete((yyvsp[-10].l));
    }
#line 11406 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 360:
#line 4684 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-8].l)) != List_Nbr((yyvsp[-4].l))){
        yymsg(0, "Number of master curves (%d) different from number of "
              "slaves (%d) ", List_Nbr((yyvsp[-4].l)), List_Nbr((yyvsp[-8].l)));
      }
      else{
        SPoint3 origin(0,0,0);
        SPoint3 axis(0,0,0);
        double  angle(0);
        SPoint3 translation((yyvsp[-1].v)[0],(yyvsp[-1].v)[1],(yyvsp[-1].v)[2]);

        std::vector<double> transfo;
        computeAffineTransformation(origin,axis,angle,translation,transfo);

        for(int i = 0; i < List_Nbr((yyvsp[-8].l)); i++){
          double d_master, d_slave;
          List_Read((yyvsp[-4].l), i, &d_master);
          List_Read((yyvsp[-8].l), i, &d_slave);
          addPeriodicEdge(d_slave,d_master,transfo);
        }
      }
      List_Delete((yyvsp[-8].l));
      List_Delete((yyvsp[-4].l));
    }
#line 11435 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 361:
#line 4710 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-8].l)) != List_Nbr((yyvsp[-4].l))){
        yymsg(0, "Number of master surfaces (%d) different from number of "
              "slaves (%d) ", List_Nbr((yyvsp[-4].l)), List_Nbr((yyvsp[-8].l)));
      }
      else{
        SPoint3 origin(0,0,0);
        SPoint3 axis(0,0,0);
        double  angle(0);
        SPoint3 translation((yyvsp[-1].v)[0],(yyvsp[-1].v)[1],(yyvsp[-1].v)[2]);

        std::vector<double> transfo;
        computeAffineTransformation(origin,axis,angle,translation,transfo);

        for(int i = 0; i < List_Nbr((yyvsp[-8].l)); i++){
          double d_master, d_slave;
          List_Read((yyvsp[-4].l), i, &d_master);
          List_Read((yyvsp[-8].l), i, &d_slave);
          addPeriodicFace(d_slave, d_master, transfo);
        }
      }
      List_Delete((yyvsp[-8].l));
      List_Delete((yyvsp[-4].l));
    }
#line 11464 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 362:
#line 4736 "Gmsh.y" /* yacc.c:1652  */
    {
      if (List_Nbr((yyvsp[-7].l)) != List_Nbr((yyvsp[-2].l))){
        yymsg(0, "Number of master surface curves (%d) different from number of "
              "slave (%d) curves", List_Nbr((yyvsp[-2].l)), List_Nbr((yyvsp[-7].l)));
      }
      else{
        int j_master = (int)(yyvsp[-4].d);
        int j_slave = (int)(yyvsp[-9].d);
        std::map<int,int> edgeCounterParts;
        for (int i = 0; i < List_Nbr((yyvsp[-7].l)); i++){
          double ds,dm;
          List_Read((yyvsp[-7].l),i,&ds);
          List_Read((yyvsp[-2].l),i,&dm);
          edgeCounterParts[(int) ds] = (int) dm;
        }
        addPeriodicFace(j_slave, j_master, edgeCounterParts);
      }
      List_Delete((yyvsp[-7].l));
      List_Delete((yyvsp[-2].l));
    }
#line 11489 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 363:
#line 4757 "Gmsh.y" /* yacc.c:1652  */
    {
      if (((yyvsp[-4].i)==2 || (yyvsp[-4].i)==3) && (yyvsp[-9].i)<(yyvsp[-4].i) ) {
        std::vector<int> tags; ListOfDouble2Vector((yyvsp[-7].l), tags);
        addEmbedded((yyvsp[-9].i), tags, (yyvsp[-4].i), (int)(yyvsp[-2].d));
      }
      else {
        yymsg(0, "GeoEntity of dim %d In GeoEntity of dim %d not allowed", (yyvsp[-9].i), (yyvsp[-4].i));
      }
      List_Delete((yyvsp[-7].l));
    }
#line 11504 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 364:
#line 4768 "Gmsh.y" /* yacc.c:1652  */
    {
      // reverse mesh constraints are stored in GEO internals in addition to
      // GModel, as they can be copied around during GEO operations
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(!(yyvsp[-1].l)){
        GModel::current()->getGEOInternals()->setReverseMesh((yyvsp[-2].i), 0);
        switch ((yyvsp[-2].i)) {
        case 1:
          for(GModel::eiter it = GModel::current()->firstEdge();
              it != GModel::current()->lastEdge(); it++){
            (*it)->meshAttributes.reverseMesh = 1;
          }
          break;
        case 2:
          for(GModel::fiter it = GModel::current()->firstFace();
              it != GModel::current()->lastFace(); it++){
            (*it)->meshAttributes.reverseMesh = 1;
          }
          break;
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
          double d;
          List_Read((yyvsp[-1].l), i, &d);
          int num = (int)d;
          GModel::current()->getGEOInternals()->setReverseMesh((yyvsp[-2].i), num);
          switch ((yyvsp[-2].i)) {
          case 1:
            {
              GEdge *ge = GModel::current()->getEdgeByTag(num);
              if(ge) ge->meshAttributes.reverseMesh = 1;
            }
            break;
          case 2:
            {
              GFace *gf = GModel::current()->getFaceByTag(num);
              if(gf) gf->meshAttributes.reverseMesh = 1;
            }
            break;
          }
        }
        List_Delete((yyvsp[-1].l));
      }
    }
#line 11556 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 365:
#line 4816 "Gmsh.y" /* yacc.c:1652  */
    {
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(GModel::current()->getGEOInternals()->getChanged())
        GModel::current()->getGEOInternals()->synchronize(GModel::current());
      if(!(yyvsp[-1].l)){
        switch ((yyvsp[-2].i)) {
        case 0:
          for(GModel::viter it = GModel::current()->firstVertex();
              it != GModel::current()->lastVertex(); it++)
            (*it)->relocateMeshVertices();
          break;
        case 1:
          for(GModel::eiter it = GModel::current()->firstEdge();
              it != GModel::current()->lastEdge(); it++)
            (*it)->relocateMeshVertices();
          break;
        case 2:
          for(GModel::fiter it = GModel::current()->firstFace();
              it != GModel::current()->lastFace(); it++)
            (*it)->relocateMeshVertices();
          break;
        }
      }
      else{
        for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
          double d;
          List_Read((yyvsp[-1].l), i, &d);
          switch ((yyvsp[-2].i)) {
          case 0:
            {
              GVertex *gv = GModel::current()->getVertexByTag((int)d);
              if(gv) gv->relocateMeshVertices();
            }
            break;
          case 1:
            {
              GEdge *ge = GModel::current()->getEdgeByTag((int)d);
              if(ge) ge->relocateMeshVertices();
            }
            break;
          case 2:
            {
              GFace *gf = GModel::current()->getFaceByTag((int)d);
              if(gf) gf->relocateMeshVertices();
            }
            break;
          }
        }
        List_Delete((yyvsp[-1].l));
      }
    }
#line 11614 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 366:
#line 4870 "Gmsh.y" /* yacc.c:1652  */
    {
      if(GModel::current()->getOCCInternals() &&
         GModel::current()->getOCCInternals()->getChanged())
        GModel::current()->getOCCInternals()->synchronize(GModel::current());
      if(GModel::current()->getGEOInternals()->getChanged())
        GModel::current()->getGEOInternals()->synchronize(GModel::current());
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
        double d;
        List_Read((yyvsp[-1].l), i, &d);
        GRegion *gr = GModel::current()->getRegionByTag((int)d);
        if(gr) gr->setOutwardOrientationMeshConstraint();
      }
      List_Delete((yyvsp[-1].l));
    }
#line 11633 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 367:
#line 4885 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
	double dnum;
	List_Read((yyvsp[-1].l), i, &dnum);
	int num = (int)dnum;
        GModel::current()->getGEOInternals()->setDegenerated(1, num);
        GEdge *ge = GModel::current()->getEdgeByTag(num);
        if(ge) ge->setTooSmall(true);
      }
      List_Delete((yyvsp[-1].l));
    }
#line 11649 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 368:
#line 4897 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-1].l), tags);
      GModel::current()->getGEOInternals()->setCompoundMesh((yyvsp[-2].i), tags);
      List_Delete((yyvsp[-1].l));
    }
#line 11659 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 369:
#line 4908 "Gmsh.y" /* yacc.c:1652  */
    {
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals())
        GModel::current()->getOCCInternals()->removeAllDuplicates();
      else
        GModel::current()->getGEOInternals()->removeAllDuplicates();
    }
#line 11670 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 370:
#line 4915 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[-1].c), "Geometry")){
        if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals())
          GModel::current()->getOCCInternals()->removeAllDuplicates();
        else
          GModel::current()->getGEOInternals()->removeAllDuplicates();
      }
      else if(!strcmp((yyvsp[-1].c), "Mesh")){
        GModel::current()->removeDuplicateMeshVertices(CTX::instance()->geom.tolerance);
      }
      else
        yymsg(0, "Unknown coherence command");
      Free((yyvsp[-1].c));
    }
#line 11689 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 371:
#line 4930 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> tags; ListOfDouble2Vector((yyvsp[-2].l), tags);
      if(gmsh_yyfactory == "OpenCASCADE" && GModel::current()->getOCCInternals())
        GModel::current()->getOCCInternals()->mergeVertices(tags);
      else
        GModel::current()->getGEOInternals()->mergeVertices(tags);
      List_Delete((yyvsp[-2].l));
    }
#line 11702 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 372:
#line 4943 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)"Homology"; }
#line 11708 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 373:
#line 4944 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)"Cohomology"; }
#line 11714 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 374:
#line 4945 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (char*)"Betti"; }
#line 11720 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 375:
#line 4950 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> domain, subdomain, dim;
      for(int i = 0; i < 4; i++) dim.push_back(i);
      GModel::current()->addHomologyRequest((yyvsp[-1].c), domain, subdomain, dim);
    }
#line 11730 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 376:
#line 4956 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> domain, subdomain, dim;
      for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
        double d;
        List_Read((yyvsp[-2].l), i, &d);
        domain.push_back((int)d);
      }
      for(int i = 0; i < 4; i++) dim.push_back(i);
      GModel::current()->addHomologyRequest((yyvsp[-4].c), domain, subdomain, dim);
      List_Delete((yyvsp[-2].l));
    }
#line 11746 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 377:
#line 4968 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> domain, subdomain, dim;
      for(int i = 0; i < List_Nbr((yyvsp[-4].l)); i++){
        double d;
        List_Read((yyvsp[-4].l), i, &d);
        domain.push_back((int)d);
      }
      for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
        double d;
        List_Read((yyvsp[-2].l), i, &d);
        subdomain.push_back((int)d);
      }
      for(int i = 0; i < 4; i++) dim.push_back(i);
      GModel::current()->addHomologyRequest((yyvsp[-6].c), domain, subdomain, dim);
      List_Delete((yyvsp[-4].l));
      List_Delete((yyvsp[-2].l));
    }
#line 11768 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 378:
#line 4986 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<int> domain, subdomain, dim;
      for(int i = 0; i < List_Nbr((yyvsp[-4].l)); i++){
        double d;
        List_Read((yyvsp[-4].l), i, &d);
        domain.push_back((int)d);
      }
      for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
        double d;
        List_Read((yyvsp[-2].l), i, &d);
        subdomain.push_back((int)d);
      }
      for(int i = 0; i < List_Nbr((yyvsp[-7].l)); i++){
        double d;
        List_Read((yyvsp[-7].l), i, &d);
        dim.push_back((int)d);
      }
      GModel::current()->addHomologyRequest((yyvsp[-9].c), domain, subdomain, dim);
      List_Delete((yyvsp[-4].l));
      List_Delete((yyvsp[-2].l));
      List_Delete((yyvsp[-7].l));
    }
#line 11795 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 379:
#line 5013 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[0].d);           }
#line 11801 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 380:
#line 5014 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-1].d);           }
#line 11807 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 381:
#line 5015 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = -(yyvsp[0].d);          }
#line 11813 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 382:
#line 5016 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[0].d);           }
#line 11819 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 383:
#line 5017 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = !(yyvsp[0].d);          }
#line 11825 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 384:
#line 5018 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) - (yyvsp[0].d);      }
#line 11831 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 385:
#line 5019 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) + (yyvsp[0].d);      }
#line 11837 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 386:
#line 5020 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) * (yyvsp[0].d);      }
#line 11843 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 387:
#line 5022 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!(yyvsp[0].d))
	yymsg(0, "Division by zero in '%g / %g'", (yyvsp[-2].d), (yyvsp[0].d));
      else
	(yyval.d) = (yyvsp[-2].d) / (yyvsp[0].d);
    }
#line 11854 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 388:
#line 5028 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (int)(yyvsp[-2].d) | (int)(yyvsp[0].d); }
#line 11860 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 389:
#line 5029 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (int)(yyvsp[-2].d) & (int)(yyvsp[0].d); }
#line 11866 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 390:
#line 5030 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (int)(yyvsp[-2].d) % (int)(yyvsp[0].d); }
#line 11872 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 391:
#line 5031 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = pow((yyvsp[-2].d), (yyvsp[0].d));  }
#line 11878 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 392:
#line 5032 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) < (yyvsp[0].d);      }
#line 11884 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 393:
#line 5033 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) > (yyvsp[0].d);      }
#line 11890 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 394:
#line 5034 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) <= (yyvsp[0].d);     }
#line 11896 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 395:
#line 5035 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) >= (yyvsp[0].d);     }
#line 11902 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 396:
#line 5036 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) == (yyvsp[0].d);     }
#line 11908 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 397:
#line 5037 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) != (yyvsp[0].d);     }
#line 11914 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 398:
#line 5038 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) && (yyvsp[0].d);     }
#line 11920 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 399:
#line 5039 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-2].d) || (yyvsp[0].d);     }
#line 11926 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 400:
#line 5040 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = ((int)(yyvsp[-2].d) >> (int)(yyvsp[0].d)); }
#line 11932 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 401:
#line 5041 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = ((int)(yyvsp[-2].d) << (int)(yyvsp[0].d)); }
#line 11938 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 402:
#line 5042 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-4].d) ? (yyvsp[-2].d) : (yyvsp[0].d); }
#line 11944 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 403:
#line 5043 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = exp((yyvsp[-1].d));      }
#line 11950 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 404:
#line 5044 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = log((yyvsp[-1].d));      }
#line 11956 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 405:
#line 5045 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = log10((yyvsp[-1].d));    }
#line 11962 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 406:
#line 5046 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = sqrt((yyvsp[-1].d));     }
#line 11968 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 407:
#line 5047 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = sin((yyvsp[-1].d));      }
#line 11974 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 408:
#line 5048 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = asin((yyvsp[-1].d));     }
#line 11980 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 409:
#line 5049 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = cos((yyvsp[-1].d));      }
#line 11986 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 410:
#line 5050 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = acos((yyvsp[-1].d));     }
#line 11992 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 411:
#line 5051 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = tan((yyvsp[-1].d));      }
#line 11998 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 412:
#line 5052 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = atan((yyvsp[-1].d));     }
#line 12004 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 413:
#line 5053 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = atan2((yyvsp[-3].d), (yyvsp[-1].d));}
#line 12010 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 414:
#line 5054 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = sinh((yyvsp[-1].d));     }
#line 12016 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 415:
#line 5055 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = cosh((yyvsp[-1].d));     }
#line 12022 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 416:
#line 5056 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = tanh((yyvsp[-1].d));     }
#line 12028 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 417:
#line 5057 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = fabs((yyvsp[-1].d));     }
#line 12034 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 418:
#line 5058 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = std::abs((yyvsp[-1].d)); }
#line 12040 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 419:
#line 5059 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = floor((yyvsp[-1].d));    }
#line 12046 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 420:
#line 5060 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = ceil((yyvsp[-1].d));     }
#line 12052 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 421:
#line 5061 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = floor((yyvsp[-1].d) + 0.5); }
#line 12058 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 422:
#line 5062 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = fmod((yyvsp[-3].d), (yyvsp[-1].d)); }
#line 12064 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 423:
#line 5063 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = fmod((yyvsp[-3].d), (yyvsp[-1].d)); }
#line 12070 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 424:
#line 5064 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = sqrt((yyvsp[-3].d) * (yyvsp[-3].d) + (yyvsp[-1].d) * (yyvsp[-1].d)); }
#line 12076 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 425:
#line 5065 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[-1].d) * (double)rand() / (double)RAND_MAX; }
#line 12082 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 426:
#line 5074 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[0].d); }
#line 12088 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 427:
#line 5075 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = 3.141592653589793; }
#line 12094 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 428:
#line 5076 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (double)ImbricatedTest; }
#line 12100 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 429:
#line 5077 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = Msg::GetCommRank(); }
#line 12106 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 430:
#line 5078 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = Msg::GetCommSize(); }
#line 12112 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 431:
#line 5079 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = GetGmshMajorVersion(); }
#line 12118 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 432:
#line 5080 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = GetGmshMinorVersion(); }
#line 12124 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 433:
#line 5081 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = GetGmshPatchVersion(); }
#line 12130 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 434:
#line 5082 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = Cpu(); }
#line 12136 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 435:
#line 5083 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = GetMemoryUsage()/1024./1024.; }
#line 12142 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 436:
#line 5084 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = TotalRam(); }
#line 12148 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 437:
#line 5089 "Gmsh.y" /* yacc.c:1652  */
    { init_options(); }
#line 12154 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 438:
#line 5091 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<double> val(1, (yyvsp[-3].d));
      Msg::ExchangeOnelabParameter("", val, floatOptions, charOptions);
      (yyval.d) = val[0];
    }
#line 12164 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 439:
#line 5097 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[0].d); }
#line 12170 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 440:
#line 5099 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = Msg::GetOnelabNumber((yyvsp[-1].c));
      Free((yyvsp[-1].c));
    }
#line 12179 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 441:
#line 5104 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = Msg::GetOnelabNumber((yyvsp[-3].c), (yyvsp[-1].d));
      Free((yyvsp[-3].c));
    }
#line 12188 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 442:
#line 5109 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_Float((yyvsp[0].c2).char1, (yyvsp[0].c2).char2);
    }
#line 12196 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 443:
#line 5114 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_Float(nullptr, (yyvsp[-3].c), 2, (int)(yyvsp[-1].d));
    }
#line 12204 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 444:
#line 5119 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_Float(nullptr, (yyvsp[-3].c), 2, (int)(yyvsp[-1].d));
    }
#line 12212 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 445:
#line 5123 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_Float((yyvsp[-1].c2).char1, (yyvsp[-1].c2).char2, 1, 0, 0., 1);
    }
#line 12220 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 446:
#line 5127 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float((yyvsp[-3].c2).char1, (yyvsp[-3].c2).char2, (yyvsp[-1].c), 0, 0., 1);
    }
#line 12228 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 447:
#line 5131 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_Float((yyvsp[-2].c2).char1, (yyvsp[-2].c2).char2, 1, 0, (yyvsp[-1].d), 2);
    }
#line 12236 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 448:
#line 5135 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float((yyvsp[-4].c2).char1, (yyvsp[-4].c2).char2, (yyvsp[-2].c), 0, (yyvsp[-1].d), 2);
    }
#line 12244 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 449:
#line 5139 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_Float((yyvsp[-5].c2).char1, (yyvsp[-5].c2).char2, 2, (int)(yyvsp[-3].d), (yyvsp[-1].d), 2);
    }
#line 12252 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 450:
#line 5143 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float((yyvsp[-7].c2).char1, (yyvsp[-7].c2).char2, (yyvsp[-5].c), (int)(yyvsp[-3].d), (yyvsp[-1].d), 2);
    }
#line 12260 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 451:
#line 5147 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
      (yyval.d) = !StatFile(tmp);
      Free((yyvsp[-1].c));
    }
#line 12270 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 452:
#line 5153 "Gmsh.y" /* yacc.c:1652  */
    {
      if(gmsh_yysymbols.count((yyvsp[-2].c))){
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-2].c)]);
	(yyval.d) = s.value.size();
      }
      else if(gmsh_yystringsymbols.count((yyvsp[-2].c))){
	(yyval.d) = gmsh_yystringsymbols[(yyvsp[-2].c)].size();
      }
      else{
        yymsg(0, "Unknown variable '%s'", (yyvsp[-2].c));
	(yyval.d) = 0.;
      }
      Free((yyvsp[-2].c));
    }
#line 12289 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 453:
#line 5169 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float_getDim((yyvsp[-4].c2).char1, (yyvsp[-4].c2).char2, (yyvsp[-2].c));
    }
#line 12297 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 454:
#line 5174 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string struct_namespace((yyvsp[-1].c));
      (yyval.d) = (double)gmsh_yynamespaces[struct_namespace].size();
      Free((yyvsp[-1].c));
    }
#line 12307 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 455:
#line 5180 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string struct_namespace(std::string(""));
      (yyval.d) = (double)gmsh_yynamespaces[struct_namespace].size();
    }
#line 12316 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 456:
#line 5186 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!gmsh_yysymbols.count((yyvsp[-1].c))){
	yymsg(0, "Unknown variable '%s'", (yyvsp[-1].c));
	(yyval.d) = 0.;
      }
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-1].c)]);
        if(s.value.empty()){
          yymsg(0, "Uninitialized variable '%s'", (yyvsp[-1].c));
          (yyval.d) = 0.;
        }
        else{
          (yyval.d) = s.value[0];
          s.value[0] += (yyvsp[0].i);
        }
      }
      Free((yyvsp[-1].c));
    }
#line 12339 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 457:
#line 5205 "Gmsh.y" /* yacc.c:1652  */
    {
      int index = (int)(yyvsp[-2].d);
      if(!gmsh_yysymbols.count((yyvsp[-4].c))){
	yymsg(0, "Unknown variable '%s'", (yyvsp[-4].c));
	(yyval.d) = 0.;
      }
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-4].c)]);
        if((int)s.value.size() < index + 1){
          yymsg(0, "Uninitialized variable '%s[%d]'", (yyvsp[-4].c), index);
          (yyval.d) = 0.;
        }
        else{
          (yyval.d) = s.value[index];
          s.value[index] += (yyvsp[0].i);
        }
      }
      Free((yyvsp[-4].c));
    }
#line 12363 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 458:
#line 5226 "Gmsh.y" /* yacc.c:1652  */
    {
      int index = (int)(yyvsp[-2].d);
      if(!gmsh_yysymbols.count((yyvsp[-4].c))){
	yymsg(0, "Unknown variable '%s'", (yyvsp[-4].c));
	(yyval.d) = 0.;
      }
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-4].c)]);
        if((int)s.value.size() < index + 1){
          yymsg(0, "Uninitialized variable '%s[%d]'", (yyvsp[-4].c), index);
          (yyval.d) = 0.;
        }
        else{
          (yyval.d) = s.value[index];
          s.value[index] += (yyvsp[0].i);
        }
      }
      Free((yyvsp[-4].c));
    }
#line 12387 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 459:
#line 5259 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float(nullptr, (yyvsp[-2].c), (yyvsp[0].c));
    }
#line 12395 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 460:
#line 5263 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float((yyvsp[-4].c), (yyvsp[-2].c), (yyvsp[0].c));
    }
#line 12403 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 461:
#line 5268 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float(nullptr, (yyvsp[-5].c), (yyvsp[-3].c), (int)(yyvsp[-1].d));
    }
#line 12411 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 462:
#line 5272 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float((yyvsp[-7].c), (yyvsp[-5].c), (yyvsp[-3].c), (int)(yyvsp[-1].d));
    }
#line 12419 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 463:
#line 5276 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float(nullptr, (yyvsp[-5].c), (yyvsp[-3].c), (int)(yyvsp[-1].d));
    }
#line 12427 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 464:
#line 5280 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = treat_Struct_FullName_dot_tSTRING_Float((yyvsp[-7].c), (yyvsp[-5].c), (yyvsp[-3].c), (int)(yyvsp[-1].d));
    }
#line 12435 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 465:
#line 5285 "Gmsh.y" /* yacc.c:1652  */
    {
      NumberOption(GMSH_GET, (yyvsp[-5].c), (int)(yyvsp[-3].d), (yyvsp[0].c), (yyval.d));
      Free((yyvsp[-5].c)); Free((yyvsp[0].c));
    }
#line 12444 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 466:
#line 5290 "Gmsh.y" /* yacc.c:1652  */
    {
      double d = 0.;
      if(NumberOption(GMSH_GET, (yyvsp[-3].c), 0, (yyvsp[-1].c), d)){
	d += (yyvsp[0].i);
	NumberOption(GMSH_SET|GMSH_GUI, (yyvsp[-3].c), 0, (yyvsp[-1].c), d);
	(yyval.d) = d;
      }
      Free((yyvsp[-3].c)); Free((yyvsp[-1].c));
    }
#line 12458 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 467:
#line 5300 "Gmsh.y" /* yacc.c:1652  */
    {
      double d = 0.;
      if(NumberOption(GMSH_GET, (yyvsp[-6].c), (int)(yyvsp[-4].d), (yyvsp[-1].c), d)){
	d += (yyvsp[0].i);
	NumberOption(GMSH_SET|GMSH_GUI, (yyvsp[-6].c), (int)(yyvsp[-4].d), (yyvsp[-1].c), d);
	(yyval.d) = d;
      }
      Free((yyvsp[-6].c)); Free((yyvsp[-1].c));
    }
#line 12472 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 468:
#line 5310 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = Msg::GetValue((yyvsp[-3].c), (yyvsp[-1].d));
      Free((yyvsp[-3].c));
    }
#line 12481 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 469:
#line 5315 "Gmsh.y" /* yacc.c:1652  */
    {
      int matches = 0;
      for(int i = 0; i < List_Nbr((yyvsp[-3].l)); i++){
        double d;
        List_Read((yyvsp[-3].l), i, &d);
        matches += List_Search((yyvsp[-1].l), &d, fcmp_double);
      }
      (yyval.d) = matches;
      Free((yyvsp[-3].l)); Free((yyvsp[-1].l));
    }
#line 12496 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 470:
#line 5326 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string s((yyvsp[-3].c)), substr((yyvsp[-1].c));
      if(s.find(substr) != std::string::npos)
        (yyval.d) = 1.;
      else
        (yyval.d) = 0.;
      Free((yyvsp[-3].c)); Free((yyvsp[-1].c));
    }
#line 12509 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 471:
#line 5335 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = strlen((yyvsp[-1].c));
      Free((yyvsp[-1].c));
    }
#line 12518 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 472:
#line 5340 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.d) = strcmp((yyvsp[-3].c), (yyvsp[-1].c));
      Free((yyvsp[-3].c)); Free((yyvsp[-1].c));
    }
#line 12527 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 473:
#line 5345 "Gmsh.y" /* yacc.c:1652  */
    {
      int align = 0, font = 0, fontsize = CTX::instance()->glFontSize;
      if(List_Nbr((yyvsp[-1].l)) % 2){
        yymsg(0, "Number of text attributes should be even");
      }
      else{
        for(int i = 0 ; i < List_Nbr((yyvsp[-1].l)); i += 2){
          char *s1, *s2; List_Read((yyvsp[-1].l), i, &s1); List_Read((yyvsp[-1].l), i + 1, &s2);
          std::string key(s1), val(s2);
          Free(s1); Free(s2);
#if defined(HAVE_OPENGL)
          if(key == "Font")
            font = drawContext::global()->getFontIndex(val.c_str());
          else if(key == "FontSize")
            fontsize = atoi(val.c_str());
          else if(key == "Align")
            align = drawContext::global()->getFontAlign(val.c_str());
#endif
        }
      }
      List_Delete((yyvsp[-1].l));
      (yyval.d) = (double)((align<<16)|(font<<8)|(fontsize));
    }
#line 12555 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 474:
#line 5372 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = 0.; }
#line 12561 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 475:
#line 5374 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.d) = (yyvsp[0].d);}
#line 12567 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 476:
#line 5379 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = nullptr; }
#line 12573 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 477:
#line 5381 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (yyvsp[0].c);}
#line 12579 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 478:
#line 5386 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string struct_namespace((yyvsp[-1].c2).char1? (yyvsp[-1].c2).char1 : std::string("")),
        struct_name((yyvsp[-1].c2).char2);
      init_options
        (gmsh_yynamespaces.getMember_ValMax(struct_namespace, struct_name));
    }
#line 12590 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 479:
#line 5393 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string struct_namespace((yyvsp[-5].c2).char1? (yyvsp[-5].c2).char1 : std::string("")),
        struct_name((yyvsp[-5].c2).char2);
      Free((yyvsp[-5].c2).char1); Free((yyvsp[-5].c2).char2);
      int tag_out;
      if (gmsh_yynamespaces.defStruct(struct_namespace, struct_name,
                                      floatOptions, charOptions,
                                      tag_out, member_ValMax, (yyvsp[-4].i)))
        yymsg(0, "Redefinition of Struct '%s::%s'",
              struct_namespace.c_str(), struct_name.c_str());
      (yyval.d) = (double)tag_out;
    }
#line 12607 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 480:
#line 5409 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c2).char1 = nullptr; (yyval.c2).char2 = (yyvsp[0].c); }
#line 12613 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 481:
#line 5411 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c2).char1 = (yyvsp[-2].c); (yyval.c2).char2 = (yyvsp[0].c); }
#line 12619 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 482:
#line 5416 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (yyvsp[0].c); flag_tSTRING_alloc = 1; }
#line 12625 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 483:
#line 5425 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 99; }
#line 12631 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 484:
#line 5427 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = (int)(yyvsp[0].d); }
#line 12637 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 485:
#line 5432 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = 0; }
#line 12643 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 486:
#line 5434 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.i) = (yyvsp[-1].i); }
#line 12649 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 487:
#line 5439 "Gmsh.y" /* yacc.c:1652  */
    {
      memcpy((yyval.v), (yyvsp[0].v), 5*sizeof(double));
    }
#line 12657 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 488:
#line 5443 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < 5; i++) (yyval.v)[i] = -(yyvsp[0].v)[i];
    }
#line 12665 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 489:
#line 5447 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < 5; i++) (yyval.v)[i] = (yyvsp[0].v)[i];
    }
#line 12673 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 490:
#line 5451 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < 5; i++) (yyval.v)[i] = (yyvsp[-2].v)[i] - (yyvsp[0].v)[i];
    }
#line 12681 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 491:
#line 5455 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < 5; i++) (yyval.v)[i] = (yyvsp[-2].v)[i] + (yyvsp[0].v)[i];
    }
#line 12689 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 492:
#line 5462 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.v)[0] = (yyvsp[-9].d);  (yyval.v)[1] = (yyvsp[-7].d);  (yyval.v)[2] = (yyvsp[-5].d);  (yyval.v)[3] = (yyvsp[-3].d); (yyval.v)[4] = (yyvsp[-1].d);
    }
#line 12697 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 493:
#line 5466 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.v)[0] = (yyvsp[-7].d);  (yyval.v)[1] = (yyvsp[-5].d);  (yyval.v)[2] = (yyvsp[-3].d);  (yyval.v)[3] = (yyvsp[-1].d); (yyval.v)[4] = 1.0;
    }
#line 12705 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 494:
#line 5470 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.v)[0] = (yyvsp[-5].d);  (yyval.v)[1] = (yyvsp[-3].d);  (yyval.v)[2] = (yyvsp[-1].d);  (yyval.v)[3] = 0.0; (yyval.v)[4] = 1.0;
    }
#line 12713 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 495:
#line 5474 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.v)[0] = (yyvsp[-5].d);  (yyval.v)[1] = (yyvsp[-3].d);  (yyval.v)[2] = (yyvsp[-1].d);  (yyval.v)[3] = 0.0; (yyval.v)[4] = 1.0;
    }
#line 12721 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 496:
#line 5481 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(List_T*));
      List_Add((yyval.l), &((yyvsp[0].l)));
    }
#line 12730 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 497:
#line 5486 "Gmsh.y" /* yacc.c:1652  */
    {
      List_Add((yyval.l), &((yyvsp[0].l)));
    }
#line 12738 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 498:
#line 5493 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(double));
      List_Add((yyval.l), &((yyvsp[0].d)));
    }
#line 12747 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 499:
#line 5498 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[0].l);
    }
#line 12755 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 500:
#line 5502 "Gmsh.y" /* yacc.c:1652  */
    {
      // creates an empty list
      (yyval.l) = List_Create(2, 1, sizeof(double));
    }
#line 12764 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 501:
#line 5507 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
    }
#line 12772 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 502:
#line 5511 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
      for(int i = 0; i < List_Nbr((yyval.l)); i++){
	double *pd = (double*)List_Pointer((yyval.l), i);
	(*pd) = - (*pd);
      }
    }
#line 12784 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 503:
#line 5519 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
      for(int i = 0; i < List_Nbr((yyval.l)); i++){
	double *pd = (double*)List_Pointer((yyval.l), i);
	(*pd) *= (yyvsp[-4].d);
      }
    }
#line 12796 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 504:
#line 5530 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[0].l);
    }
#line 12804 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 505:
#line 5534 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = 0;
    }
#line 12812 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 506:
#line 5538 "Gmsh.y" /* yacc.c:1652  */
    {
      if(!strcmp((yyvsp[0].c), "*") || !strcmp((yyvsp[0].c), "all")){
        (yyval.l) = 0;
      }
      else{
        yyerror("Unknown special string for list replacement");
        (yyval.l) = List_Create(2, 1, sizeof(double));
      }
      Free((yyvsp[0].c));
    }
#line 12827 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 507:
#line 5552 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[0].l);
      for(int i = 0; i < List_Nbr((yyval.l)); i++){
	double *pd = (double*)List_Pointer((yyval.l), i);
	(*pd) = - (*pd);
      }
    }
#line 12839 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 508:
#line 5560 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[0].l);
      for(int i = 0; i < List_Nbr((yyval.l)); i++){
	double *pd = (double*)List_Pointer((yyval.l), i);
	(*pd) *= (yyvsp[-2].d);
      }
    }
#line 12851 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 509:
#line 5568 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(double));
      for(double d = (yyvsp[-2].d); ((yyvsp[-2].d) < (yyvsp[0].d)) ? (d <= (yyvsp[0].d)) : (d >= (yyvsp[0].d));
          ((yyvsp[-2].d) < (yyvsp[0].d)) ? (d += 1.) : (d -= 1.))
	List_Add((yyval.l), &d);
    }
#line 12862 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 510:
#line 5575 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(double));
      if(!(yyvsp[0].d)){  //|| ($1 < $3 && $5 < 0) || ($1 > $3 && $5 > 0)
        yymsg(0, "Wrong increment in '%g:%g:%g'", (yyvsp[-4].d), (yyvsp[-2].d), (yyvsp[0].d));
      }
      else
	for(double d = (yyvsp[-4].d); ((yyvsp[0].d) > 0) ? (d <= (yyvsp[-2].d)) : (d >= (yyvsp[-2].d)); d += (yyvsp[0].d))
	  List_Add((yyval.l), &d);
   }
#line 12876 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 511:
#line 5585 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(3, 1, sizeof(double));
      int tag = (int)(yyvsp[-1].d);
      double x = 0., y = 0., z = 0.;
      bool found = GModel::current()->getGEOInternals()->getVertex(tag, x, y, z);
      if(!found && GModel::current()->getOCCInternals())
        found = GModel::current()->getOCCInternals()->getVertex(tag, x, y, z);
      if(!found){
        GVertex *gv = GModel::current()->getVertexByTag(tag);
        if(gv){
          x = gv->x();
          y = gv->y();
          z = gv->z();
        }
        else{
          yymsg(0, "Unknown model point with tag %d", tag);
        }
      }
      List_Add((yyval.l), &x);
      List_Add((yyval.l), &y);
      List_Add((yyval.l), &z);
    }
#line 12903 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 512:
#line 5608 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getAllElementaryTags(0, (yyval.l));
    }
#line 12912 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 513:
#line 5613 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getAllElementaryTags(0, (yyval.l));
      Free((yyvsp[0].c));
    }
#line 12922 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 514:
#line 5619 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getAllElementaryTags((yyvsp[-3].i), (yyval.l));
    }
#line 12931 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 515:
#line 5624 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getAllElementaryTags((yyvsp[-1].i), (yyval.l));
      Free((yyvsp[0].c));
    }
#line 12941 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 516:
#line 5630 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      if(!(yyvsp[0].l)){
        getAllPhysicalTags((yyvsp[-1].i), (yyval.l));
      }
      else{
        getElementaryTagsForPhysicalGroups((yyvsp[-1].i), (yyvsp[0].l), (yyval.l));
        List_Delete((yyvsp[0].l));
      }
    }
#line 12956 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 517:
#line 5641 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getParentTags((yyvsp[-1].i), (yyvsp[0].l), (yyval.l));
      List_Delete((yyvsp[0].l));
    }
#line 12966 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 518:
#line 5648 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getElementaryTagsInBoundingBox((yyvsp[-15].i), (yyvsp[-11].d), (yyvsp[-9].d), (yyvsp[-7].d), (yyvsp[-5].d), (yyvsp[-3].d), (yyvsp[-1].d), (yyval.l));
    }
#line 12975 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 519:
#line 5653 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(10, 10, sizeof(double));
      getBoundingBox((yyvsp[-3].i), (yyvsp[-1].l), (yyval.l));
      List_Delete((yyvsp[-1].l));
    }
#line 12985 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 520:
#line 5659 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(List_Nbr((yyvsp[0].l)), 1, sizeof(double));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
	Shape *s = (Shape*) List_Pointer((yyvsp[0].l), i);
	double d = s->Num;
	List_Add((yyval.l), &d);
      }
      List_Delete((yyvsp[0].l));
    }
#line 12999 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 521:
#line 5669 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(List_Nbr((yyvsp[0].l)), 1, sizeof(double));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
	Shape *s = (Shape*) List_Pointer((yyvsp[0].l), i);
	double d = s->Num;
	List_Add((yyval.l), &d);
      }
      List_Delete((yyvsp[0].l));
    }
#line 13013 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 522:
#line 5679 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(List_Nbr((yyvsp[0].l)), 1, sizeof(double));
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
	Shape *s = (Shape*) List_Pointer((yyvsp[0].l), i);
	double d = s->Num;
	List_Add((yyval.l), &d);
      }
      List_Delete((yyvsp[0].l));
    }
#line 13027 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 523:
#line 5689 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20, 20, sizeof(double));
      if(!gmsh_yysymbols.count((yyvsp[-2].c)))
	yymsg(0, "Unknown variable '%s'", (yyvsp[-2].c));
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-2].c)]);
	for(std::size_t i = 0; i < s.value.size(); i++)
	  List_Add((yyval.l), &s.value[i]);
      }
      Free((yyvsp[-2].c));
    }
#line 13043 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 524:
#line 5701 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = treat_Struct_FullName_dot_tSTRING_ListOfFloat(nullptr, (yyvsp[-4].c), (yyvsp[-2].c));
    }
#line 13051 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 525:
#line 5705 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = treat_Struct_FullName_dot_tSTRING_ListOfFloat((yyvsp[-6].c), (yyvsp[-4].c), (yyvsp[-2].c));
    }
#line 13059 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 526:
#line 5710 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(double));
      if(!gmsh_yysymbols.count((yyvsp[-1].c)))
	yymsg(0, "Unknown variable '%s'", (yyvsp[-1].c));
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-1].c)]);
	for(std::size_t i = 0; i < s.value.size(); i++)
	  List_Add((yyval.l), &s.value[i]);
      }
      Free((yyvsp[-1].c));
    }
#line 13075 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 527:
#line 5722 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
    }
#line 13083 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 528:
#line 5726 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
    }
#line 13091 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 529:
#line 5730 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-2].l);
    }
#line 13099 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 530:
#line 5734 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(double));
      if(!gmsh_yysymbols.count((yyvsp[-5].c)))
	yymsg(0, "Unknown variable '%s'", (yyvsp[-5].c));
      else{
        gmsh_yysymbol &s(gmsh_yysymbols[(yyvsp[-5].c)]);
	for(int i = 0; i < List_Nbr((yyvsp[-2].l)); i++){
	  int index = (int)(*(double*)List_Pointer_Fast((yyvsp[-2].l), i));
	  if((int)s.value.size() < index + 1)
	    yymsg(0, "Uninitialized variable '%s[%d]'", (yyvsp[-5].c), index);
	  else
	    List_Add((yyval.l), &s.value[index]);
	}
      }
      Free((yyvsp[-5].c));
      List_Delete((yyvsp[-2].l));
    }
#line 13121 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 531:
#line 5752 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20,20,sizeof(double));
      for(int i = 0; i < (int)(yyvsp[-1].d); i++) {
	double d = (yyvsp[-5].d) + ((yyvsp[-3].d)-(yyvsp[-5].d))*(double)i/((yyvsp[-1].d)-1);
	List_Add((yyval.l), &d);
      }
    }
#line 13133 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 532:
#line 5760 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20,20,sizeof(double));
      for(int i = 0; i < (int)(yyvsp[-1].d); i++) {
	double d = pow(10,(yyvsp[-5].d) + ((yyvsp[-3].d)-(yyvsp[-5].d))*(double)i/((yyvsp[-1].d)-1));
	List_Add((yyval.l), &d);
      }
    }
#line 13145 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 533:
#line 5768 "Gmsh.y" /* yacc.c:1652  */
    {
      Msg::Barrier();
      FILE *File;
      (yyval.l) = List_Create(100, 100, sizeof(double));
      std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
      if(!(File = Fopen(tmp.c_str(), "rb"))){
        yymsg(0, "Could not open file '%s'", (yyvsp[-1].c));
      }
      else{
	double d;
	while(!feof(File)){
          int ret = fscanf(File, "%lf", &d);
	  if(ret == 1){
	    List_Add((yyval.l), &d);
          }
          else if(ret == EOF){
            break;
          }
          else{
            char dummy[1024];
            fscanf(File, "%s", dummy);
            yymsg(0, "Ignoring '%s' in file '%s'", dummy, (yyvsp[-1].c));
          }
        }
	fclose(File);
      }
      Free((yyvsp[-1].c));
    }
#line 13178 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 534:
#line 5797 "Gmsh.y" /* yacc.c:1652  */
    {
      double x0 = (yyvsp[-11].d), x1 = (yyvsp[-9].d), y0 = (yyvsp[-7].d), y1 = (yyvsp[-5].d), ys = (yyvsp[-3].d);
      int N = (int)(yyvsp[-1].d);
      std::vector<double> y(N);
      if(!catenary(x0, x1, y0, y1, ys, N, &y[0]))
        Msg::Warning("Catenary did not converge, using linear interpolation");
      (yyval.l) = List_Create(N,10,sizeof(double));
      for(int i = 0; i < N; i++) List_Add((yyval.l), &y[i]);
    }
#line 13192 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 535:
#line 5807 "Gmsh.y" /* yacc.c:1652  */
    {
      std::vector<double> tmp;
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
        double d; List_Read((yyvsp[-1].l), i, &d);
        tmp.push_back(d);
      }
      std::sort(tmp.begin(), tmp.end());
      std::vector<double>::iterator last = std::unique(tmp.begin(), tmp.end());
      tmp.erase(last, tmp.end());
      (yyval.l) = (yyvsp[-1].l);
      List_Reset((yyval.l));
      for(std::size_t i = 0; i < tmp.size(); i++){
        List_Add((yyval.l), &tmp[i]);
      }
    }
#line 13212 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 536:
#line 5823 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
        double *d = (double*)List_Pointer((yyvsp[-1].l), i);
        *d = std::abs(*d);
      }
      (yyval.l) = (yyvsp[-1].l);
    }
#line 13224 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 537:
#line 5834 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(2, 1, sizeof(double));
      List_Add((yyval.l), &((yyvsp[0].d)));
    }
#line 13233 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 538:
#line 5839 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[0].l);
    }
#line 13241 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 539:
#line 5843 "Gmsh.y" /* yacc.c:1652  */
    {
      List_Add((yyval.l), &((yyvsp[0].d)));
    }
#line 13249 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 540:
#line 5847 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
	double d;
	List_Read((yyvsp[0].l), i, &d);
	List_Add((yyval.l), &d);
      }
      List_Delete((yyvsp[0].l));
    }
#line 13262 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 541:
#line 5859 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.u) = CTX::instance()->packColor((int)(yyvsp[-7].d), (int)(yyvsp[-5].d), (int)(yyvsp[-3].d), (int)(yyvsp[-1].d));
    }
#line 13270 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 542:
#line 5863 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.u) = CTX::instance()->packColor((int)(yyvsp[-5].d), (int)(yyvsp[-3].d), (int)(yyvsp[-1].d), 255);
    }
#line 13278 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 543:
#line 5875 "Gmsh.y" /* yacc.c:1652  */
    {
      int flag = 0;
      if(gmsh_yystringsymbols.count((yyvsp[0].c))){
        if(gmsh_yystringsymbols[(yyvsp[0].c)].size()){
          (yyval.u) = GetColorForString(-1, gmsh_yystringsymbols[(yyvsp[0].c)][0].c_str(), &flag);
        }
        else{
          yymsg(0, "Unknown color '%s'", (yyvsp[0].c));
          (yyval.u) = 0;
        }
      }
      else
        (yyval.u) = GetColorForString(-1, (yyvsp[0].c), &flag);
      if(flag) yymsg(0, "Unknown color '%s'", (yyvsp[0].c));
      Free((yyvsp[0].c));
    }
#line 13299 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 544:
#line 5892 "Gmsh.y" /* yacc.c:1652  */
    {
      unsigned int val = 0;
      ColorOption(GMSH_GET, (yyvsp[-4].c), 0, (yyvsp[0].c), val);
      (yyval.u) = val;
      Free((yyvsp[-4].c)); Free((yyvsp[0].c));
    }
#line 13310 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 545:
#line 5902 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = (yyvsp[-1].l);
    }
#line 13318 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 546:
#line 5906 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(256, 10, sizeof(unsigned int));
      GmshColorTable *ct = GetColorTable((int)(yyvsp[-3].d));
      if(!ct)
	yymsg(0, "View[%d] does not exist", (int)(yyvsp[-3].d));
      else{
	for(int i = 0; i < ct->size; i++)
	  List_Add((yyval.l), &ct->table[i]);
      }
      Free((yyvsp[-5].c));
    }
#line 13334 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 547:
#line 5921 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(256, 10, sizeof(unsigned int));
      List_Add((yyval.l), &((yyvsp[0].u)));
    }
#line 13343 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 548:
#line 5926 "Gmsh.y" /* yacc.c:1652  */
    {
      List_Add((yyval.l), &((yyvsp[0].u)));
    }
#line 13351 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 549:
#line 5933 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (yyvsp[0].c);
    }
#line 13359 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 550:
#line 5937 "Gmsh.y" /* yacc.c:1652  */
    {
      // No need to extend to Struct_FullName (a Tag is not a String)
      (yyval.c) = treat_Struct_FullName_String(nullptr, (yyvsp[0].c));
    }
#line 13368 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 551:
#line 5942 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string val;
      int j = (int)(yyvsp[-1].d);
      if(!gmsh_yystringsymbols.count((yyvsp[-3].c)))
        yymsg(0, "Unknown string variable '%s'", (yyvsp[-3].c));
      else if(j >= 0 && j < (int)gmsh_yystringsymbols[(yyvsp[-3].c)].size())
        val = gmsh_yystringsymbols[(yyvsp[-3].c)][j];
      else
        yymsg(0, "Index %d out of range", j);
      (yyval.c) = (char *)Malloc((val.size() + 1) * sizeof(char));
      strcpy((yyval.c), val.c_str());
      Free((yyvsp[-3].c));
    }
#line 13386 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 552:
#line 5956 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string val;
      int j = (int)(yyvsp[-1].d);
      if(!gmsh_yystringsymbols.count((yyvsp[-3].c)))
        yymsg(0, "Unknown string variable '%s'", (yyvsp[-3].c));
      else if(j >= 0 && j < (int)gmsh_yystringsymbols[(yyvsp[-3].c)].size())
        val = gmsh_yystringsymbols[(yyvsp[-3].c)][j];
      else
        yymsg(0, "Index %d out of range", j);
      (yyval.c) = (char *)Malloc((val.size() + 1) * sizeof(char));
      strcpy((yyval.c), val.c_str());
      Free((yyvsp[-3].c));
    }
#line 13404 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 553:
#line 5972 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = treat_Struct_FullName_dot_tSTRING_String(nullptr, (yyvsp[-2].c), (yyvsp[0].c));
    }
#line 13412 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 554:
#line 5976 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = treat_Struct_FullName_dot_tSTRING_String((yyvsp[-4].c), (yyvsp[-2].c), (yyvsp[0].c));
    }
#line 13420 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 555:
#line 5980 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = treat_Struct_FullName_dot_tSTRING_String(nullptr, (yyvsp[-5].c), (yyvsp[-3].c), (int)(yyvsp[-1].d));
    }
#line 13428 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 556:
#line 5984 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = treat_Struct_FullName_dot_tSTRING_String((yyvsp[-7].c), (yyvsp[-5].c), (yyvsp[-3].c), (int)(yyvsp[-1].d));
    }
#line 13436 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 557:
#line 5988 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string out;
      StringOption(GMSH_GET, (yyvsp[-5].c), (int)(yyvsp[-3].d), (yyvsp[0].c), out);
      (yyval.c) = (char*)Malloc((out.size() + 1) * sizeof(char));
      strcpy((yyval.c), out.c_str());
      Free((yyvsp[-5].c)); Free((yyvsp[0].c));
    }
#line 13448 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 558:
#line 5996 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string name = GModel::current()->getElementaryName((yyvsp[-3].i), (int)(yyvsp[-1].d));
      (yyval.c) = (char*)Malloc((name.size() + 1) * sizeof(char));
      strcpy((yyval.c), name.c_str());
    }
#line 13458 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 559:
#line 6002 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string name = GModel::current()->getPhysicalName((yyvsp[-3].i), (int)(yyvsp[-1].d));
      (yyval.c) = (char*)Malloc((name.size() + 1) * sizeof(char));
      strcpy((yyval.c), name.c_str());
    }
#line 13468 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 560:
#line 6011 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (yyvsp[0].c);
    }
#line 13476 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 561:
#line 6015 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (yyvsp[-1].c);
    }
#line 13484 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 562:
#line 6019 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (char *)Malloc(32 * sizeof(char));
      time_t now;
      time(&now);
      strcpy((yyval.c), ctime(&now));
      (yyval.c)[strlen((yyval.c)) - 1] = '\0';
    }
#line 13496 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 563:
#line 6027 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string exe = Msg::GetExecutableName();
      (yyval.c) = (char *)Malloc(exe.size() + 1);
      strcpy((yyval.c), exe.c_str());
    }
#line 13506 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 564:
#line 6033 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string action = Msg::GetOnelabAction();
      (yyval.c) = (char *)Malloc(action.size() + 1);
      strcpy((yyval.c), action.c_str());
    }
#line 13516 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 565:
#line 6039 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = strsave((char*)"Gmsh");
    }
#line 13524 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 566:
#line 6043 "Gmsh.y" /* yacc.c:1652  */
    {
      const char *env = GetEnvironmentVar((yyvsp[-1].c));
      if(!env) env = "";
      (yyval.c) = (char *)Malloc((sizeof(env) + 1) * sizeof(char));
      strcpy((yyval.c), env);
      Free((yyvsp[-1].c));
    }
#line 13536 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 567:
#line 6051 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string s = Msg::GetString((yyvsp[-3].c), (yyvsp[-1].c));
      (yyval.c) = (char *)Malloc((s.size() + 1) * sizeof(char));
      strcpy((yyval.c), s.c_str());
      Free((yyvsp[-3].c));
      Free((yyvsp[-1].c));
    }
#line 13548 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 568:
#line 6059 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string s = Msg::GetOnelabString((yyvsp[-1].c));
      (yyval.c) = (char *)Malloc((s.size() + 1) * sizeof(char));
      strcpy((yyval.c), s.c_str());
      Free((yyvsp[-1].c));
    }
#line 13559 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 569:
#line 6066 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string s = Msg::GetOnelabString((yyvsp[-3].c), (yyvsp[-1].c));
      (yyval.c) = (char *)Malloc((s.size() + 1) * sizeof(char));
      strcpy((yyval.c), s.c_str());
      Free((yyvsp[-3].c));
      Free((yyvsp[-1].c));
    }
#line 13571 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 570:
#line 6075 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = treat_Struct_FullName_String(nullptr, (yyvsp[-2].c2).char2, 1, 0, (yyvsp[-1].c), 2);
    }
#line 13579 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 571:
#line 6079 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = treat_Struct_FullName_dot_tSTRING_String((yyvsp[-4].c2).char1, (yyvsp[-4].c2).char2, (yyvsp[-2].c), 0, (yyvsp[-1].c), 2);
    }
#line 13587 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 572:
#line 6083 "Gmsh.y" /* yacc.c:1652  */
    {
      int size = 1;
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++)
        size += strlen(*(char**)List_Pointer((yyvsp[-1].l), i)) + 1;
      (yyval.c) = (char*)Malloc(size * sizeof(char));
      (yyval.c)[0] = '\0';
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
        char *s;
        List_Read((yyvsp[-1].l), i, &s);
        strcat((yyval.c), s);
        Free(s);
      }
      List_Delete((yyvsp[-1].l));
    }
#line 13606 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 573:
#line 6098 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (char *)Malloc((strlen((yyvsp[-1].c)) + 1) * sizeof(char));
      int i;
      for(i = strlen((yyvsp[-1].c)) - 1; i >= 0; i--){
	if((yyvsp[-1].c)[i] == '.'){
	  strncpy((yyval.c), (yyvsp[-1].c), i);
	  (yyval.c)[i]='\0';
	  break;
	}
      }
      if(i <= 0) strcpy((yyval.c), (yyvsp[-1].c));
      Free((yyvsp[-1].c));
    }
#line 13624 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 574:
#line 6112 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (char *)Malloc((strlen((yyvsp[-1].c)) + 1) * sizeof(char));
      int i;
      for(i = strlen((yyvsp[-1].c)) - 1; i >= 0; i--){
	if((yyvsp[-1].c)[i] == '/' || (yyvsp[-1].c)[i] == '\\')
	  break;
      }
      if(i <= 0)
	strcpy((yyval.c), (yyvsp[-1].c));
      else
	strcpy((yyval.c), &(yyvsp[-1].c)[i+1]);
      Free((yyvsp[-1].c));
    }
#line 13642 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 575:
#line 6126 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string input = (yyvsp[-5].c);
      std::string substr_old = (yyvsp[-3].c);
      std::string substr_new = (yyvsp[-1].c);
      std::string ret = ReplaceSubString(substr_old, substr_new, input);
      (yyval.c) = (char *)Malloc((ret.size() + 1) * sizeof(char));
      strcpy((yyval.c), ret.c_str());
      Free((yyvsp[-5].c));
      Free((yyvsp[-3].c));
      Free((yyvsp[-1].c));
    }
#line 13658 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 576:
#line 6138 "Gmsh.y" /* yacc.c:1652  */
    {
      int size = 1;
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++)
        size += strlen(*(char**)List_Pointer((yyvsp[-1].l), i)) + 1;
      (yyval.c) = (char*)Malloc(size * sizeof(char));
      (yyval.c)[0] = '\0';
      for(int i = 0; i < List_Nbr((yyvsp[-1].l)); i++){
        char *s;
        List_Read((yyvsp[-1].l), i, &s);
        strcat((yyval.c), s);
        Free(s);
        if(i != List_Nbr((yyvsp[-1].l)) - 1) strcat((yyval.c), "\n");
      }
      List_Delete((yyvsp[-1].l));
    }
#line 13678 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 577:
#line 6154 "Gmsh.y" /* yacc.c:1652  */
    {
      int i = 0;
      while ((yyvsp[-1].c)[i]) {
        (yyvsp[-1].c)[i] = toupper((yyvsp[-1].c)[i]);
        i++;
      }
      (yyval.c) = (yyvsp[-1].c);
    }
#line 13691 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 578:
#line 6163 "Gmsh.y" /* yacc.c:1652  */
    {
      int i = 0;
      while ((yyvsp[-1].c)[i]) {
        (yyvsp[-1].c)[i] = tolower((yyvsp[-1].c)[i]);
        i++;
      }
      (yyval.c) = (yyvsp[-1].c);
    }
#line 13704 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 579:
#line 6172 "Gmsh.y" /* yacc.c:1652  */
    {
      int i = 0;
      while ((yyvsp[-1].c)[i]) {
        if (i > 0 && (yyvsp[-1].c)[i-1] != '_')
          (yyvsp[-1].c)[i] = tolower((yyvsp[-1].c)[i]);
        i++;
      }
      (yyval.c) = (yyvsp[-1].c);
    }
#line 13718 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 580:
#line 6182 "Gmsh.y" /* yacc.c:1652  */
    {
      if((yyvsp[-5].d)){
        (yyval.c) = (yyvsp[-3].c);
        Free((yyvsp[-1].c));
      }
      else{
        (yyval.c) = (yyvsp[-1].c);
        Free((yyvsp[-3].c));
      }
    }
#line 13733 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 581:
#line 6193 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string in = (yyvsp[-5].c);
      std::string out = in.substr((int)(yyvsp[-3].d), (int)(yyvsp[-1].d));
      (yyval.c) = (char *)Malloc((out.size() + 1) * sizeof(char));
      strcpy((yyval.c), out.c_str());
      Free((yyvsp[-5].c));
    }
#line 13745 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 582:
#line 6201 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string in = (yyvsp[-3].c);
      std::string out = in.substr((int)(yyvsp[-1].d), std::string::npos);
      (yyval.c) = (char *)Malloc((out.size() + 1) * sizeof(char));
      strcpy((yyval.c), out.c_str());
      Free((yyvsp[-3].c));
    }
#line 13757 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 583:
#line 6209 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.c) = (yyvsp[-1].c);
    }
#line 13765 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 584:
#line 6213 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstring[5000];
      int i = printListOfDouble((yyvsp[-3].c), (yyvsp[-1].l), tmpstring);
      if(i < 0){
	yymsg(0, "Too few arguments in Sprintf");
	(yyval.c) = (yyvsp[-3].c);
      }
      else if(i > 0){
	yymsg(0, "%d extra argument%s in Sprintf", i, (i > 1) ? "s" : "");
	(yyval.c) = (yyvsp[-3].c);
      }
      else{
	(yyval.c) = (char*)Malloc((strlen(tmpstring) + 1) * sizeof(char));
	strcpy((yyval.c), tmpstring);
	Free((yyvsp[-3].c));
      }
      List_Delete((yyvsp[-1].l));
    }
#line 13788 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 585:
#line 6232 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp = FixRelativePath(gmsh_yyname, (yyvsp[-1].c));
      (yyval.c) = (char*)Malloc((tmp.size() + 1) * sizeof(char));
      strcpy((yyval.c), tmp.c_str());
      Free((yyvsp[-1].c));
    }
#line 13799 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 586:
#line 6239 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp = SplitFileName(GetAbsolutePath(gmsh_yyname))[0];
      (yyval.c) = (char*)Malloc((tmp.size() + 1) * sizeof(char));
      strcpy((yyval.c), tmp.c_str());
    }
#line 13809 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 587:
#line 6245 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp = SplitFileName((yyvsp[-1].c))[0];
      (yyval.c) = (char*)Malloc((tmp.size() + 1) * sizeof(char));
      strcpy((yyval.c), tmp.c_str());
      Free((yyvsp[-1].c));
    }
#line 13820 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 588:
#line 6252 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string tmp = GetAbsolutePath((yyvsp[-1].c));
      (yyval.c) = (char*)Malloc((tmp.size() + 1) * sizeof(char));
      strcpy((yyval.c), tmp.c_str());
      Free((yyvsp[-1].c));
    }
#line 13831 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 589:
#line 6259 "Gmsh.y" /* yacc.c:1652  */
    { init_options(); }
#line 13837 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 590:
#line 6261 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string val((yyvsp[-3].c));
      Msg::ExchangeOnelabParameter("", val, floatOptions, charOptions);
      (yyval.c) = (char*)Malloc((val.size() + 1) * sizeof(char));
      strcpy((yyval.c), val.c_str());
      Free((yyvsp[-3].c));
    }
#line 13849 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 591:
#line 6269 "Gmsh.y" /* yacc.c:1652  */
    {
      std::string out;
      const std::string * key_struct = nullptr;
      switch (gmsh_yynamespaces.get_key_struct_from_tag(struct_namespace,
                                                        (int)(yyvsp[-1].d), key_struct)) {
      case 0:
        out = *key_struct;
        break;
      case 1:
        yymsg(1, "Unknown NameSpace '%s' of Struct", struct_namespace.c_str());
        break;
      case 2:
        yymsg(1, "Unknown Struct of index %d", (int)(yyvsp[-1].d));
        break;
      default:
        break;
      }
      (yyval.c) = (char*)Malloc((out.size() + 1) * sizeof(char));
      strcpy((yyval.c), out.c_str());
    }
#line 13874 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 592:
#line 6293 "Gmsh.y" /* yacc.c:1652  */
    { struct_namespace = std::string(""); (yyval.d) = (yyvsp[0].d); }
#line 13880 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 593:
#line 6295 "Gmsh.y" /* yacc.c:1652  */
    { struct_namespace = (yyvsp[-3].c); Free((yyvsp[-3].c)); (yyval.d) = (yyvsp[0].d); }
#line 13886 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 594:
#line 6301 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[-1].l); }
#line 13892 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 595:
#line 6306 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[0].l); }
#line 13898 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 596:
#line 6308 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[0].l); }
#line 13904 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 597:
#line 6313 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[-1].l); }
#line 13910 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 598:
#line 6318 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20,20,sizeof(char*));
      List_Add((yyval.l), &((yyvsp[0].c)));
    }
#line 13919 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 599:
#line 6323 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.l) = (yyvsp[0].l); }
#line 13925 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 600:
#line 6325 "Gmsh.y" /* yacc.c:1652  */
    {
      List_Add((yyval.l), &((yyvsp[0].c)));
    }
#line 13933 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 601:
#line 6329 "Gmsh.y" /* yacc.c:1652  */
    {
      for(int i = 0; i < List_Nbr((yyvsp[0].l)); i++){
	char* c;
	List_Read((yyvsp[0].l), i, &c);
	List_Add((yyval.l), &c);
      }
      List_Delete((yyvsp[0].l));
    }
#line 13946 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 602:
#line 6341 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = List_Create(20, 20, sizeof(char *));
      if(!gmsh_yystringsymbols.count((yyvsp[-2].c)))
	yymsg(0, "Unknown string variable '%s'", (yyvsp[-2].c));
      else{
        std::vector<std::string> &s(gmsh_yystringsymbols[(yyvsp[-2].c)]);
	for(std::size_t i = 0; i < s.size(); i++) {
          char * val_ = strsave((char*)s.at(i).c_str());
	  List_Add((yyval.l), &val_);
        }
      }
      Free((yyvsp[-2].c));
    }
#line 13964 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 603:
#line 6355 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = treat_Struct_FullName_dot_tSTRING_ListOfString(nullptr, (yyvsp[-4].c), (yyvsp[-2].c));
    }
#line 13972 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 604:
#line 6359 "Gmsh.y" /* yacc.c:1652  */
    {
      (yyval.l) = treat_Struct_FullName_dot_tSTRING_ListOfString((yyvsp[-6].c), (yyvsp[-4].c), (yyvsp[-2].c));
    }
#line 13980 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 605:
#line 6366 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstr[256];
      sprintf(tmpstr, "_%d", (int)(yyvsp[-1].d));
      (yyval.c) = (char *)Malloc((strlen((yyvsp[-4].c))+strlen(tmpstr)+1)*sizeof(char));
      strcpy((yyval.c), (yyvsp[-4].c)); strcat((yyval.c), tmpstr);
      Free((yyvsp[-4].c));
    }
#line 13992 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 606:
#line 6374 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstr[256];
      sprintf(tmpstr, "_%d", (int)(yyvsp[-1].d));
      (yyval.c) = (char *)Malloc((strlen((yyvsp[-4].c))+strlen(tmpstr)+1)*sizeof(char)) ;
      strcpy((yyval.c), (yyvsp[-4].c)) ; strcat((yyval.c), tmpstr) ;
      Free((yyvsp[-4].c));
    }
#line 14004 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 607:
#line 6382 "Gmsh.y" /* yacc.c:1652  */
    {
      char tmpstr[256];
      sprintf(tmpstr, "_%d", (int)(yyvsp[-1].d));
      (yyval.c) = (char *)Malloc((strlen((yyvsp[-5].c))+strlen(tmpstr)+1)*sizeof(char));
      strcpy((yyval.c), (yyvsp[-5].c)); strcat((yyval.c), tmpstr);
      Free((yyvsp[-5].c));
    }
#line 14016 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 608:
#line 6393 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (yyvsp[0].c); }
#line 14022 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 609:
#line 6395 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (yyvsp[0].c); }
#line 14028 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;

  case 610:
#line 6398 "Gmsh.y" /* yacc.c:1652  */
    { (yyval.c) = (yyvsp[-1].c); }
#line 14034 "Gmsh.tab.cpp" /* yacc.c:1652  */
    break;


#line 14038 "Gmsh.tab.cpp" /* yacc.c:1652  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;


#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif


/*-----------------------------------------------------.
| yyreturn -- parsing is finished, return the result.  |
`-----------------------------------------------------*/
yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 6401 "Gmsh.y" /* yacc.c:1918  */


void assignVariable(const std::string &name, int index, int assignType,
                    double value)
{
  if(!gmsh_yysymbols.count(name)){
    if(!assignType){
      gmsh_yysymbol &s(gmsh_yysymbols[name]);
      s.list = true;
      s.value.resize(index + 1, 0.);
      s.value[index] = value;
    }
    else
      yymsg(0, "Unknown variable '%s'", name.c_str());
  }
  else{
    gmsh_yysymbol &s(gmsh_yysymbols[name]);
    if(s.list){
      if((int)s.value.size() < index + 1) s.value.resize(index + 1, 0.);
      switch(assignType){
      case 0 : s.value[index] = value; break;
      case 1 : s.value[index] += value; break;
      case 2 : s.value[index] -= value; break;
      case 3 : s.value[index] *= value; break;
      case 4 :
        if(value) s.value[index] /= value;
        else yymsg(0, "Division by zero in '%s[%d] /= %g'",
                   name.c_str(), index, value);
        break;
      }
    }
    else
      yymsg(0, "Variable '%s' is not a list", name.c_str());
  }
}

void assignVariables(const std::string &name, List_T *indices, int assignType,
                     List_T *values)
{
  if(List_Nbr(indices) != List_Nbr(values)){
    yymsg(0, "Incompatible array dimensions in affectation");
  }
  else{
    if(!gmsh_yysymbols.count(name)){
      if(!assignType){
        gmsh_yysymbol &s(gmsh_yysymbols[name]);
        s.list = true;
        for(int i = 0; i < List_Nbr(indices); i++){
          int index = (int)(*(double*)List_Pointer(indices, i));
          s.value.resize(index + 1, 0.);
          s.value[index] = *(double*)List_Pointer(values, i);
        }
      }
      else
        yymsg(0, "Unknown variable '%s'", name.c_str());
    }
    else{
      gmsh_yysymbol &s(gmsh_yysymbols[name]);
      if(s.list){
        for(int i = 0; i < List_Nbr(indices); i++){
          int index = (int)(*(double*)List_Pointer(indices, i));
          double d = *(double*)List_Pointer(values, i);
          if((int)s.value.size() < index + 1) s.value.resize(index + 1, 0.);
          switch(assignType){
          case 0 : s.value[index] = d; break;
          case 1 : s.value[index] += d; break;
          case 2 : s.value[index] -= d; break;
          case 3 : s.value[index] *= d; break;
          case 4 :
            if(d) s.value[index] /= d;
            else yymsg(0, "Division by zero in '%s[%d] /= %g'", name.c_str(), index, d);
            break;
          }
        }
      }
      else
        yymsg(0, "Variable '%s' is not a list", name.c_str());
    }
  }
}

void incrementVariable(const std::string &name, int index, double value)
{
  if(!gmsh_yysymbols.count(name))
    yymsg(0, "Unknown variable '%s'", name.c_str());
  else{
    gmsh_yysymbol &s(gmsh_yysymbols[name]);
    if(s.list){
      if((int)s.value.size() < index + 1) s.value.resize(index + 1, 0.);
      s.value[index] += value;
    }
    else
      yymsg(0, "Variable '%s' is not a list", name.c_str());
  }
}

int printListOfDouble(char *format, List_T *list, char *buffer)
{
  // if format does not contain formatting characters, dump the list (useful for
  // quick debugging of lists)
  int numFormats = 0;
  for(std::size_t i = 0; i < strlen(format); i++)
    if(format[i] == '%') numFormats++;
  if(!numFormats){
    strcpy(buffer, format);
    for(int i = 0; i < List_Nbr(list); i++){
      double d;
      List_Read(list, i, &d);
      char tmp[256];
      sprintf(tmp, " [%d]%g", i, d);
      strcat(buffer, tmp);
    }
    return 0;
  }

  char tmp1[256], tmp2[256];
  int j = 0, k = 0;
  buffer[j] = '\0';

  while(j < (int)strlen(format) && format[j] != '%') j++;
  strncpy(buffer, format, j);
  buffer[j] = '\0';
  for(int i = 0; i < List_Nbr(list); i++){
    k = j;
    j++;
    if(j < (int)strlen(format)){
      if(format[j] == '%'){
	strcat(buffer, "%");
	j++;
      }
      while(j < (int)strlen(format) && format[j] != '%') j++;
      if(k != j){
	strncpy(tmp1, &(format[k]), j-k);
	tmp1[j-k] = '\0';
	sprintf(tmp2, tmp1, *(double*)List_Pointer(list, i));
	strcat(buffer, tmp2);
      }
    }
    else
      return List_Nbr(list) - i;
  }
  if(j != (int)strlen(format))
    return -1;
  return 0;
}

void PrintParserSymbols(bool help, std::vector<std::string> &vec)
{
  if(help){
    vec.push_back("//");
    vec.push_back("// Numbers");
    vec.push_back("//");
  }
  for(std::map<std::string, gmsh_yysymbol>::iterator it = gmsh_yysymbols.begin();
      it != gmsh_yysymbols.end(); it++){
    gmsh_yysymbol s(it->second);
    std::ostringstream sstream;
    sstream.precision(12);
    sstream << it->first;
    if(s.list){
      sstream << "[] = {";
      for(std::size_t i = 0; i < s.value.size(); i++){
        if(i) sstream << ", ";
        sstream << s.value[i];
      }
      sstream << "}";
    }
    else
      sstream << " = " << s.value[0];
    sstream << ";";
    vec.push_back(sstream.str());
  }
  if(help){
    vec.push_back("//");
    vec.push_back("// Strings");
    vec.push_back("//");
  }
  for(std::map<std::string, std::vector<std::string> >::iterator it =
        gmsh_yystringsymbols.begin(); it != gmsh_yystringsymbols.end(); it++){
    if(it->second.size() == 1)
      vec.push_back(it->first + " = \"" + it->second[0] + "\";");
    else{
      std::string s = it->first + "[] = Str({";
      for(std::size_t i = 0; i < it->second.size(); i++){
        if(i) s += ", ";
        s += std::string("\"") + it->second[i] + "\"";
      }
      s += "});";
      vec.push_back(s);
    }
  }
  if (gmsh_yynamespaces.size()){
    if(help){
      vec.push_back("//");
      vec.push_back("// Structures");
      vec.push_back("//");
    }
    std::vector<std::string> strs;
    gmsh_yynamespaces.sprint(strs);
    vec.insert(vec.end(), strs.begin(), strs.end());
  }
}

fullMatrix<double> ListOfListOfDouble2Matrix(List_T *list)
{
  // Warning: this returns a fullMatrix copy, and deletes the input list
  int M = List_Nbr(list);
  int N = 0;
  for(int i = 0; i < M; i++){
    List_T *line = *(List_T**)List_Pointer_Fast(list, i);
    N = std::max(N, List_Nbr(line));
  }
  fullMatrix<double> mat(M, N);
  for(int i = 0; i < M; i++){
    List_T *line = *(List_T**)List_Pointer_Fast(list, i);
    for(int j = 0; j < List_Nbr(line); j++){
      double val;
      List_Read(line, j, &val);
      mat(i, j) = val;
    }
  }
  for(int i = 0; i < List_Nbr(list); i++)
    List_Delete(*(List_T**)List_Pointer(list, i));
  List_Delete(list);
  return mat;
}

void ListOfDouble2Vector(List_T *list, std::vector<int> &v)
{
  v.clear();
  if(!list) return;
  v.reserve(List_Nbr(list));
  for(int i = 0; i < List_Nbr(list); i++){
    double d;
    List_Read(list, i, &d);
    v.push_back((int)d);
  }
}

void ListOfDouble2Vector(List_T *list, std::vector<double> &v)
{
  v.clear();
  if(!list) return;
  v.reserve(List_Nbr(list));
  for(int i = 0; i < List_Nbr(list); i++){
    double d;
    List_Read(list, i, &d);
    v.push_back(d);
  }
}

void ListOfShapes2VectorOfPairs(List_T *list, std::vector<std::pair<int, int> > &v)
{
  for(int i = 0; i < List_Nbr(list); i++){
    Shape s;
    List_Read(list, i, &s);
    int dim = s.Type / 100 - 1;
    if(dim >= 0 && dim <= 3) v.push_back(std::pair<int, int>(dim, s.Num));
  }
}

void VectorOfPairs2ListOfShapes(const std::vector<std::pair<int, int> > &v, List_T *list)
{
  for(std::size_t i = 0; i < v.size(); i++){
    int dim = v[i].first;
    int tag = v[i].second;
    Shape s;
    s.Type = (dim == 3) ? MSH_VOLUME : (dim == 2) ? MSH_SURF_PLAN :
      (dim == 1) ? MSH_SEGM_LINE : MSH_POINT;
    s.Num = tag;
    List_Add(list, &s);
  }
}

void yyerror(const char *s)
{
  Msg::Error("'%s', line %d : %s (%s)", gmsh_yyname.c_str(), gmsh_yylineno - 1,
             s, gmsh_yytext);
  gmsh_yyerrorstate++;
}

void yymsg(int level, const char *fmt, ...)
{
  va_list args;
  char tmp[1024];

  va_start(args, fmt);
  vsprintf(tmp, fmt, args);
  va_end(args);

  if(level == 0){
    Msg::Error("'%s', line %d : %s", gmsh_yyname.c_str(), gmsh_yylineno - 1, tmp);
    gmsh_yyerrorstate++;
  }
  else if(level == 1){
    Msg::Warning("'%s', line %d : %s", gmsh_yyname.c_str(), gmsh_yylineno - 1, tmp);
  }
  else{
    Msg::Info("'%s', line %d : %s", gmsh_yyname.c_str(), gmsh_yylineno - 1, tmp);
  }
}

void addPeriodicFace(int iTarget, int iSource,
                     const std::vector<double>& affineTransform)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  GFace *target = GModel::current()->getFaceByTag(std::abs(iTarget));
  GFace *source = GModel::current()->getFaceByTag(std::abs(iSource));
  if (!target || !source) {
    Msg::Error("Could not find curve slave %d or master %d for periodic copy",
               iTarget, iSource);
  }
  else target->setMeshMaster(source, affineTransform);
}

void addPeriodicFace(int iTarget, int iSource,
                     const std::map<int,int>& edgeCounterparts)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  Msg::Info("Encoding periodic connection between %d and %d", iTarget, iSource);
  std::map<int,int>::const_iterator sIter = edgeCounterparts.begin();
  for (; sIter != edgeCounterparts.end(); ++sIter) {
    Msg::Info("%d - %d", sIter->first, sIter->second);
  }

  GFace *target = GModel::current()->getFaceByTag(std::abs(iTarget));
  GFace *source = GModel::current()->getFaceByTag(std::abs(iSource));
  if (!target || !source) {
    Msg::Error("Could not find surface slave %d or master %d for periodic copy",
               iTarget,iSource);
  }
  else target->setMeshMaster(source, edgeCounterparts);
}

void addPeriodicEdge(int iTarget,int iSource,
                     const std::vector<double>& affineTransform)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  GEdge *target = GModel::current()->getEdgeByTag(std::abs(iTarget));
  GEdge *source = GModel::current()->getEdgeByTag(std::abs(iSource));
  if (!target || !source)
    Msg::Error("Could not find surface %d or %d for periodic copy",
               iTarget,iSource);
  if (affineTransform.size() >= 12) {
    target->setMeshMaster(source, affineTransform);
  }
  else {
    target->setMeshMaster(source, iSource * iTarget < 0 ? -1 : 1);
  }
}

void computeAffineTransformation(SPoint3& origin, SPoint3& axis,
                                 double angle, SPoint3& translation,
                                 std::vector<double>& tfo)
{
  tfo.resize(16,0.0);

  double ca = cos(angle);
  double sa = sin(angle);

  double ux = axis.x();
  double uy = axis.y();
  double uz = axis.z();

  tfo.resize(16);

  tfo[0*4+0] = ca + ux*ux*(1.-ca);
  tfo[0*4+1] = ux*uy*(1.-ca) - uz * sa;
  tfo[0*4+2] = ux*uz*(1.-ca) + uy * sa;

  tfo[1*4+0] = ux*uy*(1.-ca) + uz * sa;
  tfo[1*4+1] = ca + uy*uy*(1.-ca);
  tfo[1*4+2] = uy*uz*(1.-ca) - ux * sa;

  tfo[2*4+0] = ux*uz*(1.-ca) - uy * sa;
  tfo[2*4+1] = uy*uz*(1.-ca) + ux * sa;
  tfo[2*4+2] = ca + uz*uz*(1.-ca);

  int idx = 0;
  for (size_t i = 0; i < 3; i++,idx++) {
    int tIdx = i*4+3;
    tfo[tIdx] = origin[i] + translation[i];
    for (int j = 0; j < 3; j++,idx++) tfo[tIdx] -= tfo[idx] * origin[j];
  }

  for (int i = 0; i < 4; i++) tfo[12+i] = 0;
  tfo[15] = 1;
}

void addEmbedded(int dim, std::vector<int> tags, int dim2, int tag2)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  if(dim2 == 2){
    GFace *gf = GModel::current()->getFaceByTag(tag2);
    if(!gf){
      yymsg(0, "Unknown model surface with tag %d", tag2);
      return;
    }
    for(std::size_t i = 0; i < tags.size(); i++){
      if(dim == 0){
        GVertex *gv = GModel::current()->getVertexByTag(tags[i]);
        if(gv)
          gf->addEmbeddedVertex(gv);
        else
          yymsg(0, "Unknown model point %d", tags[i]);
      }
      else if(dim == 1){
        GEdge *ge = GModel::current()->getEdgeByTag(tags[i]);
        if(ge)
          gf->addEmbeddedEdge(ge);
        else
          yymsg(0, "Unknown model curve %d", tags[i]);
      }
    }
  }
  else if(dim2 == 3){
    GRegion *gr = GModel::current()->getRegionByTag(tag2);
    if(!gr){
      yymsg(0, "Unknown model volume with tag %d", tag2);
      return;
    }
    for(std::size_t i = 0; i < tags.size(); i++){
      if(dim == 0){
        GVertex *gv = GModel::current()->getVertexByTag(tags[i]);
        if(gv)
          gr->addEmbeddedVertex(gv);
        else
          yymsg(0, "Unknown model point with tag %d", tags[i]);
      }
      else if(dim == 1){
        GEdge *ge = GModel::current()->getEdgeByTag(tags[i]);
        if(ge)
          gr->addEmbeddedEdge(ge);
        else
          yymsg(0, "Unknown model curve with tag %d", tags[i]);
      }
      else if(dim == 2){
        GFace *gf = GModel::current()->getFaceByTag(tags[i]);
        if(gf)
          gr->addEmbeddedFace(gf);
        else
          yymsg(0, "Unknown model surface with tag %d", tags[i]);
      }
    }
  }
}

void getAllElementaryTags(int dim, List_T *out)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  std::vector<GEntity*> entities;
  GModel::current()->getEntities(entities, dim);
  for(std::size_t i = 0; i < entities.size(); i++){
    double tag = entities[i]->tag();
    List_Add(out, &tag);
  }
}

void getAllPhysicalTags(int dim, List_T *out)
{
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  std::map<int, std::vector<GEntity*> > groups;
  GModel::current()->getPhysicalGroups(dim, groups);
  for(std::map<int, std::vector<GEntity*> >::iterator it = groups.begin();
      it != groups.end(); it++){
    double d = it->first;
    List_Add(out, &d);
  }
}

void getElementaryTagsForPhysicalGroups(int dim, List_T *in, List_T *out)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  std::map<int, std::vector<GEntity*> > groups;
  GModel::current()->getPhysicalGroups(dim, groups);
  for(int i = 0; i < List_Nbr(in); i++){
    double num;
    List_Read(in, i, &num);
    std::map<int, std::vector<GEntity*> >::iterator it = groups.find(num);
    if(it != groups.end()){
      for(unsigned j = 0; j < it->second.size(); j++){
        double d = it->second[j]->tag();
        List_Add(out, &d);
      }
    }
  }
}

void getElementaryTagsInBoundingBox(int dim, double x1, double y1, double z1,
                                    double x2, double y2, double z2, List_T *out)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  SBoundingBox3d box(x1, y1, z1, x2, y2, z2);
  std::vector<GEntity*> entities;
  GModel::current()->getEntitiesInBox(entities, box, dim);
  for(std::size_t i = 0; i < entities.size(); i++){
    double d = entities[i]->tag();
    List_Add(out, &d);
  }
}

void getParentTags(int dim, List_T *in, List_T *out)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  for(int i = 0; i < List_Nbr(in); i++){
    double num;
    List_Read(in, i, &num);
    GEntity *ge = GModel::current()->getEntityByTag(dim, (int)num);
    if(ge){
      GEntity *parent = ge->getParentEntity();
      if(parent){
        double tag = parent->tag();
        List_Add(out, &tag);
      }
    }
  }
}

void getBoundingBox(int dim, List_T *in, List_T *out)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  SBoundingBox3d box;
  for(int i = 0; i < List_Nbr(in); i++){
    double num;
    List_Read(in, i, &num);
    GEntity *ge = GModel::current()->getEntityByTag(dim, (int)num);
    if(ge) box += ge->bounds();
  }
  if(!box.empty()){
    double b[6] = {box.min().x(), box.min().y(), box.min().z(),
                   box.max().x(), box.max().y(), box.max().z()};
    for(int i = 0; i < 6; i++)
      List_Add(out, &b[i]);
  }
}

void setVisibility(int dim, int visible, bool recursive)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  std::vector<GEntity*> entities;
  GModel::current()->getEntities(entities, dim);
  for(std::size_t i = 0; i < entities.size(); i++){
    entities[i]->setVisibility(visible);
  }
}

void setVisibility(const std::vector<std::pair<int, int> > &dimTags,
                   int visible, bool recursive)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  for(std::size_t i = 0; i < dimTags.size(); i++){
    GEntity *ge = GModel::current()->getEntityByTag
      (dimTags[i].first, std::abs(dimTags[i].second));
    if(ge) ge->setVisibility(visible, recursive);
  }
}

void setColor(const std::vector<std::pair<int, int> > &dimTags,
              unsigned int val, bool recursive)
{
  if(GModel::current()->getOCCInternals() &&
     GModel::current()->getOCCInternals()->getChanged())
    GModel::current()->getOCCInternals()->synchronize(GModel::current());
  if(GModel::current()->getGEOInternals()->getChanged())
    GModel::current()->getGEOInternals()->synchronize(GModel::current());

  for(std::size_t i = 0; i < dimTags.size(); i++){
    GEntity *ge = GModel::current()->getEntityByTag
      (dimTags[i].first, std::abs(dimTags[i].second));
    if(ge) ge->setColor(val, recursive);
  }
}

double treat_Struct_FullName_Float
(char* c1, char* c2, int type_var, int index, double val_default, int type_treat)
{
  double out;
  if(!c1 && gmsh_yysymbols.count(c2)){
    if (type_treat == 1) out = 1.; // Exists (type_treat == 1)
    else { // Get (0) or GetForced (2)
      if (type_var == 1) {
        gmsh_yysymbol &s(gmsh_yysymbols[c2]);
        if(s.value.empty()){
          out = val_default;
          if (type_treat == 0) yymsg(0, "Uninitialized variable '%s'", c2);
        }
        else
          out = s.value[0];
      }
      else if (type_var == 2) {
        gmsh_yysymbol &s(gmsh_yysymbols[c2]);
        if(index < 0 || (int)s.value.size() < index + 1){
          out = val_default;
          if (type_treat == 0) yymsg(0, "Uninitialized variable '%s[%d]'", c2, index);
        }
        else{
          out = s.value[index];
        }
      }
      else {
        out = val_default;
      }
    }
  }
  else if(!c1 && type_treat == 1 && gmsh_yystringsymbols.count(c2)) {
    out = 1.;
  }
  else{
    if (type_var == 1) {
      std::string struct_namespace(c1? c1 : std::string("")), struct_name(c2);
      if(gmsh_yynamespaces.getTag(struct_namespace, struct_name, out)) {
        out = val_default;
        if (type_treat == 0) yymsg(0, "Unknown variable '%s'", struct_name.c_str());
      }
    }
    else {
      out = val_default;
      if (type_treat == 0) yymsg(0, "Unknown variable '%s(.)'", c2);
    }
  }
  Free(c1); Free(c2);
  return out;
}

double treat_Struct_FullName_dot_tSTRING_Float
(char* c1, char* c2, char* c3, int index, double val_default, int type_treat)
{
  double out;
  std::string struct_namespace(c1? c1 : std::string("")), struct_name(c2);
  std::string key_member(c3);
  switch (gmsh_yynamespaces.getMember
          (struct_namespace, struct_name, key_member, out, index)) {
  case 0:
    if (type_treat == 1) out = 1.; // Exists (type_treat == 1)
    break;
  case 1:
    if (!NumberOption(GMSH_GET, c2, 0, c3, out, type_treat==0))
      out = val_default;
    break;
  case 2:
    if (type_treat != 0) {
      const std::string * out_dummy = nullptr;
      out = (gmsh_yynamespaces.getMember
             (struct_namespace, struct_name, key_member, out_dummy))?
        val_default : 1.;
    }
    else {
      out = val_default;
      if (type_treat == 0)
        yymsg(0, "Unknown member '%s' of Struct %s", c3, struct_name.c_str());
    }
    break;
  case 3:
    out = val_default;
    if (type_treat == 0)
      yymsg(0, "Index %d out of range", index);
    break;
  }
  Free(c1); Free(c2);
  if (flag_tSTRING_alloc) Free(c3);
  return out;
}

List_T * treat_Struct_FullName_dot_tSTRING_ListOfFloat
(char* c1, char* c2, char* c3)
{
  List_T * out, * val_default = nullptr;
  const std::vector<double> * out_vector; double val_;
  std::string struct_namespace(c1? c1 : std::string("")), struct_name(c2);
  std::string key_member(c3);
  switch (gmsh_yynamespaces.getMember_Vector
          (struct_namespace, struct_name, key_member, out_vector)) {
  case 0:
    out = List_Create(out_vector->size(), 1, sizeof(double));
    for(std::size_t i = 0; i < out_vector->size(); i++) {
      val_ = out_vector->at(i);
      List_Add(out, &val_);
    }
    break;
  case 1:
    yymsg(0, "Unknown Struct: %s", struct_name.c_str());
    out = val_default;
    break;
  case 2:
    out = val_default;
    yymsg(0, "Unknown member '%s' of Struct %s", c3, struct_name.c_str());
    break;
  }
  Free(c1); Free(c2);
  if (flag_tSTRING_alloc) Free(c3);
  return out;
}

int treat_Struct_FullName_dot_tSTRING_Float_getDim
(char* c1, char* c2, char* c3)
{
  int out;
  std::string struct_namespace(c1? c1 : std::string("")), struct_name(c2);
  std::string key_member(c3);
  switch (gmsh_yynamespaces.getMember_Dim
          (struct_namespace, struct_name, key_member, out)) {
  case 0:
    break;
  case 1:
    out = 0;
    break;
  case 2:
    out = 0;
    yymsg(0, "Unknown member '%s' of Struct %s", c3, struct_name.c_str());
    break;
  }
  Free(c1); Free(c2);
  if (flag_tSTRING_alloc) Free(c3);
  return out;
}

char * treat_Struct_FullName_String
(char* c1, char* c2, int type_var, int index, char * val_default, int type_treat)
{
  std::string string_default(val_default? val_default : std::string(""));
  const std::string * out = nullptr;
  std::string out_tmp;
  if(!c1 && gmsh_yystringsymbols.count(c2)){
    // Get (0) or GetForced (2)
    if(gmsh_yystringsymbols[c2].size() != 1){
      out = &string_default;
      if (type_treat == 0)
        yymsg(0, "Expected single valued string variable '%s'", c2);
    }
    else {
      out_tmp = gmsh_yystringsymbols[c2][0];
      out = &out_tmp;
    }
  }
  else{
    out = &string_default;
    if (type_treat == 0) yymsg(0, "Unknown string variable '%s'", c2);
  }
  char* out_c = (char*)Malloc((out->size() + 1) * sizeof(char));
  strcpy(out_c, out->c_str());
  Free(c1); Free(c2);
  return out_c;
}

char* treat_Struct_FullName_dot_tSTRING_String
(char* c1, char* c2, char* c3, int index, char * val_default, int type_treat)
{
  std::string string_default(val_default? val_default : std::string(""));
  const std::string * out = nullptr;
  std::string out_tmp; // PD: we should avoid that -> StringOption() to be changed
  std::string struct_namespace(c1? c1 : std::string("")), struct_name(c2);
  std::string key_member(c3);
  switch (gmsh_yynamespaces.getMember
          (struct_namespace, struct_name, key_member, out, index)) {
  case 0:
    break;
  case 1:
    if (StringOption(GMSH_GET, c2, 0, c3, out_tmp, type_treat==0))
      out = &out_tmp;
    else
      out = &string_default;
    break;
  case 2:
    out = &string_default;
    if (type_treat == 0)
      yymsg(0, "Unknown member '%s' of Struct %s", c3, struct_name.c_str());
    break;
  case 3:
    out = &string_default;
    if (type_treat == 0)
      yymsg(0, "Index %d out of range", index);
    break;
  }
  char* out_c = (char*)Malloc((out->size() + 1) * sizeof(char));
  strcpy(out_c, out->c_str());
  Free(c1); Free(c2);
  if (flag_tSTRING_alloc) Free(c3);
  return out_c;
}

List_T * treat_Struct_FullName_dot_tSTRING_ListOfString
(char* c1, char* c2, char* c3)
{
  List_T * out, * val_default = nullptr;
  const std::vector<std::string> * out_vector; char * val_;
  std::string struct_namespace(c1? c1 : std::string("")), struct_name(c2);
  std::string key_member(c3);
  switch (gmsh_yynamespaces.getMember_Vector
          (struct_namespace, struct_name, key_member, out_vector)) {
  case 0:
    out = List_Create(out_vector->size(), 1, sizeof(char *));
    for(std::size_t i = 0; i < out_vector->size(); i++) {
      val_ = strsave((char*)out_vector->at(i).c_str());
      List_Add(out, &val_);
    }
    break;
  case 1:
    yymsg(0, "Unknown Struct: %s", struct_name.c_str());
    out = val_default;
    break;
  case 2:
    out = val_default;
    yymsg(0, "Unknown member '%s' of Struct %s", c3, struct_name.c_str());
    break;
  }
  Free(c1); Free(c2);
  if (flag_tSTRING_alloc) Free(c3);
  return out;
}
