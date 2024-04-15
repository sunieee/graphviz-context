/* A Bison parser, made by GNU Bison 3.5.1.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2020 Free Software Foundation,
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
#define YYBISON_VERSION "3.5.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Substitute the type names.  */
#define YYSTYPE         EX_STYPE
/* Substitute the variable and function names.  */
#define yyparse         ex_parse
#define yylex           ex_lex
#define yyerror         ex_error
#define yydebug         ex_debug
#define yynerrs         ex_nerrs
#define yylval          ex_lval
#define yychar          ex_char

/* First part of user prologue.  */
#line 19 "../../lib/expr/exparse.y"


/*
 * Glenn Fowler
 * AT&T Research
 *
 * expression library grammar and compiler
 */

#include <assert.h>
#include <ctype.h>
#include <expr/exop.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ast/ast.h>


#line 97 "exparse.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
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

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
#ifndef YY_EX_EXPARSE_H_INCLUDED
# define YY_EX_EXPARSE_H_INCLUDED
/* Debug traces.  */
#ifndef EX_DEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define EX_DEBUG 1
#  else
#   define EX_DEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define EX_DEBUG 1
# endif /* ! defined YYDEBUG */
#endif  /* ! defined EX_DEBUG */
#if EX_DEBUG
extern int ex_debug;
#endif

/* Token type.  */
#ifndef EX_TOKENTYPE
# define EX_TOKENTYPE
  enum ex_tokentype
  {
    MINTOKEN = 258,
    INTEGER = 259,
    UNSIGNED = 260,
    CHARACTER = 261,
    FLOATING = 262,
    STRING = 263,
    VOIDTYPE = 264,
    STATIC = 265,
    ADDRESS = 266,
    ARRAY = 267,
    BREAK = 268,
    CALL = 269,
    CASE = 270,
    CONSTANT = 271,
    CONTINUE = 272,
    DECLARE = 273,
    DEFAULT = 274,
    DYNAMIC = 275,
    ELSE = 276,
    EXIT = 277,
    FOR = 278,
    FUNCTION = 279,
    GSUB = 280,
    ITERATE = 281,
    ITERATER = 282,
    ID = 283,
    IF = 284,
    LABEL = 285,
    MEMBER = 286,
    NAME = 287,
    POS = 288,
    PRAGMA = 289,
    PRE = 290,
    PRINT = 291,
    PRINTF = 292,
    PROCEDURE = 293,
    QUERY = 294,
    RAND = 295,
    RETURN = 296,
    SCANF = 297,
    SPLIT = 298,
    SPRINTF = 299,
    SRAND = 300,
    SSCANF = 301,
    SUB = 302,
    SUBSTR = 303,
    SWITCH = 304,
    TOKENS = 305,
    UNSET = 306,
    WHILE = 307,
    F2I = 308,
    F2S = 309,
    I2F = 310,
    I2S = 311,
    S2B = 312,
    S2F = 313,
    S2I = 314,
    F2X = 315,
    I2X = 316,
    S2X = 317,
    X2F = 318,
    X2I = 319,
    X2S = 320,
    X2X = 321,
    XPRINT = 322,
    OR = 323,
    AND = 324,
    EQ = 325,
    NE = 326,
    LE = 327,
    GE = 328,
    LSH = 329,
    RSH = 330,
    IN_OP = 331,
    UNARY = 332,
    INC = 333,
    DEC = 334,
    CAST = 335,
    MAXTOKEN = 336
  };
#endif
/* Tokens.  */
#define MINTOKEN 258
#define INTEGER 259
#define UNSIGNED 260
#define CHARACTER 261
#define FLOATING 262
#define STRING 263
#define VOIDTYPE 264
#define STATIC 265
#define ADDRESS 266
#define ARRAY 267
#define BREAK 268
#define CALL 269
#define CASE 270
#define CONSTANT 271
#define CONTINUE 272
#define DECLARE 273
#define DEFAULT 274
#define DYNAMIC 275
#define ELSE 276
#define EXIT 277
#define FOR 278
#define FUNCTION 279
#define GSUB 280
#define ITERATE 281
#define ITERATER 282
#define ID 283
#define IF 284
#define LABEL 285
#define MEMBER 286
#define NAME 287
#define POS 288
#define PRAGMA 289
#define PRE 290
#define PRINT 291
#define PRINTF 292
#define PROCEDURE 293
#define QUERY 294
#define RAND 295
#define RETURN 296
#define SCANF 297
#define SPLIT 298
#define SPRINTF 299
#define SRAND 300
#define SSCANF 301
#define SUB 302
#define SUBSTR 303
#define SWITCH 304
#define TOKENS 305
#define UNSET 306
#define WHILE 307
#define F2I 308
#define F2S 309
#define I2F 310
#define I2S 311
#define S2B 312
#define S2F 313
#define S2I 314
#define F2X 315
#define I2X 316
#define S2X 317
#define X2F 318
#define X2I 319
#define X2S 320
#define X2X 321
#define XPRINT 322
#define OR 323
#define AND 324
#define EQ 325
#define NE 326
#define LE 327
#define GE 328
#define LSH 329
#define RSH 330
#define IN_OP 331
#define UNARY 332
#define INC 333
#define DEC 334
#define CAST 335
#define MAXTOKEN 336

/* Value type.  */
#if ! defined EX_STYPE && ! defined EX_STYPE_IS_DECLARED
union EX_STYPE
{
#line 39 "../../lib/expr/exparse.y"

	struct Exnode_s*expr;
	double		floating;
	struct Exref_s*	reference;
	struct Exid_s*	id;
	Sflong_t	integer;
	int		op;
	char*		string;

#line 329 "exparse.c"

};
typedef union EX_STYPE EX_STYPE;
# define EX_STYPE_IS_TRIVIAL 1
# define EX_STYPE_IS_DECLARED 1
#endif


extern EX_STYPE ex_lval;

int ex_parse (void);

#endif /* !YY_EX_EXPARSE_H_INCLUDED  */

/* Second part of user prologue.  */
#line 169 "../../lib/expr/exparse.y"


#include <expr/exgram.h>

void ex_error(const char *message);


#line 353 "exparse.c"


#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))

/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

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

#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                            \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
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

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

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
         || (defined EX_STYPE_IS_TRIVIAL && EX_STYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
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
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1112

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  107
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  44
/* YYNRULES -- Number of rules.  */
#define YYNRULES  142
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  286

#define YYUNDEFTOK  2
#define YYMAXUTOK   336


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    91,     2,    93,     2,    90,    76,     2,
      98,   103,    88,    85,    68,    86,   106,    89,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    71,   102,
      79,    69,    80,    70,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   104,     2,   105,    75,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   100,    74,   101,    92,     2,     2,     2,
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
      65,    66,    67,    72,    73,    77,    78,    81,    82,    83,
      84,    87,    94,    95,    96,    97,    99
};

#if EX_DEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   179,   179,   200,   201,   204,   204,   244,   247,   262,
     266,   270,   270,   270,   275,   285,   298,   313,   326,   334,
     345,   355,   355,   366,   378,   382,   395,   424,   427,   458,
     459,   462,   481,   488,   491,   497,   498,   505,   505,   554,
     555,   556,   557,   560,   561,   565,   568,   575,   578,   581,
     585,   589,   642,   646,   650,   654,   658,   662,   666,   670,
     674,   678,   682,   686,   690,   694,   698,   702,   715,   719,
     729,   729,   729,   770,   790,   797,   801,   805,   809,   813,
     817,   821,   825,   829,   833,   837,   841,   847,   851,   855,
     861,   866,   870,   895,   931,   951,   959,   967,   978,   982,
     986,   989,   990,   992,  1000,  1005,  1010,  1015,  1022,  1023,
    1024,  1027,  1028,  1031,  1035,  1055,  1068,  1071,  1075,  1089,
    1092,  1099,  1102,  1110,  1115,  1122,  1125,  1131,  1134,  1138,
    1149,  1149,  1162,  1165,  1178,  1199,  1203,  1209,  1212,  1219,
    1220,  1237,  1220
};
#endif

#if EX_DEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "MINTOKEN", "INTEGER", "UNSIGNED",
  "CHARACTER", "FLOATING", "STRING", "VOIDTYPE", "STATIC", "ADDRESS",
  "ARRAY", "BREAK", "CALL", "CASE", "CONSTANT", "CONTINUE", "DECLARE",
  "DEFAULT", "DYNAMIC", "ELSE", "EXIT", "FOR", "FUNCTION", "GSUB",
  "ITERATE", "ITERATER", "ID", "IF", "LABEL", "MEMBER", "NAME", "POS",
  "PRAGMA", "PRE", "PRINT", "PRINTF", "PROCEDURE", "QUERY", "RAND",
  "RETURN", "SCANF", "SPLIT", "SPRINTF", "SRAND", "SSCANF", "SUB",
  "SUBSTR", "SWITCH", "TOKENS", "UNSET", "WHILE", "F2I", "F2S", "I2F",
  "I2S", "S2B", "S2F", "S2I", "F2X", "I2X", "S2X", "X2F", "X2I", "X2S",
  "X2X", "XPRINT", "','", "'='", "'?'", "':'", "OR", "AND", "'|'", "'^'",
  "'&'", "EQ", "NE", "'<'", "'>'", "LE", "GE", "LSH", "RSH", "'+'", "'-'",
  "IN_OP", "'*'", "'/'", "'%'", "'!'", "'~'", "'#'", "UNARY", "INC", "DEC",
  "CAST", "'('", "MAXTOKEN", "'{'", "'}'", "';'", "')'", "'['", "']'",
  "'.'", "$accept", "program", "action_list", "action", "$@1",
  "statement_list", "statement", "$@2", "$@3", "$@4", "switch_list",
  "switch_item", "case_list", "case_item", "static", "dcl_list",
  "dcl_item", "$@5", "dcl_name", "name", "else_opt", "expr_opt", "expr",
  "$@6", "$@7", "splitop", "constant", "print", "scan", "variable",
  "array", "index", "args", "arg_list", "formals", "formal_list",
  "formal_item", "$@8", "members", "member", "assign", "initialize", "$@9",
  "$@10", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,    44,    61,
      63,    58,   323,   324,   124,    94,    38,   325,   326,    60,
      62,   327,   328,   329,   330,    43,    45,   331,    42,    47,
      37,    33,   126,    35,   332,   333,   334,   335,    40,   336,
     123,   125,    59,    41,    91,    93,    46
};
# endif

#define YYPACT_NINF (-144)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-127)

#define yytable_value_is_error(Yyn) \
  ((Yyn) == YYTABLE_NINF)

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    -144,     9,   200,  -144,  -144,  -144,  -144,  -144,  -144,   -89,
     691,  -144,   691,   -80,   -71,   -64,   -63,   -43,   -35,   -27,
     -11,  -144,    11,  -144,    16,  -144,    20,   691,  -144,  -144,
    -144,    23,  -144,    31,    39,    41,  -144,    56,    57,     1,
     691,   691,   691,   691,    79,     1,     1,   596,  -144,    92,
    -144,  -144,    49,   872,    58,  -144,    60,    61,   -37,   691,
      65,    66,   691,   -27,   691,   691,   691,   691,     1,   -12,
    -144,  -144,   691,   691,   691,    59,    68,    88,   691,   691,
     691,   140,   691,  -144,  -144,  -144,  -144,  -144,  -144,  -144,
    -144,    72,   284,   299,   100,  -144,   158,  -144,   691,  -144,
     691,   691,   691,   691,   691,   691,   691,   691,   691,   691,
     691,   691,   691,   691,   691,   143,   691,   691,   691,   691,
     691,   691,   691,  -144,  -144,  -144,   916,    77,   109,  -144,
    -144,   185,  -144,   383,    83,   -56,    84,    85,    95,    93,
    -144,   482,    97,    98,  -144,  -144,  -144,   577,   106,   108,
     872,   -51,   672,   691,  -144,  -144,  -144,  -144,   916,   691,
     934,   951,   967,   982,   996,  1010,  1010,  1022,  1022,  1022,
    1022,   107,   107,    53,    53,  -144,  -144,  -144,  -144,   895,
     111,   112,   916,  -144,   691,  -144,  -144,   691,   497,  -144,
    -144,   497,    29,  -144,   497,  -144,  -144,  -144,  -144,  -144,
     116,   691,  -144,   497,  -144,  -144,    87,   849,   786,  -144,
    -144,   916,   124,  -144,  -144,  -144,   168,    90,   767,  -144,
     200,  -144,  -144,  -144,  -144,   -49,  -144,  -144,  -144,   -54,
     691,   497,  -144,  -144,  -144,    87,  -144,   126,   691,   691,
    -144,   128,  -144,    -7,  -144,   -16,   -44,   916,   815,   497,
     145,   162,  -144,  -144,    86,  -144,   129,  -144,  -144,  -144,
    -144,  -144,  -144,   183,  -144,   200,  -144,  -144,   238,  -144,
     174,  -144,   210,  -144,   -10,   176,   262,  -144,  -144,  -144,
     181,  -144,  -144,  -144,   398,  -144
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       7,     0,     3,     1,   105,   107,   104,   106,    34,     0,
      47,   103,    47,   119,     0,     0,     0,     0,     0,   132,
       0,   115,     0,   108,     0,   109,     0,    47,   111,   101,
     110,     0,   112,     0,     0,     0,   102,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     7,     2,
       8,    11,     0,    48,     0,   100,     0,     0,   137,   121,
       0,     0,     0,   132,     0,    47,   121,   121,     0,     0,
     113,   133,     0,   121,   121,     0,     0,     0,   121,   121,
       0,     0,     0,    78,    77,    76,    73,    75,    74,    95,
      98,     0,     0,    33,     0,     4,     0,    10,     0,    70,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     121,   121,     0,    96,    99,    94,   123,     0,   122,    23,
      24,     0,   114,     0,     0,   137,     0,     0,     0,   135,
     136,     0,     0,     0,    87,    25,    88,     0,     0,     0,
      21,     0,     0,     0,    49,     9,     5,    12,    69,     0,
      68,    67,    64,    65,    63,    61,    62,    51,    58,    59,
      60,    56,    57,    66,    52,    97,    53,    54,    55,     0,
       0,     0,   138,    79,     0,   120,    86,    47,    33,    80,
      81,    33,     0,   134,    33,    91,    90,    89,    82,    83,
       0,     0,    18,    33,    50,     7,     0,     0,     0,    92,
      93,   124,     0,    15,    17,   135,    45,     0,     0,    20,
       6,    40,    42,    41,    39,     0,    35,    37,    71,   119,
      47,    33,    14,    26,    19,     0,    13,   116,     0,     0,
      84,     0,    46,     0,    36,     0,   137,    72,     0,    33,
       0,     0,    22,    27,     7,    29,     0,   117,   140,   139,
      38,    85,    16,     0,    32,    28,    30,   118,   125,    31,
     130,   141,   127,   128,     0,     0,     0,    44,    43,   131,
       0,   130,   129,     7,    33,   142
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -144,  -144,  -144,  -144,  -144,   -48,  -143,  -144,  -144,  -144,
    -144,  -144,  -144,    28,  -144,  -144,    48,  -144,  -144,  -144,
    -144,    -9,   -36,  -144,  -144,  -144,    34,  -144,  -144,   101,
    -144,  -144,    24,  -144,  -144,  -144,    12,  -144,   224,   150,
      51,  -144,  -144,  -144
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,    49,    95,   205,     2,    50,    96,   206,   200,
     243,   253,   254,   255,    51,   225,   226,   237,   227,   279,
     232,    52,    53,   159,   238,    54,    55,    56,    57,    58,
     246,    63,   127,   128,   271,   272,   273,   274,    70,    71,
     125,   260,   268,   275
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      93,    60,   256,    61,    84,    85,    86,    87,   250,     3,
     277,    92,   251,   122,   239,    59,   139,   201,    76,   235,
     140,    13,   278,   126,    62,   122,   131,    64,   133,    19,
     126,   126,   122,    21,    65,    66,   141,   126,   126,   123,
     124,   147,   126,   126,   150,   213,   152,   188,   214,   240,
      62,   216,   202,   236,   258,    67,   134,   215,   123,   124,
     219,   140,   158,    68,   160,   161,   162,   163,   164,   165,
     166,   167,   168,   169,   170,   171,   172,   173,   174,    69,
     176,   177,   178,   179,   126,   126,   182,    72,   242,   257,
     136,   137,     4,     5,   252,     6,     7,   142,   143,    88,
       9,   250,   148,   149,    11,   251,   262,   221,    13,    73,
      14,   222,    16,    17,    74,   223,    19,   204,    75,   224,
      21,    77,    94,   207,    22,    23,    24,    25,    26,    78,
      28,    29,    30,    31,    32,    33,    34,    79,    36,    80,
      83,   116,   117,   118,   180,   181,    89,    90,   211,     4,
       5,    97,     6,     7,    81,    82,   119,   220,   120,   121,
     151,    11,   144,   175,    39,   218,   135,   129,   130,   138,
     145,   156,   158,    40,    41,   153,   157,   184,   212,    42,
      43,    44,   183,    45,    46,   187,    47,   189,   190,   231,
     233,   146,   113,   114,   115,   116,   117,   118,   191,   192,
     195,   196,   247,   248,     4,     5,   265,     6,     7,   198,
       8,   199,     9,    10,   209,   210,    11,    12,   -33,   217,
      13,   241,    14,    15,    16,    17,   230,    18,    19,    20,
     245,   249,    21,   264,   267,   284,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    98,   269,    99,   270,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   109,   110,   111,   112,
     113,   114,   115,   116,   117,   118,    39,  -126,   276,   280,
     281,   283,   266,   244,   263,    40,    41,   132,   282,   193,
     185,    42,    43,    44,     0,    45,    46,   259,    47,     0,
      48,     0,   -47,     4,     5,     0,     6,     7,     0,     8,
       0,     9,    10,     0,     0,    11,    12,     0,     0,    13,
       0,    14,    15,    16,    17,     0,    18,    19,    20,     0,
       0,    21,     0,     0,     0,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    98,     0,    99,     0,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,   112,   113,
     114,   115,   116,   117,   118,    39,     0,     0,     0,     0,
       0,     0,     0,     0,    40,    41,     0,   154,     0,     0,
      42,    43,    44,     0,    45,    46,     0,    47,     0,    48,
     155,   -47,     4,     5,     0,     6,     7,     0,     8,     0,
       9,    10,     0,     0,    11,    12,     0,     0,    13,     0,
      14,    15,    16,    17,     0,    18,    19,    20,     0,     0,
      21,     0,     0,     0,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    98,     0,    99,     0,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,    39,     0,     0,     0,     0,     0,
       0,     0,     0,    40,    41,     0,   186,     0,     0,    42,
      43,    44,     0,    45,    46,     0,    47,     0,    48,   285,
     -47,     4,     5,     0,     6,     7,     0,     8,     0,     9,
      10,     0,     0,    11,    12,     0,     0,    13,     0,    14,
      15,    16,    17,     0,    18,    19,    20,     0,     0,    21,
       0,     0,     0,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      98,     0,    99,     0,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,    39,     0,     0,     0,     0,     0,     0,
       0,     0,    40,    41,     0,   194,     0,     0,    42,    43,
      44,     0,    45,    46,     0,    47,     0,    48,     0,   -47,
       4,     5,     0,     6,     7,     0,     0,     0,     9,     0,
       0,     0,    11,     0,    91,     0,    13,     0,    14,     0,
      16,    17,     0,     0,    19,     0,     0,     0,    21,     0,
       0,     0,    22,    23,    24,    25,    26,     0,    28,    29,
      30,    31,    32,    33,    34,    98,    36,    99,     0,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   109,   110,
     111,   112,   113,   114,   115,   116,   117,   118,     0,     0,
       0,     0,    39,     0,     0,     0,     0,     0,     0,     0,
     197,    40,    41,     0,     0,     0,     0,    42,    43,    44,
       0,    45,    46,     0,    47,     4,     5,     0,     6,     7,
       0,     0,     0,     9,     0,     0,     0,    11,     0,     0,
       0,    13,     0,    14,     0,    16,    17,     0,     0,    19,
       0,     0,     0,    21,     0,     0,     0,    22,    23,    24,
      25,    26,     0,    28,    29,    30,    31,    32,    33,    34,
      98,    36,    99,     0,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,     0,     0,     0,     0,    39,     0,     0,
       0,     0,     0,     0,     0,   203,    40,    41,     0,     0,
       0,     0,    42,    43,    44,     0,    45,    46,     0,    47,
       4,     5,     0,     6,     7,     0,     0,     0,     9,     0,
       0,     0,    11,     0,     0,     0,   229,     0,    14,     0,
      16,    17,     0,     0,    19,     0,     0,     0,    21,     0,
       0,     0,    22,    23,    24,    25,    26,     0,    28,    29,
      30,    31,    32,    33,    34,    98,    36,    99,     0,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   109,   110,
     111,   112,   113,   114,   115,   116,   117,   118,     0,     0,
       0,     0,    39,     0,     0,     0,     0,     0,     0,     0,
     234,    40,    41,     0,     0,     0,     0,    42,    43,    44,
       0,    45,    46,    98,    47,    99,     0,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   109,   110,   111,   112,
     113,   114,   115,   116,   117,   118,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    98,   261,    99,
     228,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
      98,     0,    99,     0,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   208,     0,    99,     0,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   109,   110,   111,   112,
     113,   114,   115,   116,   117,   118,    99,     0,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,   112,   113,
     114,   115,   116,   117,   118,   102,   103,   104,   105,   106,
     107,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,   103,   104,   105,   106,   107,   108,   109,   110,
     111,   112,   113,   114,   115,   116,   117,   118,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   105,   106,   107,   108,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,  -127,  -127,   107,
     108,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,  -127,  -127,  -127,  -127,   111,   112,   113,   114,   115,
     116,   117,   118
};

static const yytype_int16 yycheck[] =
{
      48,    10,    18,    12,    40,    41,    42,    43,    15,     0,
      20,    47,    19,    69,    68,   104,    28,    68,    27,    68,
      32,    20,    32,    59,   104,    69,    62,    98,    64,    28,
      66,    67,    69,    32,    98,    98,    72,    73,    74,    95,
      96,    77,    78,    79,    80,   188,    82,   103,   191,   103,
     104,   194,   103,   102,    98,    98,    65,    28,    95,    96,
     203,    32,    98,    98,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   106,
     116,   117,   118,   119,   120,   121,   122,    98,   231,   105,
      66,    67,     4,     5,   101,     7,     8,    73,    74,    20,
      12,    15,    78,    79,    16,    19,   249,    20,    20,    98,
      22,    24,    24,    25,    98,    28,    28,   153,    98,    32,
      32,    98,    30,   159,    36,    37,    38,    39,    40,    98,
      42,    43,    44,    45,    46,    47,    48,    98,    50,    98,
      39,    88,    89,    90,   120,   121,    45,    46,   184,     4,
       5,   102,     7,     8,    98,    98,    98,   205,    98,    98,
      20,    16,   103,    20,    76,   201,    65,   102,   102,    68,
     102,    71,   208,    85,    86,   103,    18,    68,   187,    91,
      92,    93,   105,    95,    96,   102,    98,   103,   103,    21,
     100,   103,    85,    86,    87,    88,    89,    90,   103,   106,
     103,   103,   238,   239,     4,     5,   254,     7,     8,   103,
      10,   103,    12,    13,   103,   103,    16,    17,    18,   103,
      20,   230,    22,    23,    24,    25,   102,    27,    28,    29,
     104,   103,    32,    71,   105,   283,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    68,    71,    70,    18,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    76,   103,    68,   103,
      18,   100,   254,   235,   250,    85,    86,    63,   276,   139,
     105,    91,    92,    93,    -1,    95,    96,   246,    98,    -1,
     100,    -1,   102,     4,     5,    -1,     7,     8,    -1,    10,
      -1,    12,    13,    -1,    -1,    16,    17,    -1,    -1,    20,
      -1,    22,    23,    24,    25,    -1,    27,    28,    29,    -1,
      -1,    32,    -1,    -1,    -1,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    68,    -1,    70,    -1,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    76,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    85,    86,    -1,   103,    -1,    -1,
      91,    92,    93,    -1,    95,    96,    -1,    98,    -1,   100,
     101,   102,     4,     5,    -1,     7,     8,    -1,    10,    -1,
      12,    13,    -1,    -1,    16,    17,    -1,    -1,    20,    -1,
      22,    23,    24,    25,    -1,    27,    28,    29,    -1,    -1,
      32,    -1,    -1,    -1,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    68,    -1,    70,    -1,    72,    73,    74,    75,    76,
      77,    78,    79,    80,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    76,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    85,    86,    -1,   103,    -1,    -1,    91,
      92,    93,    -1,    95,    96,    -1,    98,    -1,   100,   101,
     102,     4,     5,    -1,     7,     8,    -1,    10,    -1,    12,
      13,    -1,    -1,    16,    17,    -1,    -1,    20,    -1,    22,
      23,    24,    25,    -1,    27,    28,    29,    -1,    -1,    32,
      -1,    -1,    -1,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      68,    -1,    70,    -1,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    76,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    85,    86,    -1,   103,    -1,    -1,    91,    92,
      93,    -1,    95,    96,    -1,    98,    -1,   100,    -1,   102,
       4,     5,    -1,     7,     8,    -1,    -1,    -1,    12,    -1,
      -1,    -1,    16,    -1,    18,    -1,    20,    -1,    22,    -1,
      24,    25,    -1,    -1,    28,    -1,    -1,    -1,    32,    -1,
      -1,    -1,    36,    37,    38,    39,    40,    -1,    42,    43,
      44,    45,    46,    47,    48,    68,    50,    70,    -1,    72,
      73,    74,    75,    76,    77,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    -1,    -1,
      -1,    -1,    76,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     103,    85,    86,    -1,    -1,    -1,    -1,    91,    92,    93,
      -1,    95,    96,    -1,    98,     4,     5,    -1,     7,     8,
      -1,    -1,    -1,    12,    -1,    -1,    -1,    16,    -1,    -1,
      -1,    20,    -1,    22,    -1,    24,    25,    -1,    -1,    28,
      -1,    -1,    -1,    32,    -1,    -1,    -1,    36,    37,    38,
      39,    40,    -1,    42,    43,    44,    45,    46,    47,    48,
      68,    50,    70,    -1,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    -1,    -1,    -1,    -1,    76,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   103,    85,    86,    -1,    -1,
      -1,    -1,    91,    92,    93,    -1,    95,    96,    -1,    98,
       4,     5,    -1,     7,     8,    -1,    -1,    -1,    12,    -1,
      -1,    -1,    16,    -1,    -1,    -1,    20,    -1,    22,    -1,
      24,    25,    -1,    -1,    28,    -1,    -1,    -1,    32,    -1,
      -1,    -1,    36,    37,    38,    39,    40,    -1,    42,    43,
      44,    45,    46,    47,    48,    68,    50,    70,    -1,    72,
      73,    74,    75,    76,    77,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    -1,    -1,
      -1,    -1,    76,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     103,    85,    86,    -1,    -1,    -1,    -1,    91,    92,    93,
      -1,    95,    96,    68,    98,    70,    -1,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    68,   103,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      68,    -1,    70,    -1,    72,    73,    74,    75,    76,    77,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    68,    -1,    70,    -1,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    70,    -1,    72,    73,
      74,    75,    76,    77,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    73,    74,    75,
      76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    74,    75,    76,    77,    78,
      79,    80,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    75,    76,    77,    78,    79,    80,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    76,    77,
      78,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    77,    78,    79,    80,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    77,    78,    79,
      80,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    79,    80,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,   108,   112,     0,     4,     5,     7,     8,    10,    12,
      13,    16,    17,    20,    22,    23,    24,    25,    27,    28,
      29,    32,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    76,
      85,    86,    91,    92,    93,    95,    96,    98,   100,   109,
     113,   121,   128,   129,   132,   133,   134,   135,   136,   104,
     128,   128,   104,   138,    98,    98,    98,    98,    98,   106,
     145,   146,    98,    98,    98,    98,   128,    98,    98,    98,
      98,    98,    98,   136,   129,   129,   129,   129,    20,   136,
     136,    18,   129,   112,    30,   110,   114,   102,    68,    70,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    98,
      98,    98,    69,    95,    96,   147,   129,   139,   140,   102,
     102,   129,   145,   129,   128,   136,   139,   139,   136,    28,
      32,   129,   139,   139,   103,   102,   103,   129,   139,   139,
     129,    20,   129,   103,   103,   101,    71,    18,   129,   130,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   129,   129,   129,   129,    20,   129,   129,   129,   129,
     139,   139,   129,   105,    68,   105,   103,   102,   103,   103,
     103,   103,   106,   146,   103,   103,   103,   103,   103,   103,
     116,    68,   103,   103,   129,   111,   115,   129,    68,   103,
     103,   129,   128,   113,   113,    28,   113,   103,   129,   113,
     112,    20,    24,    28,    32,   122,   123,   125,    71,    20,
     102,    21,   127,   100,   103,    68,   102,   124,   131,    68,
     103,   128,   113,   117,   123,   104,   137,   129,   129,   103,
      15,    19,   101,   118,   119,   120,    18,   105,    98,   147,
     148,   103,   113,   133,    71,   112,   120,   105,   149,    71,
      18,   141,   142,   143,   144,   150,    68,    20,    32,   126,
     103,    18,   143,   100,   112,   101
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   107,   108,   109,   109,   111,   110,   112,   112,   113,
     113,   114,   115,   113,   113,   113,   113,   113,   113,   113,
     113,   116,   113,   113,   113,   113,   117,   117,   118,   119,
     119,   120,   120,   121,   121,   122,   122,   124,   123,   125,
     125,   125,   125,   126,   126,   127,   127,   128,   128,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     130,   131,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   132,   132,   133,   133,   133,   133,   133,   134,   134,
     134,   135,   135,   136,   136,   136,   137,   137,   137,   138,
     138,   139,   139,   140,   140,   141,   141,   141,   142,   142,
     144,   143,   145,   145,   145,   146,   146,   147,   147,   148,
     149,   150,   148
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     0,     2,     0,     4,     0,     2,     3,
       2,     0,     0,     6,     6,     5,     9,     5,     4,     6,
       5,     0,     8,     3,     3,     3,     0,     2,     2,     1,
       2,     3,     2,     0,     1,     1,     3,     0,     4,     1,
       1,     1,     1,     1,     1,     0,     2,     0,     1,     3,
       4,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       0,     0,     7,     2,     2,     2,     2,     2,     2,     4,
       4,     4,     4,     4,     6,     8,     4,     3,     3,     4,
       4,     4,     4,     4,     2,     2,     2,     3,     2,     2,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     2,     3,     1,     0,     2,     3,     0,
       3,     0,     1,     1,     3,     0,     1,     1,     1,     3,
       0,     3,     0,     1,     3,     2,     2,     0,     2,     1,
       0,     0,     8
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
#if EX_DEBUG

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
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
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
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
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
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[+yyssp[yyi + 1 - yynrhs]],
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
#else /* !EX_DEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !EX_DEBUG */


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
#   define yystrlen(S) (YY_CAST (YYPTRDIFF_T, strlen (S)))
#  else
/* Return the length of YYSTR.  */
static YYPTRDIFF_T
yystrlen (const char *yystr)
{
  YYPTRDIFF_T yylen;
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
static YYPTRDIFF_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYPTRDIFF_T yyn = 0;
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

  if (yyres)
    return yystpcpy (yyres, yystr) - yyres;
  else
    return yystrlen (yystr);
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
yysyntax_error (YYPTRDIFF_T *yymsg_alloc, char **yymsg,
                yy_state_t *yyssp, int yytoken)
{
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat: reported tokens (one for the "unexpected",
     one per "expected"). */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Actual size of YYARG. */
  int yycount = 0;
  /* Cumulated lengths of YYARG.  */
  YYPTRDIFF_T yysize = 0;

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
      int yyn = yypact[+*yyssp];
      YYPTRDIFF_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
      yysize = yysize0;
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
                  YYPTRDIFF_T yysize1
                    = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
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
    /* Don't count the "%s"s in the final size, but reserve room for
       the terminator.  */
    YYPTRDIFF_T yysize1 = yysize + (yystrlen (yyformat) - 2 * yycount) + 1;
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
          ++yyp;
          ++yyformat;
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
    yy_state_fast_t yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss;
    yy_state_t *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYPTRDIFF_T yystacksize;

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
  YYPTRDIFF_T yymsg_alloc = sizeof yymsgbuf;
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
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
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
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
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

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

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
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
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
  case 2:
#line 180 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[-1].expr) && !(expr.program->disc->flags & EX_STRICT))
			{
				if (expr.program->main.value && !(expr.program->disc->flags & EX_RETAIN))
					exfreenode(expr.program, expr.program->main.value);
				if ((yyvsp[-1].expr)->op == S2B)
				{
					Exnode_t*	x;

					x = (yyvsp[-1].expr);
					(yyvsp[-1].expr) = x->data.operand.left;
					x->data.operand.left = 0;
					exfreenode(expr.program, x);
				}
				expr.program->main.lex = PROCEDURE;
				expr.program->main.value = exnewnode(expr.program, PROCEDURE, 1, (yyvsp[-1].expr)->type, NULL, (yyvsp[-1].expr));
			}
		}
#line 1911 "exparse.c"
    break;

  case 5:
#line 204 "../../lib/expr/exparse.y"
                                  {
				Dtdisc_t*	disc;

				if (expr.procedure)
					exerror("no nested function definitions");
				(yyvsp[-1].id)->lex = PROCEDURE;
				expr.procedure = (yyvsp[-1].id)->value = exnewnode(expr.program, PROCEDURE, 1, (yyvsp[-1].id)->type, NULL, NULL);
				expr.procedure->type = INTEGER;
				if (!(disc = calloc(1, sizeof(Dtdisc_t))))
					exnospace();
				disc->key = offsetof(Exid_t, name);
				if (expr.assigned && !streq((yyvsp[-1].id)->name, "begin"))
				{
					if (!(expr.procedure->data.procedure.frame = dtopen(disc, Dtset)) || !dtview(expr.procedure->data.procedure.frame, expr.program->symbols))
						exnospace();
					expr.program->symbols = expr.program->frame = expr.procedure->data.procedure.frame;
				}
			}
#line 1934 "exparse.c"
    break;

  case 6:
#line 222 "../../lib/expr/exparse.y"
                {
			expr.procedure = 0;
			if (expr.program->frame)
			{
				expr.program->symbols = expr.program->frame->view;
				dtview(expr.program->frame, NULL);
				expr.program->frame = 0;
			}
			if ((yyvsp[0].expr) && (yyvsp[0].expr)->op == S2B)
			{
				Exnode_t*	x;

				x = (yyvsp[0].expr);
				(yyvsp[0].expr) = x->data.operand.left;
				x->data.operand.left = 0;
				exfreenode(expr.program, x);
			}
			(yyvsp[-3].id)->value->data.operand.right = excast(expr.program, (yyvsp[0].expr), (yyvsp[-3].id)->type, NULL, 0);
		}
#line 1958 "exparse.c"
    break;

  case 7:
#line 244 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 1966 "exparse.c"
    break;

  case 8:
#line 248 "../../lib/expr/exparse.y"
                {
			if (!(yyvsp[-1].expr))
				(yyval.expr) = (yyvsp[0].expr);
			else if (!(yyvsp[0].expr))
				(yyval.expr) = (yyvsp[-1].expr);
			else if ((yyvsp[-1].expr)->op == CONSTANT)
			{
				exfreenode(expr.program, (yyvsp[-1].expr));
				(yyval.expr) = (yyvsp[0].expr);
			}
			else (yyval.expr) = exnewnode(expr.program, ';', 1, (yyvsp[0].expr)->type, (yyvsp[-1].expr), (yyvsp[0].expr));
		}
#line 1983 "exparse.c"
    break;

  case 9:
#line 263 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[-1].expr);
		}
#line 1991 "exparse.c"
    break;

  case 10:
#line 267 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = ((yyvsp[-1].expr) && (yyvsp[-1].expr)->type == STRING) ? exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[-1].expr), NULL) : (yyvsp[-1].expr);
		}
#line 1999 "exparse.c"
    break;

  case 11:
#line 270 "../../lib/expr/exparse.y"
                               {expr.instatic=(yyvsp[0].integer);}
#line 2005 "exparse.c"
    break;

  case 12:
#line 270 "../../lib/expr/exparse.y"
                                                           {expr.declare=(yyvsp[0].id)->type;}
#line 2011 "exparse.c"
    break;

  case 13:
#line 271 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[-1].expr);
			expr.declare = 0;
		}
#line 2020 "exparse.c"
    break;

  case 14:
#line 276 "../../lib/expr/exparse.y"
                {
			if (exisAssign ((yyvsp[-3].expr)))
				exwarn ("assignment used as boolean in if statement");
			if ((yyvsp[-3].expr)->type == STRING)
				(yyvsp[-3].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[-3].expr), NULL);
			else if (!INTEGRAL((yyvsp[-3].expr)->type))
				(yyvsp[-3].expr) = excast(expr.program, (yyvsp[-3].expr), INTEGER, NULL, 0);
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-5].id)->index, 1, INTEGER, (yyvsp[-3].expr), exnewnode(expr.program, ':', 1, (yyvsp[-1].expr) ? (yyvsp[-1].expr)->type : 0, (yyvsp[-1].expr), (yyvsp[0].expr)));
		}
#line 2034 "exparse.c"
    break;

  case 15:
#line 286 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ITERATE, 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.generate.array = (yyvsp[-2].expr);
			if (!(yyvsp[-2].expr)->data.variable.index || (yyvsp[-2].expr)->data.variable.index->op != DYNAMIC)
				exerror("simple index variable expected");
			(yyval.expr)->data.generate.index = (yyvsp[-2].expr)->data.variable.index->data.variable.symbol;
			if ((yyvsp[-2].expr)->op == ID && (yyval.expr)->data.generate.index->type != INTEGER)
				exerror("integer index variable expected");
			exfreenode(expr.program, (yyvsp[-2].expr)->data.variable.index);
			(yyvsp[-2].expr)->data.variable.index = 0;
			(yyval.expr)->data.generate.statement = (yyvsp[0].expr);
		}
#line 2051 "exparse.c"
    break;

  case 16:
#line 299 "../../lib/expr/exparse.y"
                {
			if (!(yyvsp[-4].expr))
			{
				(yyvsp[-4].expr) = exnewnode(expr.program, CONSTANT, 0, INTEGER, NULL, NULL);
				(yyvsp[-4].expr)->data.constant.value.integer = 1;
			}
			else if ((yyvsp[-4].expr)->type == STRING)
				(yyvsp[-4].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[-4].expr), NULL);
			else if (!INTEGRAL((yyvsp[-4].expr)->type))
				(yyvsp[-4].expr) = excast(expr.program, (yyvsp[-4].expr), INTEGER, NULL, 0);
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-8].id)->index, 1, INTEGER, (yyvsp[-4].expr), exnewnode(expr.program, ';', 1, 0, (yyvsp[-2].expr), (yyvsp[0].expr)));
			if ((yyvsp[-6].expr))
				(yyval.expr) = exnewnode(expr.program, ';', 1, INTEGER, (yyvsp[-6].expr), (yyval.expr));
		}
#line 2070 "exparse.c"
    break;

  case 17:
#line 314 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ITERATER, 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.generate.array = (yyvsp[-2].expr);
			if (!(yyvsp[-2].expr)->data.variable.index || (yyvsp[-2].expr)->data.variable.index->op != DYNAMIC)
				exerror("simple index variable expected");
			(yyval.expr)->data.generate.index = (yyvsp[-2].expr)->data.variable.index->data.variable.symbol;
			if ((yyvsp[-2].expr)->op == ID && (yyval.expr)->data.generate.index->type != INTEGER)
				exerror("integer index variable expected");
			exfreenode(expr.program, (yyvsp[-2].expr)->data.variable.index);
			(yyvsp[-2].expr)->data.variable.index = 0;
			(yyval.expr)->data.generate.statement = (yyvsp[0].expr);
		}
#line 2087 "exparse.c"
    break;

  case 18:
#line 327 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[-1].id)->local.pointer == 0)
              			exerror("cannot apply unset to non-array %s", (yyvsp[-1].id)->name);
			(yyval.expr) = exnewnode(expr.program, UNSET, 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.variable.symbol = (yyvsp[-1].id);
			(yyval.expr)->data.variable.index = NULL;
		}
#line 2099 "exparse.c"
    break;

  case 19:
#line 335 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[-3].id)->local.pointer == 0)
              			exerror("cannot apply unset to non-array %s", (yyvsp[-3].id)->name);
			if (((yyvsp[-3].id)->index_type > 0) && ((yyvsp[-1].expr)->type != (yyvsp[-3].id)->index_type))
            		    exerror("%s indices must have type %s, not %s", 
				(yyvsp[-3].id)->name, extypename(expr.program, (yyvsp[-3].id)->index_type),extypename(expr.program, (yyvsp[-1].expr)->type));
			(yyval.expr) = exnewnode(expr.program, UNSET, 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.variable.symbol = (yyvsp[-3].id);
			(yyval.expr)->data.variable.index = (yyvsp[-1].expr);
		}
#line 2114 "exparse.c"
    break;

  case 20:
#line 346 "../../lib/expr/exparse.y"
                {
			if (exisAssign ((yyvsp[-2].expr)))
				exwarn ("assignment used as boolean in while statement");
			if ((yyvsp[-2].expr)->type == STRING)
				(yyvsp[-2].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[-2].expr), NULL);
			else if (!INTEGRAL((yyvsp[-2].expr)->type))
				(yyvsp[-2].expr) = excast(expr.program, (yyvsp[-2].expr), INTEGER, NULL, 0);
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-4].id)->index, 1, INTEGER, (yyvsp[-2].expr), exnewnode(expr.program, ';', 1, 0, NULL, (yyvsp[0].expr)));
		}
#line 2128 "exparse.c"
    break;

  case 21:
#line 355 "../../lib/expr/exparse.y"
                                        {expr.declare=(yyvsp[0].expr)->type;}
#line 2134 "exparse.c"
    break;

  case 22:
#line 356 "../../lib/expr/exparse.y"
                {
			Switch_t*	sw = expr.swstate;

			(yyval.expr) = exnewnode(expr.program, (yyvsp[-7].id)->index, 1, INTEGER, (yyvsp[-5].expr), exnewnode(expr.program, DEFAULT, 1, 0, sw->defcase, sw->firstcase));
			expr.swstate = expr.swstate->prev;
			free(sw->base);
			if (sw != &swstate)
				free(sw);
			expr.declare = 0;
		}
#line 2149 "exparse.c"
    break;

  case 23:
#line 367 "../../lib/expr/exparse.y"
                {
		loopop:
			if (!(yyvsp[-1].expr))
			{
				(yyvsp[-1].expr) = exnewnode(expr.program, CONSTANT, 0, INTEGER, NULL, NULL);
				(yyvsp[-1].expr)->data.constant.value.integer = 1;
			}
			else if (!INTEGRAL((yyvsp[-1].expr)->type))
				(yyvsp[-1].expr) = excast(expr.program, (yyvsp[-1].expr), INTEGER, NULL, 0);
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-2].id)->index, 1, INTEGER, (yyvsp[-1].expr), NULL);
		}
#line 2165 "exparse.c"
    break;

  case 24:
#line 379 "../../lib/expr/exparse.y"
                {
			goto loopop;
		}
#line 2173 "exparse.c"
    break;

  case 25:
#line 383 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[-1].expr))
			{
				if (expr.procedure && !expr.procedure->type)
					exerror("return in void function");
				(yyvsp[-1].expr) = excast(expr.program, (yyvsp[-1].expr), expr.procedure ? expr.procedure->type : INTEGER, NULL, 0);
			}
			(yyval.expr) = exnewnode(expr.program, RETURN, 1, (yyvsp[-1].expr) ? (yyvsp[-1].expr)->type : 0, (yyvsp[-1].expr), NULL);
		}
#line 2187 "exparse.c"
    break;

  case 26:
#line 395 "../../lib/expr/exparse.y"
                {
			Switch_t*		sw;

			if (expr.swstate)
			{
				if (!(sw = calloc(1, sizeof(Switch_t))))
				{
					exnospace();
					sw = &swstate;
				}
				sw->prev = expr.swstate;
			}
			else
				sw = &swstate;
			expr.swstate = sw;
			sw->type = expr.declare;
			sw->firstcase = 0;
			sw->lastcase = 0;
			sw->defcase = 0;
			sw->def = 0;
			size_t n = 8;
			if (!(sw->base = calloc(n, sizeof(Extype_t*))))
			{
				exnospace();
				n = 0;
			}
			sw->cur = sw->base;
			sw->last = sw->base + n;
		}
#line 2221 "exparse.c"
    break;

  case 28:
#line 428 "../../lib/expr/exparse.y"
                {
			Switch_t*	sw = expr.swstate;

			(yyval.expr) = exnewnode(expr.program, CASE, 1, 0, (yyvsp[0].expr), NULL);
			if (sw->cur > sw->base)
			{
				if (sw->lastcase)
					sw->lastcase->data.select.next = (yyval.expr);
				else
					sw->firstcase = (yyval.expr);
				sw->lastcase = (yyval.expr);
				size_t n = (size_t)(sw->cur - sw->base);
				sw->cur = sw->base;
				(yyval.expr)->data.select.constant = exalloc(expr.program, (n + 1) * sizeof(Extype_t*));
				memcpy((yyval.expr)->data.select.constant, sw->base, n * sizeof(Extype_t*));
				(yyval.expr)->data.select.constant[n] = 0;
			}
			else
				(yyval.expr)->data.select.constant = 0;
			if (sw->def)
			{
				sw->def = 0;
				if (sw->defcase)
					exerror("duplicate default in switch");
				else
					sw->defcase = (yyvsp[0].expr);
			}
		}
#line 2254 "exparse.c"
    break;

  case 31:
#line 463 "../../lib/expr/exparse.y"
                {
			if (expr.swstate->cur >= expr.swstate->last)
			{
				size_t n = (size_t)(expr.swstate->cur - expr.swstate->base);
				if (!(expr.swstate->base = realloc(expr.swstate->base, sizeof(Extype_t*) * 2 * n)))
				{
					exerror("too many case labels for switch");
					n = 0;
				}
				expr.swstate->cur = expr.swstate->base + n;
				expr.swstate->last = expr.swstate->base + 2 * n;
			}
			if (expr.swstate->cur)
			{
				(yyvsp[-1].expr) = excast(expr.program, (yyvsp[-1].expr), expr.swstate->type, NULL, 0);
				*expr.swstate->cur++ = &((yyvsp[-1].expr)->data.constant.value);
			}
		}
#line 2277 "exparse.c"
    break;

  case 32:
#line 482 "../../lib/expr/exparse.y"
                {
			expr.swstate->def = 1;
		}
#line 2285 "exparse.c"
    break;

  case 33:
#line 488 "../../lib/expr/exparse.y"
                {
			(yyval.integer) = 0;
		}
#line 2293 "exparse.c"
    break;

  case 34:
#line 492 "../../lib/expr/exparse.y"
                {
			(yyval.integer) = 1;
		}
#line 2301 "exparse.c"
    break;

  case 36:
#line 499 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[0].expr))
				(yyval.expr) = (yyvsp[-2].expr) ? exnewnode(expr.program, ',', 1, (yyvsp[0].expr)->type, (yyvsp[-2].expr), (yyvsp[0].expr)) : (yyvsp[0].expr);
		}
#line 2310 "exparse.c"
    break;

  case 37:
#line 505 "../../lib/expr/exparse.y"
                                 {checkName ((yyvsp[0].id)); expr.id=(yyvsp[0].id);}
#line 2316 "exparse.c"
    break;

  case 38:
#line 506 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
			if (!(yyvsp[-3].id)->type || expr.declare)
				(yyvsp[-3].id)->type = expr.declare;
			if ((yyvsp[0].expr) && (yyvsp[0].expr)->op == PROCEDURE)
			{
				(yyvsp[-3].id)->lex = PROCEDURE;
				(yyvsp[-3].id)->type = (yyvsp[0].expr)->type;
				(yyvsp[-3].id)->value = (yyvsp[0].expr);
			}
			else
			{
				(yyvsp[-3].id)->lex = DYNAMIC;
				(yyvsp[-3].id)->value = exnewnode(expr.program, 0, 0, 0, NULL, NULL);
				if ((yyvsp[-1].integer) && !(yyvsp[-3].id)->local.pointer)
				{
					Dtdisc_t*	disc;

					if (!(disc = calloc(1, sizeof(Dtdisc_t))))
						exnospace();
					if ((yyvsp[-1].integer) == INTEGER) {
						disc->key = offsetof(Exassoc_t, key);
						disc->size = sizeof(Extype_t);
						disc->comparf = (Dtcompar_f)cmpKey;
					}
					else
						disc->key = offsetof(Exassoc_t, name);
					if (!((yyvsp[-3].id)->local.pointer = (char*)dtopen(disc, Dtoset)))
						exerror("%s: cannot initialize associative array", (yyvsp[-3].id)->name);
					(yyvsp[-3].id)->index_type = (yyvsp[-1].integer); /* -1 indicates no typechecking */
				}
				if ((yyvsp[0].expr))
				{
					if ((yyvsp[0].expr)->type != (yyvsp[-3].id)->type)
					{
						(yyvsp[0].expr)->type = (yyvsp[-3].id)->type;
						(yyvsp[0].expr)->data.operand.right = excast(expr.program, (yyvsp[0].expr)->data.operand.right, (yyvsp[-3].id)->type, NULL, 0);
					}
					(yyvsp[0].expr)->data.operand.left = exnewnode(expr.program, DYNAMIC, 0, (yyvsp[-3].id)->type, NULL, NULL);
					(yyvsp[0].expr)->data.operand.left->data.variable.symbol = (yyvsp[-3].id);
					(yyval.expr) = (yyvsp[0].expr);
				}
				else if (!(yyvsp[-1].integer))
					(yyvsp[-3].id)->value->data.value = exzero((yyvsp[-3].id)->type);
			}
		}
#line 2367 "exparse.c"
    break;

  case 45:
#line 565 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 2375 "exparse.c"
    break;

  case 46:
#line 569 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[0].expr);
		}
#line 2383 "exparse.c"
    break;

  case 47:
#line 575 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 2391 "exparse.c"
    break;

  case 49:
#line 582 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[-1].expr);
		}
#line 2399 "exparse.c"
    break;

  case 50:
#line 586 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = ((yyvsp[0].expr)->type == (yyvsp[-2].id)->type) ? (yyvsp[0].expr) : excast(expr.program, (yyvsp[0].expr), (yyvsp[-2].id)->type, NULL, 0);
		}
#line 2407 "exparse.c"
    break;

  case 51:
#line 590 "../../lib/expr/exparse.y"
                {
			int	rel;

		relational:
			rel = INTEGER;
			goto coerce;
		binary:
			rel = 0;
		coerce:
			if (!(yyvsp[-2].expr)->type)
			{
				if (!(yyvsp[0].expr)->type)
					(yyvsp[-2].expr)->type = (yyvsp[0].expr)->type = rel ? STRING : INTEGER;
				else
					(yyvsp[-2].expr)->type = (yyvsp[0].expr)->type;
			}
			else if (!(yyvsp[0].expr)->type)
				(yyvsp[0].expr)->type = (yyvsp[-2].expr)->type;
			if ((yyvsp[-2].expr)->type != (yyvsp[0].expr)->type)
			{
				if ((yyvsp[-2].expr)->type == STRING)
					(yyvsp[-2].expr) = excast(expr.program, (yyvsp[-2].expr), (yyvsp[0].expr)->type, (yyvsp[0].expr), 0);
				else if ((yyvsp[0].expr)->type == STRING)
					(yyvsp[0].expr) = excast(expr.program, (yyvsp[0].expr), (yyvsp[-2].expr)->type, (yyvsp[-2].expr), 0);
				else if ((yyvsp[-2].expr)->type == FLOATING)
					(yyvsp[0].expr) = excast(expr.program, (yyvsp[0].expr), FLOATING, (yyvsp[-2].expr), 0);
				else if ((yyvsp[0].expr)->type == FLOATING)
					(yyvsp[-2].expr) = excast(expr.program, (yyvsp[-2].expr), FLOATING, (yyvsp[0].expr), 0);
			}
			if (!rel)
				rel = ((yyvsp[-2].expr)->type == STRING) ? STRING : (((yyvsp[-2].expr)->type == UNSIGNED) ? UNSIGNED : (yyvsp[0].expr)->type);
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-1].op), 1, rel, (yyvsp[-2].expr), (yyvsp[0].expr));
			if (!expr.program->errors && (yyvsp[-2].expr)->op == CONSTANT && (yyvsp[0].expr)->op == CONSTANT)
			{
				(yyval.expr)->data.constant.value = exeval(expr.program, (yyval.expr), NULL);
				/* If a constant string, re-allocate from program heap. This is because the
				 * value was constructed from string operators, which create a value in the 
				 * temporary heap, which is cleared when exeval is called again. 
				 */
				if ((yyval.expr)->type == STRING) {
					(yyval.expr)->data.constant.value.string =
						vmstrdup(expr.program->vm, (yyval.expr)->data.constant.value.string);
				}
				(yyval.expr)->binary = 0;
				(yyval.expr)->op = CONSTANT;
				exfreenode(expr.program, (yyvsp[-2].expr));
				exfreenode(expr.program, (yyvsp[0].expr));
			}
			else if (!BUILTIN((yyvsp[-2].expr)->type) || !BUILTIN((yyvsp[0].expr)->type)) {
				checkBinary(expr.program, (yyvsp[-2].expr), (yyval.expr), (yyvsp[0].expr));
			}
		}
#line 2464 "exparse.c"
    break;

  case 52:
#line 643 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2472 "exparse.c"
    break;

  case 53:
#line 647 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2480 "exparse.c"
    break;

  case 54:
#line 651 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2488 "exparse.c"
    break;

  case 55:
#line 655 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2496 "exparse.c"
    break;

  case 56:
#line 659 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2504 "exparse.c"
    break;

  case 57:
#line 663 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2512 "exparse.c"
    break;

  case 58:
#line 667 "../../lib/expr/exparse.y"
                {
			goto relational;
		}
#line 2520 "exparse.c"
    break;

  case 59:
#line 671 "../../lib/expr/exparse.y"
                {
			goto relational;
		}
#line 2528 "exparse.c"
    break;

  case 60:
#line 675 "../../lib/expr/exparse.y"
                {
			goto relational;
		}
#line 2536 "exparse.c"
    break;

  case 61:
#line 679 "../../lib/expr/exparse.y"
                {
			goto relational;
		}
#line 2544 "exparse.c"
    break;

  case 62:
#line 683 "../../lib/expr/exparse.y"
                {
			goto relational;
		}
#line 2552 "exparse.c"
    break;

  case 63:
#line 687 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2560 "exparse.c"
    break;

  case 64:
#line 691 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2568 "exparse.c"
    break;

  case 65:
#line 695 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2576 "exparse.c"
    break;

  case 66:
#line 699 "../../lib/expr/exparse.y"
                {
			goto binary;
		}
#line 2584 "exparse.c"
    break;

  case 67:
#line 703 "../../lib/expr/exparse.y"
                {
		logical:
			if ((yyvsp[-2].expr)->type == STRING)
				(yyvsp[-2].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[-2].expr), NULL);
			else if (!BUILTIN((yyvsp[-2].expr)->type))
				(yyvsp[-2].expr) = excast(expr.program, (yyvsp[-2].expr), INTEGER, NULL, 0);
			if ((yyvsp[0].expr)->type == STRING)
				(yyvsp[0].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[0].expr), NULL);
			else if (!BUILTIN((yyvsp[0].expr)->type))
				(yyvsp[0].expr) = excast(expr.program, (yyvsp[0].expr), INTEGER, NULL, 0);
			goto binary;
		}
#line 2601 "exparse.c"
    break;

  case 68:
#line 716 "../../lib/expr/exparse.y"
                {
			goto logical;
		}
#line 2609 "exparse.c"
    break;

  case 69:
#line 720 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[-2].expr)->op == CONSTANT)
			{
				exfreenode(expr.program, (yyvsp[-2].expr));
				(yyval.expr) = (yyvsp[0].expr);
			}
			else
				(yyval.expr) = exnewnode(expr.program, ',', 1, (yyvsp[0].expr)->type, (yyvsp[-2].expr), (yyvsp[0].expr));
		}
#line 2623 "exparse.c"
    break;

  case 70:
#line 729 "../../lib/expr/exparse.y"
                                 {expr.nolabel=1;}
#line 2629 "exparse.c"
    break;

  case 71:
#line 729 "../../lib/expr/exparse.y"
                                                            {expr.nolabel=0;}
#line 2635 "exparse.c"
    break;

  case 72:
#line 730 "../../lib/expr/exparse.y"
                {
			if (!(yyvsp[-3].expr)->type)
			{
				if (!(yyvsp[0].expr)->type)
					(yyvsp[-3].expr)->type = (yyvsp[0].expr)->type = INTEGER;
				else
					(yyvsp[-3].expr)->type = (yyvsp[0].expr)->type;
			}
			else if (!(yyvsp[0].expr)->type)
				(yyvsp[0].expr)->type = (yyvsp[-3].expr)->type;
			if ((yyvsp[-6].expr)->type == STRING)
				(yyvsp[-6].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[-6].expr), NULL);
			else if (!INTEGRAL((yyvsp[-6].expr)->type))
				(yyvsp[-6].expr) = excast(expr.program, (yyvsp[-6].expr), INTEGER, NULL, 0);
			if ((yyvsp[-3].expr)->type != (yyvsp[0].expr)->type)
			{
				if ((yyvsp[-3].expr)->type == STRING || (yyvsp[0].expr)->type == STRING)
					exerror("if statement string type mismatch");
				else if ((yyvsp[-3].expr)->type == FLOATING)
					(yyvsp[0].expr) = excast(expr.program, (yyvsp[0].expr), FLOATING, NULL, 0);
				else if ((yyvsp[0].expr)->type == FLOATING)
					(yyvsp[-3].expr) = excast(expr.program, (yyvsp[-3].expr), FLOATING, NULL, 0);
			}
			if ((yyvsp[-6].expr)->op == CONSTANT)
			{
				if ((yyvsp[-6].expr)->data.constant.value.integer)
				{
					(yyval.expr) = (yyvsp[-3].expr);
					exfreenode(expr.program, (yyvsp[0].expr));
				}
				else
				{
					(yyval.expr) = (yyvsp[0].expr);
					exfreenode(expr.program, (yyvsp[-3].expr));
				}
				exfreenode(expr.program, (yyvsp[-6].expr));
			}
			else
				(yyval.expr) = exnewnode(expr.program, '?', 1, (yyvsp[-3].expr)->type, (yyvsp[-6].expr), exnewnode(expr.program, ':', 1, (yyvsp[-3].expr)->type, (yyvsp[-3].expr), (yyvsp[0].expr)));
		}
#line 2680 "exparse.c"
    break;

  case 73:
#line 771 "../../lib/expr/exparse.y"
                {
		iunary:
			if ((yyvsp[0].expr)->type == STRING)
				(yyvsp[0].expr) = exnewnode(expr.program, S2B, 1, INTEGER, (yyvsp[0].expr), NULL);
			else if (!INTEGRAL((yyvsp[0].expr)->type))
				(yyvsp[0].expr) = excast(expr.program, (yyvsp[0].expr), INTEGER, NULL, 0);
		unary:
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-1].op), 1, (yyvsp[0].expr)->type == UNSIGNED ? INTEGER : (yyvsp[0].expr)->type, (yyvsp[0].expr), NULL);
			if ((yyvsp[0].expr)->op == CONSTANT)
			{
				(yyval.expr)->data.constant.value = exeval(expr.program, (yyval.expr), NULL);
				(yyval.expr)->binary = 0;
				(yyval.expr)->op = CONSTANT;
				exfreenode(expr.program, (yyvsp[0].expr));
			}
			else if (!BUILTIN((yyvsp[0].expr)->type)) {
				checkBinary(expr.program, (yyvsp[0].expr), (yyval.expr), 0);
			}
		}
#line 2704 "exparse.c"
    break;

  case 74:
#line 791 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[0].id)->local.pointer == 0)
              			exerror("cannot apply '#' operator to non-array %s", (yyvsp[0].id)->name);
			(yyval.expr) = exnewnode(expr.program, '#', 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.variable.symbol = (yyvsp[0].id);
		}
#line 2715 "exparse.c"
    break;

  case 75:
#line 798 "../../lib/expr/exparse.y"
                {
			goto iunary;
		}
#line 2723 "exparse.c"
    break;

  case 76:
#line 802 "../../lib/expr/exparse.y"
                {
			goto unary;
		}
#line 2731 "exparse.c"
    break;

  case 77:
#line 806 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[0].expr);
		}
#line 2739 "exparse.c"
    break;

  case 78:
#line 810 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ADDRESS, 0, T((yyvsp[0].expr)->type), (yyvsp[0].expr), NULL);
		}
#line 2747 "exparse.c"
    break;

  case 79:
#line 814 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ARRAY, 1, T((yyvsp[-3].id)->type), call(0, (yyvsp[-3].id), (yyvsp[-1].expr)), (yyvsp[-1].expr));
		}
#line 2755 "exparse.c"
    break;

  case 80:
#line 818 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, FUNCTION, 1, T((yyvsp[-3].id)->type), call(0, (yyvsp[-3].id), (yyvsp[-1].expr)), (yyvsp[-1].expr));
		}
#line 2763 "exparse.c"
    break;

  case 81:
#line 822 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewsub (expr.program, (yyvsp[-1].expr), GSUB);
		}
#line 2771 "exparse.c"
    break;

  case 82:
#line 826 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewsub (expr.program, (yyvsp[-1].expr), SUB);
		}
#line 2779 "exparse.c"
    break;

  case 83:
#line 830 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewsubstr (expr.program, (yyvsp[-1].expr));
		}
#line 2787 "exparse.c"
    break;

  case 84:
#line 834 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewsplit (expr.program, (yyvsp[-5].id)->index, (yyvsp[-1].id), (yyvsp[-3].expr), NULL);
		}
#line 2795 "exparse.c"
    break;

  case 85:
#line 838 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewsplit (expr.program, (yyvsp[-7].id)->index, (yyvsp[-3].id), (yyvsp[-5].expr), (yyvsp[-1].expr));
		}
#line 2803 "exparse.c"
    break;

  case 86:
#line 842 "../../lib/expr/exparse.y"
                {
			if (!INTEGRAL((yyvsp[-1].expr)->type))
				(yyvsp[-1].expr) = excast(expr.program, (yyvsp[-1].expr), INTEGER, NULL, 0);
			(yyval.expr) = exnewnode(expr.program, EXIT, 1, INTEGER, (yyvsp[-1].expr), NULL);
		}
#line 2813 "exparse.c"
    break;

  case 87:
#line 848 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, RAND, 0, FLOATING, NULL, NULL);
		}
#line 2821 "exparse.c"
    break;

  case 88:
#line 852 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, SRAND, 0, INTEGER, NULL, NULL);
		}
#line 2829 "exparse.c"
    break;

  case 89:
#line 856 "../../lib/expr/exparse.y"
                {
			if (!INTEGRAL((yyvsp[-1].expr)->type))
				(yyvsp[-1].expr) = excast(expr.program, (yyvsp[-1].expr), INTEGER, NULL, 0);
			(yyval.expr) = exnewnode(expr.program, SRAND, 1, INTEGER, (yyvsp[-1].expr), NULL);
		}
#line 2839 "exparse.c"
    break;

  case 90:
#line 862 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, CALL, 1, (yyvsp[-3].id)->type, NULL, (yyvsp[-1].expr));
			(yyval.expr)->data.call.procedure = (yyvsp[-3].id);
		}
#line 2848 "exparse.c"
    break;

  case 91:
#line 867 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exprint(expr.program, (yyvsp[-3].id), (yyvsp[-1].expr));
		}
#line 2856 "exparse.c"
    break;

  case 92:
#line 871 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-3].id)->index, 0, (yyvsp[-3].id)->type, NULL, NULL);
			if ((yyvsp[-1].expr) && (yyvsp[-1].expr)->data.operand.left->type == INTEGER)
			{
				(yyval.expr)->data.print.descriptor = (yyvsp[-1].expr)->data.operand.left;
				(yyvsp[-1].expr) = (yyvsp[-1].expr)->data.operand.right;
			}
			else 
				switch ((yyvsp[-3].id)->index)
				{
				case QUERY:
					(yyval.expr)->data.print.descriptor = exnewnode(expr.program, CONSTANT, 0, INTEGER, NULL, NULL);
					(yyval.expr)->data.print.descriptor->data.constant.value.integer = 2;
					break;
				case PRINTF:
					(yyval.expr)->data.print.descriptor = exnewnode(expr.program, CONSTANT, 0, INTEGER, NULL, NULL);
					(yyval.expr)->data.print.descriptor->data.constant.value.integer = 1;
					break;
				case SPRINTF:
					(yyval.expr)->data.print.descriptor = 0;
					break;
				}
			(yyval.expr)->data.print.args = preprint((yyvsp[-1].expr));
		}
#line 2885 "exparse.c"
    break;

  case 93:
#line 896 "../../lib/expr/exparse.y"
                {
			Exnode_t*	x;

			(yyval.expr) = exnewnode(expr.program, (yyvsp[-3].id)->index, 0, (yyvsp[-3].id)->type, NULL, NULL);
			if ((yyvsp[-1].expr) && (yyvsp[-1].expr)->data.operand.left->type == INTEGER)
			{
				(yyval.expr)->data.scan.descriptor = (yyvsp[-1].expr)->data.operand.left;
				(yyvsp[-1].expr) = (yyvsp[-1].expr)->data.operand.right;
			}
			else 
				switch ((yyvsp[-3].id)->index)
				{
				case SCANF:
					(yyval.expr)->data.scan.descriptor = 0;
					break;
				case SSCANF:
					if ((yyvsp[-1].expr) && (yyvsp[-1].expr)->data.operand.left->type == STRING)
					{
						(yyval.expr)->data.scan.descriptor = (yyvsp[-1].expr)->data.operand.left;
						(yyvsp[-1].expr) = (yyvsp[-1].expr)->data.operand.right;
					}
					else
						exerror("%s: string argument expected", (yyvsp[-3].id)->name);
					break;
				}
			if (!(yyvsp[-1].expr) || !(yyvsp[-1].expr)->data.operand.left || (yyvsp[-1].expr)->data.operand.left->type != STRING)
				exerror("%s: format argument expected", (yyvsp[-3].id)->name);
			(yyval.expr)->data.scan.format = (yyvsp[-1].expr)->data.operand.left;
			for (x = (yyval.expr)->data.scan.args = (yyvsp[-1].expr)->data.operand.right; x; x = x->data.operand.right)
			{
				if (x->data.operand.left->op != ADDRESS)
					exerror("%s: address argument expected", (yyvsp[-3].id)->name);
				x->data.operand.left = x->data.operand.left->data.operand.left;
			}
		}
#line 2925 "exparse.c"
    break;

  case 94:
#line 932 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[0].expr))
			{
				if ((yyvsp[-1].expr)->op == ID && !expr.program->disc->setf)
					exerror("%s: variable assignment not supported", (yyvsp[-1].expr)->data.variable.symbol->name);
				else
				{
					if (!(yyvsp[-1].expr)->type)
						(yyvsp[-1].expr)->type = (yyvsp[0].expr)->type;
					else if ((yyvsp[0].expr)->type != (yyvsp[-1].expr)->type)
					{
						(yyvsp[0].expr)->type = (yyvsp[-1].expr)->type;
						(yyvsp[0].expr)->data.operand.right = excast(expr.program, (yyvsp[0].expr)->data.operand.right, (yyvsp[-1].expr)->type, NULL, 0);
					}
					(yyvsp[0].expr)->data.operand.left = (yyvsp[-1].expr);
					(yyval.expr) = (yyvsp[0].expr);
				}
			}
		}
#line 2949 "exparse.c"
    break;

  case 95:
#line 952 "../../lib/expr/exparse.y"
                {
		pre:
			if ((yyvsp[0].expr)->type == STRING)
				exerror("++ and -- invalid for string variables");
			(yyval.expr) = exnewnode(expr.program, (yyvsp[-1].op), 0, (yyvsp[0].expr)->type, (yyvsp[0].expr), NULL);
			(yyval.expr)->subop = PRE;
		}
#line 2961 "exparse.c"
    break;

  case 96:
#line 960 "../../lib/expr/exparse.y"
                {
		pos:
			if ((yyvsp[-1].expr)->type == STRING)
				exerror("++ and -- invalid for string variables");
			(yyval.expr) = exnewnode(expr.program, (yyvsp[0].op), 0, (yyvsp[-1].expr)->type, (yyvsp[-1].expr), NULL);
			(yyval.expr)->subop = POS;
		}
#line 2973 "exparse.c"
    break;

  case 97:
#line 968 "../../lib/expr/exparse.y"
                {
			if ((yyvsp[0].id)->local.pointer == 0)
              			exerror("cannot apply IN to non-array %s", (yyvsp[0].id)->name);
			if (((yyvsp[0].id)->index_type > 0) && ((yyvsp[-2].expr)->type != (yyvsp[0].id)->index_type))
            		    exerror("%s indices must have type %s, not %s", 
				(yyvsp[0].id)->name, extypename(expr.program, (yyvsp[0].id)->index_type),extypename(expr.program, (yyvsp[-2].expr)->type));
			(yyval.expr) = exnewnode(expr.program, IN_OP, 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.variable.symbol = (yyvsp[0].id);
			(yyval.expr)->data.variable.index = (yyvsp[-2].expr);
		}
#line 2988 "exparse.c"
    break;

  case 98:
#line 979 "../../lib/expr/exparse.y"
                {
			goto pre;
		}
#line 2996 "exparse.c"
    break;

  case 99:
#line 983 "../../lib/expr/exparse.y"
                {
			goto pos;
		}
#line 3004 "exparse.c"
    break;

  case 103:
#line 993 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, CONSTANT, 0, (yyvsp[0].id)->type, NULL, NULL);
			if (!expr.program->disc->reff)
				exerror("%s: identifier references not supported", (yyvsp[0].id)->name);
			else
				(yyval.expr)->data.constant.value = expr.program->disc->reff(expr.program, (yyval.expr), (yyvsp[0].id), NULL);
		}
#line 3016 "exparse.c"
    break;

  case 104:
#line 1001 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, CONSTANT, 0, FLOATING, NULL, NULL);
			(yyval.expr)->data.constant.value.floating = (yyvsp[0].floating);
		}
#line 3025 "exparse.c"
    break;

  case 105:
#line 1006 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, CONSTANT, 0, INTEGER, NULL, NULL);
			(yyval.expr)->data.constant.value.integer = (yyvsp[0].integer);
		}
#line 3034 "exparse.c"
    break;

  case 106:
#line 1011 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, CONSTANT, 0, STRING, NULL, NULL);
			(yyval.expr)->data.constant.value.string = (yyvsp[0].string);
		}
#line 3043 "exparse.c"
    break;

  case 107:
#line 1016 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, CONSTANT, 0, UNSIGNED, NULL, NULL);
			(yyval.expr)->data.constant.value.integer = (yyvsp[0].integer);
		}
#line 3052 "exparse.c"
    break;

  case 113:
#line 1032 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = makeVar(expr.program, (yyvsp[-1].id), 0, 0, (yyvsp[0].reference));
		}
#line 3060 "exparse.c"
    break;

  case 114:
#line 1036 "../../lib/expr/exparse.y"
                {
			Exnode_t*   n;

			n = exnewnode(expr.program, DYNAMIC, 0, (yyvsp[-2].id)->type, NULL, NULL);
			n->data.variable.symbol = (yyvsp[-2].id);
			n->data.variable.reference = 0;
			if (((n->data.variable.index = (yyvsp[-1].expr)) == 0) != ((yyvsp[-2].id)->local.pointer == 0))
				exerror("%s: is%s an array", (yyvsp[-2].id)->name, (yyvsp[-2].id)->local.pointer ? "" : " not");
			if ((yyvsp[-2].id)->local.pointer && ((yyvsp[-2].id)->index_type > 0)) {
				if ((yyvsp[-1].expr)->type != (yyvsp[-2].id)->index_type)
					exerror("%s: indices must have type %s, not %s", 
						(yyvsp[-2].id)->name, extypename(expr.program, (yyvsp[-2].id)->index_type),extypename(expr.program, (yyvsp[-1].expr)->type));
			}
			if ((yyvsp[0].reference)) {
				n->data.variable.dyna =exnewnode(expr.program, 0, 0, 0, NULL, NULL);
				(yyval.expr) = makeVar(expr.program, (yyvsp[-2].id), (yyvsp[-1].expr), n, (yyvsp[0].reference));
			}
			else (yyval.expr) = n;
		}
#line 3084 "exparse.c"
    break;

  case 115:
#line 1056 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ID, 0, STRING, NULL, NULL);
			(yyval.expr)->data.variable.symbol = (yyvsp[0].id);
			(yyval.expr)->data.variable.reference = 0;
			(yyval.expr)->data.variable.index = 0;
			(yyval.expr)->data.variable.dyna = 0;
			if (!(expr.program->disc->flags & EX_UNDECLARED))
				exerror("unknown identifier");
		}
#line 3098 "exparse.c"
    break;

  case 116:
#line 1068 "../../lib/expr/exparse.y"
                {
			(yyval.integer) = 0;
		}
#line 3106 "exparse.c"
    break;

  case 117:
#line 1072 "../../lib/expr/exparse.y"
                {
			(yyval.integer) = -1;
		}
#line 3114 "exparse.c"
    break;

  case 118:
#line 1076 "../../lib/expr/exparse.y"
                {
			/* If DECLARE is VOID, its type is 0, so this acts like
			 * the empty case.
			 */
			if (INTEGRAL((yyvsp[-1].id)->type))
				(yyval.integer) = INTEGER;
			else
				(yyval.integer) = (yyvsp[-1].id)->type;
				
		}
#line 3129 "exparse.c"
    break;

  case 119:
#line 1089 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 3137 "exparse.c"
    break;

  case 120:
#line 1093 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[-1].expr);
		}
#line 3145 "exparse.c"
    break;

  case 121:
#line 1099 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 3153 "exparse.c"
    break;

  case 122:
#line 1103 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = (yyvsp[0].expr)->data.operand.left;
			(yyvsp[0].expr)->data.operand.left = (yyvsp[0].expr)->data.operand.right = 0;
			exfreenode(expr.program, (yyvsp[0].expr));
		}
#line 3163 "exparse.c"
    break;

  case 123:
#line 1111 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ',', 1, 0, exnewnode(expr.program, ',', 1, (yyvsp[0].expr)->type, (yyvsp[0].expr), NULL), NULL);
			(yyval.expr)->data.operand.right = (yyval.expr)->data.operand.left;
		}
#line 3172 "exparse.c"
    break;

  case 124:
#line 1116 "../../lib/expr/exparse.y"
                {
			(yyvsp[-2].expr)->data.operand.right = (yyvsp[-2].expr)->data.operand.right->data.operand.right = exnewnode(expr.program, ',', 1, (yyvsp[-2].expr)->type, (yyvsp[0].expr), NULL);
		}
#line 3180 "exparse.c"
    break;

  case 125:
#line 1122 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 3188 "exparse.c"
    break;

  case 126:
#line 1126 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
			if ((yyvsp[0].id)->type)
				exerror("(void) expected");
		}
#line 3198 "exparse.c"
    break;

  case 128:
#line 1135 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ',', 1, (yyvsp[0].expr)->type, (yyvsp[0].expr), NULL);
		}
#line 3206 "exparse.c"
    break;

  case 129:
#line 1139 "../../lib/expr/exparse.y"
                {
			Exnode_t*	x;
			Exnode_t*	y;

			(yyval.expr) = (yyvsp[-2].expr);
			for (x = (yyvsp[-2].expr); (y = x->data.operand.right); x = y);
			x->data.operand.right = exnewnode(expr.program, ',', 1, (yyvsp[0].expr)->type, (yyvsp[0].expr), NULL);
		}
#line 3219 "exparse.c"
    break;

  case 130:
#line 1149 "../../lib/expr/exparse.y"
                                {expr.declare=(yyvsp[0].id)->type;}
#line 3225 "exparse.c"
    break;

  case 131:
#line 1150 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, ID, 0, (yyvsp[-2].id)->type, NULL, NULL);
			(yyval.expr)->data.variable.symbol = (yyvsp[0].id);
			(yyvsp[0].id)->lex = DYNAMIC;
			(yyvsp[0].id)->type = (yyvsp[-2].id)->type;
			(yyvsp[0].id)->value = exnewnode(expr.program, 0, 0, 0, NULL, NULL);
			expr.procedure->data.procedure.arity++;
			expr.declare = 0;
		}
#line 3239 "exparse.c"
    break;

  case 132:
#line 1162 "../../lib/expr/exparse.y"
                {
			(yyval.reference) = expr.refs = expr.lastref = 0;
		}
#line 3247 "exparse.c"
    break;

  case 133:
#line 1166 "../../lib/expr/exparse.y"
                {
			Exref_t*	r;

			r = ALLOCATE(expr.program, Exref_t);
			memzero(r, sizeof(*r));
			r->symbol = (yyvsp[0].id);
			expr.refs = r;
			expr.lastref = r;
			r->next = 0;
			r->index = 0;
			(yyval.reference) = expr.refs;
		}
#line 3264 "exparse.c"
    break;

  case 134:
#line 1179 "../../lib/expr/exparse.y"
                {
			Exref_t*	r;
			Exref_t*	l;

			r = ALLOCATE(expr.program, Exref_t);
			memzero(r, sizeof(*r));
			r->symbol = (yyvsp[0].id);
			r->index = 0;
			r->next = 0;
			l = ALLOCATE(expr.program, Exref_t);
			memzero(l, sizeof(*l));
			l->symbol = (yyvsp[-1].id);
			l->index = 0;
			l->next = r;
			expr.refs = l;
			expr.lastref = r;
			(yyval.reference) = expr.refs;
		}
#line 3287 "exparse.c"
    break;

  case 135:
#line 1200 "../../lib/expr/exparse.y"
                {
			(yyval.id) = (yyvsp[0].id);
		}
#line 3295 "exparse.c"
    break;

  case 136:
#line 1204 "../../lib/expr/exparse.y"
                {
			(yyval.id) = (yyvsp[0].id);
		}
#line 3303 "exparse.c"
    break;

  case 137:
#line 1209 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = 0;
		}
#line 3311 "exparse.c"
    break;

  case 138:
#line 1213 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = exnewnode(expr.program, '=', 1, (yyvsp[0].expr)->type, NULL, (yyvsp[0].expr));
			(yyval.expr)->subop = (yyvsp[-1].op);
		}
#line 3320 "exparse.c"
    break;

  case 140:
#line 1220 "../../lib/expr/exparse.y"
                            {
				Dtdisc_t*	disc;

				if (expr.procedure)
					exerror("%s: nested function definitions not supported", expr.id->name);
				expr.procedure = exnewnode(expr.program, PROCEDURE, 1, expr.declare, NULL, NULL);
				if (!(disc = calloc(1, sizeof(Dtdisc_t))))
					exnospace();
				disc->key = offsetof(Exid_t, name);
				if (!streq(expr.id->name, "begin"))
				{
					if (!(expr.procedure->data.procedure.frame = dtopen(disc, Dtset)) || !dtview(expr.procedure->data.procedure.frame, expr.program->symbols))
						exnospace();
					expr.program->symbols = expr.program->frame = expr.procedure->data.procedure.frame;
					expr.program->formals = 1;
				}
				expr.declare = 0;
			}
#line 3343 "exparse.c"
    break;

  case 141:
#line 1237 "../../lib/expr/exparse.y"
                                  {
				expr.id->lex = PROCEDURE;
				expr.id->type = expr.procedure->type;
				expr.program->formals = 0;
				expr.declare = 0;
			}
#line 3354 "exparse.c"
    break;

  case 142:
#line 1243 "../../lib/expr/exparse.y"
                {
			(yyval.expr) = expr.procedure;
			expr.procedure = 0;
			if (expr.program->frame)
			{
				expr.program->symbols = expr.program->frame->view;
				dtview(expr.program->frame, NULL);
				expr.program->frame = 0;
			}
			(yyval.expr)->data.operand.left = (yyvsp[-5].expr);
			(yyval.expr)->data.operand.right = excast(expr.program, (yyvsp[-1].expr), (yyval.expr)->type, NULL, 0);

			/*
			 * NOTE: procedure definition was slipped into the
			 *	 declaration initializer statement production,
			 *	 therefore requiring the statement terminator
			 */

			exunlex(expr.program, ';');
		}
#line 3379 "exparse.c"
    break;


#line 3383 "exparse.c"

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
            yymsg = YY_CAST (char *, YYSTACK_ALLOC (YY_CAST (YYSIZE_T, yymsg_alloc)));
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
                  yystos[+*yyssp], yyvsp);
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
#line 1265 "../../lib/expr/exparse.y"


const char *exop(size_t index) {

  /* yytname is generated by the %token-table directive */

  /* find the index of MINTOKEN */
  size_t minid;
  for (minid = 0; yytname[minid] != NULL; ++minid) {
    if (strcmp(yytname[minid], "MINTOKEN") == 0) {
      break;
    }
  }

  assert(yytname[minid] != NULL
    && "failed to find MINTOKEN; incorrect token list in exparse.y?");

  /* find the requested token */
  {
    size_t i, j;
    for (i = j = minid; yytname[i] != NULL; ++i) {

      /* if this token is not a word, skip it */
      {
        size_t k;
        for (k = 0; yytname[i][k] != '\0'; ++k) {
          if (yytname[i][k] != '_' && !isalnum((int)yytname[i][k])) {
            break;
          }
        }
        if (yytname[i][k] != '\0') {
          continue;
        }
      }

      if (j == index + minid) {
        return yytname[i];
      }
      ++j;
    }
  }

  /* failed to find the requested token */
  return NULL;
}

void ex_error(const char *message) {
  exerror("%s", message);
}

#include <expr/exgram.h>
