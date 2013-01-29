# Ritchie Lab base autoconfigure script
# Checks for standard necessary libraries and adds options to 

# Arguments are as follows:
# 1 - Boost version requirement
# 2 - Comma separated list of HAVE_* required to run 

AC_DEFUN([RL_CONFIGURE],
[

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lmysqlclient':
AC_CHECK_LIB([mysqlclient], [main])
# FIXME: This needs to be fixed to correctly check for mysql client
AM_CONDITIONAL(HAVE_MYSQL,false)

AC_CHECK_LIB([pthread], [main])
AC_CHECK_LIB([z], [main])

# Checks for libraries.
# Biobin/ biofilter use sqlite 3.5.4 (group_concat)
PKG_CHECK_MODULES([SQLITE],[sqlite3 >= 3.5.4], [AC_DEFINE([HAVE_SQLITE],[],[SQLite is available])],[AC_DEFINE([NO_SQLITE],[],[No SQLite library available])])
PKG_CHECK_MODULES([PNG],[libpng >= 1.2.10], [AC_DEFINE([HAVE_PNG],[],[Png Library Available])],[AC_DEFINE([NO_PNG],[],[No PNG library available])])
PKG_CHECK_MODULES([R], [libR >= 2.10], [AC_DEFINE([HAVE_R],[],[R Library Available])],[AC_DEFINE([NO_R],[],[No R available])])
AX_FREETYPE_BASE([2.2.0])

AX_SOCI_CORE
AX_SOCI_SQLITE

# If desired, automake will have a conditional "USE_R" which will be true
# If not, "USE_R" will be false.  Also, will define the macro "USE_R" on the 
# command line for compilation flags.

#Check for boost
AX_BOOST_BASE([$1])
AX_BOOST_THREAD
AX_BOOST_REGEX
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS

#OK, now go through and make sure that we have everything that you requested:
rl_found_all="yes"
for rl_defn in $2; do
	AC_TRY_CPP([#ifndef $rl_defn
				#error "Library not found"
				#endif
			   ],[rl_found=yes],[rl_found=no])
	if test "$rl_found" != "yes"; then
		rl_lib=`echo $rl_defn | sed 's/^HAVE_//g'`
		AC_MSG_WARN([Could not find $rl_lib library])
		rl_found_all=no
	fi
done

	

if test "$rl_found_all" != "yes"; then
	AC_MSG_ERROR([Cound not find all prerequisite libraries.  Please check the library locations])
fi



# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h fcntl.h fenv.h float.h inttypes.h limits.h malloc.h netdb.h netinet/in.h stddef.h stdint.h stdlib.h string.h strings.h sys/file.h sys/ioctl.h sys/mount.h sys/param.h sys/socket.h sys/time.h sys/timeb.h syslog.h unistd.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM
AC_C_VOLATILE
AC_TYPE_UID_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_MKTIME
AC_FUNC_REALLOC
AC_FUNC_SELECT_ARGTYPES
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([bzero clock_gettime fdatasync floor ftime ftruncate getcwd gethrtime gettimeofday localtime_r memmove memset mkdir pow select socket sqrt strchr strrchr strerror])

])
