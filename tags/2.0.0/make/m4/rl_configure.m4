# Ritchie Lab base autoconfigure script
# Checks for standard necessary libraries and adds options to 

# Arguments are as follows:
# 1 - Boost version requirement
# 2 - space separated list of HAVE_* required to run 
# 3 - Python version requirement
# 4 - comma or space separated list of python modules needed

AC_DEFUN([RL_CONFIGURE],
[
 
# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG

AC_PROG_CC([gcc cc icc])
AC_PROG_CPP
AC_PROG_CXX

# Checks for programs.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for --enable-loki
AX_CHECK_ENABLE_LOKI

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lmysqlclient':
AC_CHECK_LIB([mysqlclient], [main])
# FIXME: This needs to be fixed to correctly check for mysql client
AM_CONDITIONAL(HAVE_MYSQL,false)

AC_CHECK_LIB([pthread], [main])
AC_CHECK_LIB([z], [main])

# Checks for libraries.
AX_PATH_GSL([1.0],[AC_DEFINE([HAVE_GSL],[],[GSL Available])],[AC_DEFINE([NO_GSL],[],[No GSL Found])])

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
AX_BOOST_LIB([thread],[boost::thread_group thrds;],[])
AX_BOOST_LIB([regex],[boost::regex r();])
AX_BOOST_LIB([system],[boost::system::system_category],[],[boost/system/error_code.hpp])
AX_BOOST_LIB([filesystem],[boost::filesystem::path my_path("foo/bar/data.txt");],[$BOOST_SYSTEM_LIB])
AX_BOOST_LIB([program_options],[boost::program_options::options_description po;])

AC_LANG_PUSH(C)

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

# Check for the python version, but only if given!
if test "x$3" != "x"; then
	AC_PATH_PROG([PYTHON],[python])
	if test "x$PYTHON" = "x"; then
		AC_MSG_ERROR([Could not find python executable.  Please ensure that Python version $3 is in the PATH])
	else
		AX_PROG_PYTHON_VERSION([2.6],[rl_found_python="yes"],[rl_found_python="no"])
		if test "x$rl_found_python" != "xyes"; then
			AC_MSG_ERROR([The Python executable is too old.  Please update your python to version $3])
		fi
	fi
fi

rl_found_mods="yes"
for rl_mod in $4; do
	AX_PYTHON_MODULE([$rl_mod],[1])
done

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h fcntl.h fenv.h float.h inttypes.h limits.h malloc.h netdb.h netinet/in.h stddef.h stdint.h stdlib.h string.h strings.h sys/file.h sys/ioctl.h sys/mount.h sys/param.h sys/socket.h sys/timeb.h syslog.h unistd.h wchar.h])

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
