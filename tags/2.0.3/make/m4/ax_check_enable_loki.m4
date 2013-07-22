AC_DEFUN([AX_CHECK_ENABLE_LOKI],[
    AC_MSG_CHECKING(whether to enable installation of the LOKI databse)
    m4_define(ax_enable_loki_default,[m4_tolower(m4_normalize(ifelse([$1],,[yes],[$1])))])
    AC_ARG_ENABLE(loki,
        [AS_HELP_STRING([--enable-loki]@<:@=ax_enable_loki_default@:>@,[Download and install the LOKI databse; one of yes/no])],
        [],enable_loki=ax_enable_loki_default)
    if test "x$enable_loki" = "xno"; then
        AC_MSG_RESULT(no)
    else
        AC_MSG_RESULT(yes)
        dnl assert.h is a NOP if NDEBUG is defined, so define it.
    fi
    ax_enable_loki=$enable_loki
    AM_CONDITIONAL(LOKI_INSTALL, test x$enable_loki = xyes)
])
