AC_INIT(harvest)

AC_ARG_WITH(protobuf, [  --with-protobuf=<path/to/protobuf>     Protocol Buffers install dir (default: /usr/local/)])
AC_ARG_WITH(zlib, [  --with-zlib=<path/to/zlib>     Zlib install dir (default: /usr/local/)])

if test "$with_protobuf" == ""
then
	with_protobuf=/usr/local/
fi

if test "$with_zlib" == ""
then
	with_zlib=/usr/local/
fi

AC_LANG(C++)

AC_CHECK_PROG(protocCheck, protoc, "yes", "no" $with_protobuf/bin)

if test "$protocCheck" != "yes"
then
	AC_MSG_ERROR([Protocol Buffer compiler (protoc) not found.])
fi

CPPFLAGS="-I$with_protobuf/include -I$with_zlib/include"
AC_CHECK_HEADER(google/protobuf/stubs/common.h)
AC_CHECK_HEADER(zlib.h, [result=1], [result=0])

if test $result == 0
then
	AC_MSG_ERROR([Zlib not found.])
fi

AC_SUBST(protobuf, $with_protobuf)
AC_SUBST(zlib, $with_zlib/lib/libz.a)
AC_OUTPUT(Makefile)
