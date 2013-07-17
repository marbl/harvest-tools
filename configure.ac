AC_INIT(harvest)

AC_ARG_WITH(zlib, [  --with-zlib=<path/to/zlib>     Zlib install dir (default: /usr/local/])

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
AC_CHECK_HEADER(zlib.h)

AC_SUBST(protolib, $with_protobuf/lib/libprotobuf.a)
AC_SUBST(protoc, $with_protobuf/bin/protoc)
AC_SUBST(zlib, $with_zlib/lib/libz.a)
AC_OUTPUT(Makefile)
