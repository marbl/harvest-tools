AC_INIT(src/harvest/harvest.cpp)

AC_ARG_WITH(protobuf, [  --with-protobuf=<path/to/protobuf>     Protocol Buffers install dir (default: /usr/local/)])
AC_ARG_WITH(capnp, [  --with-capnp=<path/to/capnp>     Cap'n Proto install dir (default: /usr/local/)])

if test "$with_protobuf" == ""
then
	with_protobuf=/usr/local/
fi

if test "$with_capnp" == ""
then
	with_capnp=/usr/local/
fi

AC_LANG(C++)

AC_CHECK_PROG(protocCheck, protoc, "yes", "no", $with_protobuf/bin)

if test "$protocCheck" != "yes"
then
	AC_MSG_ERROR([Protocol Buffer compiler (protoc) not found.])
fi

AC_CHECK_PROG(capnpCheck, capnp, "yes", "no", $with_capnp/bin)

if test "$capnpCheck" != "yes"
then
	AC_MSG_ERROR([Cap'n Proto compiler (capnp) not found.])
fi

CPPFLAGS="-I$with_protobuf/include -I$with_capnp/include -std=c++14"

AC_CHECK_HEADER(google/protobuf/stubs/common.h, [result=1], [result=0])

if test $result == 0
then
	AC_MSG_ERROR([Protocol Buffers headers not found.])
fi

AC_CHECK_HEADER(capnp/common.h, [result=1], [result=0])

if test $result == 0
then
	AC_MSG_ERROR([Cap'n Proto headers not found.])
fi

AC_CHECK_HEADER(zlib.h, [result=1], [result=0])

if test $result == 0
then
	AC_MSG_ERROR([Zlib not found.])
fi

AC_SUBST(protobuf, $with_protobuf)
AC_SUBST(capnp, $with_capnp)

AC_OUTPUT(Makefile)
