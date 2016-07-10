//
// This file defines
//
//   template<typename T> T lexical_cast(const std::string &x);
//
// which converts a string to type T.  Currently the following types T are supported:
//
//   string  (trivial)
//   long
//   int
//   double
//   uint16_t
//
// Note that a more general implementation of lexical_cast is already defined in boost, and considered 
// for inclusion in C++ TR2, but it's not in C++11, so we're forced to define it here!
//

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// trivial case: convert string -> string
template<> string lexical_cast<string> (const string &x) { return x; }


inline bool is_all_spaces(const char *s)
{
    if (!s)
	throw runtime_error("fatal: NULL pointer passed to is_all_spaces()");

    for (;;) {
	if (!*s)
	    return true;
	if (!isspace(*s))
	    return false;
	s++;
    }
}


template<> long lexical_cast<long> (const string &x)
{ 
    const char *ptr = x.c_str();
    char *endptr = NULL;

    long ret = strtol(ptr, &endptr, 10);

    if (endptr == ptr)
	throw runtime_error(string("lexical_cast<long> failed: arg=\"") + x + "\"");
    if ((ret == LONG_MIN) || (ret == LONG_MAX))
	throw runtime_error(string("lexical_cast<long> failed: arg=\"") + x + "\"");
    if (!is_all_spaces(endptr))
	throw runtime_error(string("lexical_cast<long> failed: arg=\"") + x + "\"");
    
    return ret;
}


template<> int lexical_cast<int> (const string &x)
{
    long ret = lexical_cast<long> (x);

    if ((sizeof(int) != sizeof(long)) && ((ret < INT_MIN) || (ret > INT_MAX)))
	throw runtime_error(string("lexical_cast<int> failed: arg=\"") + x + "\"");

    return ret;
}


template<> uint16_t lexical_cast<uint16_t> (const string &x)
{
    long ret = lexical_cast<long> (x);
    
    if ((ret < 0) || (ret > 65535))
	throw runtime_error(string("lexical_cast<uint16_t> failed: arg=\"" ) + x + "\"");

    return ret;
}


template<> double lexical_cast<double> (const string &x)
{ 
    const char *ptr = x.c_str();
    char *endptr = NULL;

    double ret = strtod(ptr, &endptr);

    if (endptr == ptr)
	throw runtime_error(string("lexical_cast<double> failed: arg=\"") + x + "\"");
    if ((ret == LONG_MIN) || (ret == LONG_MAX))
	throw runtime_error(string("lexical_cast<double> failed: arg=\"") + x + "\"");
    if (!is_all_spaces(endptr))
	throw runtime_error(string("lexical_cast<double> failed: arg=\"") + x + "\"");
    
    return ret;
}


template<> float lexical_cast<float> (const string &x)
{
    // FIXME very minor loose end: check for overflow converting double->float
    return lexical_cast<double> (x);
}


// -------------------------------------------------------------------------------------------------
//
// Unit test


template<typename T> 
static void check_convert(const string &x, T y)
{
    T ret;

    try {
	ret = lexical_cast<T> (x);
    }
    catch (...) {
	cerr << "test_lexical_cast(): threw unexpected exception\n";
	exit(1);
    }

    if (fabs(double(ret) - double(y)) > 1.0e-5) {
	cerr << "test_lexical_cast(): didn't correctly convert\n";
	exit(1);
    }
}


template<typename T> 
static void check_convert_throws(const string &x)
{
    try {
	lexical_cast<T> (x);
    }
    catch (...) {
	return;
    }

    cerr << "test_lexical_cast(): didn't throw an exception as expected\n";
    exit(1);
}


void test_lexical_cast()
{
    check_convert<int>("0", 0);
    check_convert<int>("-0", 0);
    check_convert<int>("12", 12);
    check_convert<int>("-123", -123);
    check_convert<int>(" \t 1234  \n\t", 1234);

    check_convert_throws<int>("");
    check_convert_throws<int>("  ");
    check_convert_throws<int>("oops");
    check_convert_throws<int>(" oops ");
    check_convert_throws<int>("1234abc");
    check_convert_throws<int>("1234 abc");
    check_convert_throws<int>("0.1");

    check_convert<uint16_t>("0", 0);
    check_convert<uint16_t> ("65535", 65535);
    check_convert_throws<uint16_t> ("-1");
    check_convert_throws<uint16_t> ("65536");

    check_convert<double>("1.23", 1.23);
    check_convert<double>("-1.23e-5", -1.23e-5);
    check_convert<double>("-5", -5.0);
    check_convert<double>(".23", 0.23);
    check_convert<double>("-.034e3", -0.034e3);
    check_convert<double>("  0.03e20  ", 0.03e20);

    check_convert_throws<double>("");
    check_convert_throws<double>("  ");
    check_convert_throws<double>("oops");
    check_convert_throws<double>(" oops ");
    check_convert_throws<double>("5x");
    check_convert_throws<double>("-1.3e20x");

    cerr << "test_lexical_cast(): success\n";
}


}   // namespace ch_frb_io
