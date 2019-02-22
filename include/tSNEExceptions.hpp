#ifndef TSNEEXCEPTIONS_HPP_0DB9DFD3
#define TSNEEXCEPTIONS_HPP_0DB9DFD3
#include <exception>

class BadValueException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "This value cannot be assigned.";
	}
};
extern BadValueException badValueEx;

class NotCalculatedYetException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Processing has not reached the required level to return a meaningful value here.";
	}
};
extern NotCalculatedYetException notCalculatedEx;

class SizeMismatchException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "The dimensions of the input do not match what we want do do with it.";
	}
};
extern SizeMismatchException sizeMismatchEx;

#endif /* end of include guard: TSNEEXCEPTIONS_HPP_0DB9DFD3 */
