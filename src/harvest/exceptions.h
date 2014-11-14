
#include <stdexcept>

class BadInputFileException : public std::exception
{
public:
	
	BadInputFileException()
	{
	}
	
	virtual ~BadInputFileException() throw() {}
};

