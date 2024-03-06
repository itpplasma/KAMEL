// Copyright (C) 2012 Ivan B. Ivanov. All rights reserved.

// Contact author: www.ivi.com or navi.adler@gmail.com.

// This file is part of ZeroSolver C++ library.

// ZeroSolver is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ZeroSolver is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ZeroSolver. If not, see <http://www.gnu.org/licenses/>.

/*! \file exceptions.hpp
    \brief A simple class hierarchy for basic exceptions
*/

#ifndef NAL_EXCEPTION_H
#define NAL_EXCEPTION_H

#include <stdexcept>
#include <iostream>
#include <sstream>

namespace exc
{
/********************************************************************/

class Exception : public std::exception
{
    public:

        Exception (std::string info, int Line = 0, const char * Func = 0, const char * File = 0) : exception()
        {
            line = 0;
            func = "";
            file = "";

            if (Line == 0) { message = info; return; }

            line = Line;
            func = std::string(Func);
            file = std::string(File);

            std::stringstream strstr;
            strstr << info + ": thrown at line " << line << " in function \'" << func << "\' in file \'" << file << "\'";
            message = strstr.str();
        }

        ~Exception (void) throw()
        {
        }

        virtual const char * what() const throw()
        {
            return message.c_str();
        }

    protected:

        std::string message; //!< exception message
        int         line;    //!< line number in that the exception throws
        std::string func;    //!< name of a function in that the exception throws
        std::string file;    //!< name of a file in that the exception throws
};

/********************************************************************/

class File_Open_Exception : public Exception
{
    public:

        File_Open_Exception (std::string info, int Line = 0, const char * Func = 0, const char * File = 0) : Exception (info, Line, Func, File)
        {
            message = "Failed to open file: " + message;
        }

        ~File_Open_Exception (void) throw()
        {
        }

        virtual const char * what() const throw()
        {
            return message.c_str();
        }
};

/********************************************************************/

class File_Read_Exception : public Exception
{
    public:

        File_Read_Exception (std::string info, int Line = 0, const char * Func = 0, const char * File = 0) : Exception (info, Line, Func, File)
        {
            message = "Failed to read from file (possibly bad format): " + message;
        }

        ~File_Read_Exception (void) throw()
        {
        }

        virtual const char * what() const throw()
        {
            return message.c_str();
        }
};

/********************************************************************/

class File_Write_Exception : public Exception
{
    public:

        File_Write_Exception (std::string info, int Line = 0, const char * Func = 0, const char * File = 0) : Exception (info, Line, Func, File)
        {
            message = "Failed to write to file: " + message;
        }

        ~File_Write_Exception (void) throw()
        {
        }

        virtual const char * what() const throw()
        {
            return message.c_str();
        }
};

/********************************************************************/

class Algorithm_Error_Exception : public Exception
{
    public:

        Algorithm_Error_Exception (std::string info, int Line = 0, const char * Func = 0, const char * File = 0) : Exception (info, Line, Func, File)
        {
            message = "Error: the internal check failed: " + message;
        }

        ~Algorithm_Error_Exception (void) throw()
        {
        }

        virtual const char * what() const throw()
        {
            return message.c_str();
        }
};

/********************************************************************/
} //end of namespace

#define throw_Exception(msg)                  throw exc::Exception((msg), __LINE__, __FUNCTION__, __FILE__)

#define throw_File_Open_Exception(msg)        throw exc::File_Open_Exception((msg), __LINE__, __FUNCTION__, __FILE__)

#define throw_File_Read_Exception(msg)        throw exc::File_Read_Exception((msg), __LINE__, __FUNCTION__, __FILE__)

#define throw_File_Write_Exception(msg)       throw exc::File_Write_Exception((msg), __LINE__, __FUNCTION__, __FILE__)

#define throw_Algorithm_Error_Exception(msg)  throw exc::Algorithm_Error_Exception((msg), __LINE__, __FUNCTION__, __FILE__)

#endif
