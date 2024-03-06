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

/*! \file
    \brief The declaration and definition of List class.
*/

#ifndef TYPE_LIST_HPP
#define TYPE_LIST_HPP

#include <iostream>
#include <cassert>

/********************************************************************/

namespace type
{
/********************************************************************/

template <typename T> class List_Node
{
    public:

        List_Node (const T & Data, List_Node<T> * Prev = 0, List_Node<T> * Next = 0) : data(Data), prev(Prev), next(Next)
        {
            if (prev) prev->next = this;
            if (next) next->prev = this;
        }

        ~List_Node (void)
        {
        }

        List_Node<T> * get_prev (void) const { return prev; }
        List_Node<T> * get_next (void) const { return next; }
            const T  & get_data (void) const { return data; }
                  T  & get_data (void)       { return data; }

        void set_data (const T & Data) { data = Data; }

        void write (std::ostream & output, int index) const
        {
            output << "\n" << index << "\t" << data;
        }

        void destruct_data (void) { delete data; data = 0; }

    private:

        T data;

        List_Node<T> * prev;
        List_Node<T> * next;
};

/********************************************************************/

template <typename T> class List
{
    public:

        List (void) : size(0), first(0), last(0), is_alloc(true)
        {
        }

        List (const List_Node<T> * first_, const List_Node<T> * last_, size_t size_ = 0)
        :
        size(size_), first(const_cast<List_Node<T> *>(first_)), last(const_cast<List_Node<T> *>(last_)), is_alloc(false)
        {
            if (size) return;

            List_Node<T> * node = first;

            while (node)
            {
                ++size;
                node = node->get_next();
            }
        }

        ~List (void);

        List_Node<T> * push_front (const T & Data);
        List_Node<T> * push_back  (const T & Data);

        List_Node<T> * insert_before (const T & Data, List_Node<T> * N);
        List_Node<T> * insert_after  (const T & Data, List_Node<T> * N);

        size_t get_size (void) const { return size; }

        List_Node<T> * get_first (void) const { return first; }
        List_Node<T> * get_last  (void) const { return last;  }

        friend std::ostream & operator << (std::ostream & output, const List<T> & A)
        {
            A.write(output);
            return output;
        }

        void write (std::ostream & output) const;

        void destruct_data (void);

    private:

        size_t size;

        List_Node<T> * first;
        List_Node<T> * last;

        bool is_alloc;
};

/********************************************************************/

template <typename T> List<T> :: ~List (void)
{
if (is_alloc)
{
    List_Node<T> * node = first;

    while (node)
    {
        List_Node<T> * next = node->get_next();
        delete node;
        node = next;
    }
}
}

/********************************************************************/

template <typename T> inline List_Node<T> * List<T> :: push_front (const T & Data)
{
List_Node<T> * address = new List_Node<T>(Data, 0, first);

if (!size) last = address;

++size;

return first = address;
}

/********************************************************************/

template <typename T> inline List_Node<T> * List<T> :: push_back (const T & Data)
{
List_Node<T> * address = new List_Node<T>(Data, last, 0);

if (!size) first = address;

++size;

return last = address;
}

/********************************************************************/

template <typename T> inline List_Node<T> * List<T> :: insert_before (const T & Data, List_Node<T> * N)
{
assert(N);

List_Node<T> * address = new List_Node<T>(Data, N->get_prev(), N);

++size;

if (N == first) first = address;

return address;
}

/********************************************************************/

template <typename T> inline List_Node<T> * List<T> :: insert_after (const T & Data, List_Node<T> * N)
{
assert(N);

List_Node<T> * address = new List_Node<T>(Data, N, N->get_next());

++size;

if (N == last) last = address;

return address;
}

/********************************************************************/

template <typename T> inline void List<T> :: write (std::ostream & output) const
{
List_Node<T> * node = first;

size_t index = 0;

while (node)
{
    node->write(output, index++);
    node = node->get_next();
}

output << std::endl;
}

/********************************************************************/

template <typename T> inline void List<T> :: destruct_data (void)
{
List_Node<T> * node = first;

while (node)
{
    node->destruct_data();
    node = node->get_next();
}
}

/********************************************************************/

} //end of namespace

/********************************************************************/

#endif
