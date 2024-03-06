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
    \brief The declaration and definition of Tree class.
*/

#ifndef TYPE_TREE_HPP
#define TYPE_TREE_HPP

#include <iostream>
#include <cassert>

/********************************************************************/

namespace type
{
/********************************************************************/

template <typename T> class Tree_Node
{
    public:

        Tree_Node (void) : data(0), level(0), root(0), first_child(0), last_child(0), prev_sibling(0), next_sibling(0)
        {
        }

        Tree_Node (const T & Data,
                   int Level = 0,
                   Tree_Node<T> * Root = 0,
                   Tree_Node<T> * First_Child  = 0, Tree_Node<T> * Last_Child   = 0,
                   Tree_Node<T> * Prev_Sibling = 0, Tree_Node<T> * Next_Sibling = 0
                  )
                  :
                  data(Data),
                  level(Level),
                  root(Root),
                  first_child(First_Child),
                  last_child(Last_Child),
                  prev_sibling(Prev_Sibling),
                  next_sibling(Next_Sibling)
        {
        }

        ~Tree_Node (void);

        const T & get_data (void) const { return data; }
              T & get_data (void)       { return data; }

        int get_level (void) const { return level; }

        Tree_Node<T> *         get_root (void) const { return root;         }
        Tree_Node<T> *  get_first_child (void) const { return first_child;  }
        Tree_Node<T> *   get_last_child (void) const { return last_child;   }
        Tree_Node<T> * get_prev_sibling (void) const { return prev_sibling; }
        Tree_Node<T> * get_next_sibling (void) const { return next_sibling; }

        void set_data (const T & Data) { data = Data; }

        Tree_Node<T> * add_child (const T & Data);

        void write (std::ostream & output, int Iparent, int Ichild) const;

        void destruct_data (void);

    private:

        T data;

        int level;

        Tree_Node<T> * root;
        Tree_Node<T> * first_child;
        Tree_Node<T> * last_child;
        Tree_Node<T> * prev_sibling;
        Tree_Node<T> * next_sibling;
};

/********************************************************************/

template <typename T> Tree_Node<T> :: ~Tree_Node (void)
{
Tree_Node<T> * child = first_child;

while (child)
{
    Tree_Node<T> * next = child->next_sibling;
    delete child;
    child = next;
}
}

/********************************************************************/

template <typename T> inline Tree_Node<T> * Tree_Node<T> :: add_child (const T & Data)
{
Tree_Node<T> * address = new Tree_Node<T>(Data, level+1, this, 0, 0, last_child, 0);

if (last_child)
{
    last_child->next_sibling = address;
}
else
{
    first_child = address;
}

return last_child = address;
}

/********************************************************************/

template <typename T> inline void Tree_Node<T> :: write (std::ostream & output, int Iparent, int Ichild) const
{
output << "\n" << "(level, Iparent, Ichild) = (" << level << ", " << Iparent << ", " << Ichild << "): " << data;

int i = 0;

Tree_Node<T> * ptr = first_child;

while (ptr)
{
    ptr->write (output, Ichild, i++);
    ptr = ptr->next_sibling;
}
}

/********************************************************************/

template <typename T> inline void Tree_Node<T> :: destruct_data (void)
{
delete data; data = 0;

Tree_Node<T> * ptr = first_child;

while (ptr)
{
    ptr->destruct_data();
    ptr = ptr->next_sibling;
}
}

/********************************************************************/

template <typename T> class Tree
{
    public:

        Tree (void) : root(new Tree_Node<T>()), is_alloc(true)
        {
        }

        Tree (const T & Data) : root(new Tree_Node<T>(Data)), is_alloc(true)
        {
        }

        Tree (const Tree_Node<T> * R) : root(const_cast<Tree_Node<T> *>(R)), is_alloc(false)
        {
        }

        ~Tree (void)
        {
            if (is_alloc) delete root;
        }

        Tree_Node<T> * get_root (void) const { return root; }

        Tree_Node<T> * add_child (const T & Data) { return root->add_child(Data); }

        friend std::ostream & operator << (std::ostream & output, const Tree<T> & A)
        {
            A.write_simple(output);

            return output;
        }

        void write_simple (std::ostream & output) const
        {
            root->write(output, 0, 0);
            output << std::endl;
        }

        void destruct_data (void)
        {
            root->destruct_data();
        }

    private:

        Tree_Node<T> * root;

        bool is_alloc;
};

/********************************************************************/
} //end of namespace

/********************************************************************/

#endif
