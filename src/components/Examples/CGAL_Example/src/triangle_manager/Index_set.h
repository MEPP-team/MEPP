///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////

#ifndef INDEX_SET_H_INCLUDED
#define INDEX_SET_H_INCLUDED

#include <set>
#include <vector>
#include <iostream>

class index_set // this class is used to build a unique set of constraint!!
{
    public:
    std::set< unsigned int > m_indices; // one constraint

    public:
    // constructor
    index_set(const std::set< unsigned int >& indices):m_indices(indices) { }
    // destructor
    ~index_set(){ m_indices.clear(); }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline bool operator< (const index_set& id) const
    {
        // we firstly put the shortest constraints
        if(m_indices.size()<id.m_indices.size()) return true;
        else if(m_indices.size()>id.m_indices.size()) return false;

        if(m_indices==id.m_indices) return false;

        std::set< unsigned int >::const_iterator its(m_indices.begin()),ite(m_indices.end());
        std::set< unsigned int >::const_iterator itsC(id.m_indices.begin()),iteC(id.m_indices.end());
        for(;its!=ite && itsC!=iteC;++its, ++itsC)
        {
            if(*its <*itsC) return true;
            else if(*itsC<*its) return false;
        }

        // if the sets are equals at this point:
        return false;
    };
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static void from_set_set_to_vector_set( const std::set< index_set >& set_triangle_constraints,
                                            std::vector< std::set< unsigned int > >& vector_triangle_constraints);

    static void from_set_set_to_vector_set( const std::set< index_set >& set_triangle_constraints,
                                            std::vector< std::set< unsigned int >* >& vector_triangle_constraints);

    static void from_set_set_to_vector_vector(  const std::set< index_set >& set_triangle_constraints,
                                                std::vector< std::vector< unsigned int > >& vector_triangle_constraints);

    static void from_set_set_to_vector_vector(  const std::set< index_set >& set_triangle_constraints,
                                                std::vector< std::vector< int > >& vector_triangle_constraints);

    static void from_set_set_to_set(const std::set< index_set >& set_triangle_constraints,
                                    std::set< unsigned int > & triangles);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void print() const;
};

#endif // INDEX_SET_H_INCLUDED
