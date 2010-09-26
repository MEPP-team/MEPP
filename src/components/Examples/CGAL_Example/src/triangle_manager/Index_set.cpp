///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////

#include "Index_set.h"

void index_set::from_set_set_to_vector_set( const std::set< index_set >& set_triangle_constraints,
                                            std::vector< std::set< unsigned int > >& vector_triangle_constraints)
{
    std::set< index_set >::const_iterator its(set_triangle_constraints.begin()), ite(set_triangle_constraints.end());
    for(;its!=ite;++its)
    {
        vector_triangle_constraints.push_back(its->m_indices);
    }
};

void index_set::from_set_set_to_vector_set( const std::set< index_set >& set_triangle_constraints,
                                            std::vector< std::set< unsigned int >* >& vector_triangle_constraints)
{
    std::set< index_set >::const_iterator its(set_triangle_constraints.begin()), ite(set_triangle_constraints.end());
    for(;its!=ite;++its)
    {
        vector_triangle_constraints.push_back(new std::set< unsigned int >(its->m_indices));
    }
};

void index_set::from_set_set_to_vector_vector(  const std::set< index_set >& set_triangle_constraints,
                                                std::vector< std::vector< unsigned int > >& vector_triangle_constraints)
{
    std::set< index_set >::const_iterator its(set_triangle_constraints.begin()), ite(set_triangle_constraints.end());
    for(;its!=ite;++its)
    {
        std::vector< unsigned int > one_line_constraints;
        std::set< unsigned int >::const_iterator its2(its->m_indices.begin()), ite2(its->m_indices.end());
        for(;its2!=ite2;++its2)
        {
            one_line_constraints.push_back(*its2);
        }
        vector_triangle_constraints.push_back(one_line_constraints);
    }
};

void index_set::from_set_set_to_vector_vector(  const std::set< index_set >& set_triangle_constraints,
                                                std::vector< std::vector< int > >& vector_triangle_constraints)
{
    std::set< index_set >::const_iterator its(set_triangle_constraints.begin()), ite(set_triangle_constraints.end());
    for(;its!=ite;++its)
    {
        std::vector< int > one_line_constraints;
        std::set< unsigned int >::const_iterator its2(its->m_indices.begin()), ite2(its->m_indices.end());
        for(;its2!=ite2;++its2)
        {
            one_line_constraints.push_back(*its2);
        }
        vector_triangle_constraints.push_back(one_line_constraints);
    }
};

void index_set::from_set_set_to_set(const std::set< index_set >& set_triangle_constraints,
                                    std::set< unsigned int > & triangles)
{
    std::set< index_set >::const_iterator its(set_triangle_constraints.begin()), ite(set_triangle_constraints.end());
    for(;its!=ite;++its)
    {
        std::set< unsigned int >::const_iterator its2(its->m_indices.begin()), ite2(its->m_indices.end());
        for(;its2!=ite2;++its2)
        {
            triangles.insert(*its2);
        }
    }
}
//////////////////////////////////////////////////////////////////
void index_set::print() const
{
    std::set< unsigned int >::const_iterator itsind(m_indices.begin()),iteind(m_indices.end());
    for(;itsind!=iteind;++itsind)
    {
        std::cout << *itsind << " ";
    }
    std::cout << std::endl;
};
