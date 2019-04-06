//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2012 Jean-Marie Aubry
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef DEPENDENCYNODE_H_
#define DEPENDENCYNODE_H_

#include "../Definitions.h"
#include <list>

//#define VERBOSE_DEPENDENCY_NODE

/**
 * @brief Base class for the template family DependencyNode
 *
 * DependencyNodes are non-copyable by design, as the dependency relationship
 * is established in the constructor.
 */
class DependencyBase
{
public:
  DependencyBase(const DependencyBase&) = delete;
  DependencyBase& operator=(const DependencyBase&) = delete;
  
  DependencyBase() :
  m_dirty( true )
  {}
  
  virtual ~DependencyBase()
  {}
  
  bool isDirty() const
  {
    return m_dirty;
  }
  
  void setDirty()
  {
    if( !m_dirty )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
      std::cout << "Dirtying " << name() << ' ' << this << '\n';
#endif
      m_dirty = true;
      // Unlike Maya, we also dirty transitively
      setDependentsDirty();
    }
    // NB if this was dirty, we can assume that it's dependents were dirty also because this method
    // and the constructor are the only places where m_dirty can be set to true; in both cases the
    // dependents are also dirtied. So no need to call setDependentsDirty() if we are already dirty here.
  }
  
  void setDependentsDirty()
  {
    for( auto dep = m_dependents.begin(); dep != m_dependents.end(); ++dep )
    {
      ( *dep )->setDirty();
    }
  }
  
  void setClean()
  {
#ifdef VERBOSE_DEPENDENCY_NODE
    std::cout << "Setting clean " << name() << ' ' << this << '\n';
#endif
    m_dirty = false;
  }
  
  virtual const char* name() const = 0;
  
  void addDependent( DependencyBase* dependent )
  {
#ifdef VERBOSE_DEPENDENCY_NODE
    std::cout << name() << ' ' << this << " adding " << dependent->name() << ' ' << dependent
    << " as dependent" << '\n';
#endif
    m_dependents.push_back( dependent );
  }
  
protected:
  
  void setDirtyWithoutPropagating()
  {
    m_dirty = true;
  }
  
  virtual void compute() = 0;
  std::list<DependencyBase*> m_dependents;
  
private:
  bool m_dirty;
  
};

/**
 * @brief Template class to contain computable quantities and links to their dependencies.
 *
 * To register dependencies, a class derived from DependencyNode<ValueT> must have those dependencies
 * as reference member variables (thus needing their initialization in the constructor). Each dependency
 * must also call addDependent(this) in the constructor.
 */
template<typename ValueT>
class DependencyNode: public DependencyBase
{
public:
  DependencyNode( const ValueT& value ) :
  m_value( value )
  {}
  
  virtual ~DependencyNode()
  {}
  
  virtual const ValueT& get()
  {
    if( isDirty() )
    {
      compute();
#ifdef VERBOSE_DEPENDENCY_NODE
      std::cout << "Computed " << name() << ' ' << this << '\n';
#endif
      setClean();
    }
    
    return m_value;
  }
  
  virtual void set( const ValueT& value )
  {
    setDependentsDirty();
    m_value = value;
  }
  
  /**
   * @brief Erase m_value and free the memory
   */
  void free()
  {
    ValueT empty;
    std::swap( empty, m_value );
    // Do not propagate dirtiness, as the dependents are still valid
    // ( Free as not been called to signal that some input changed, but just to save memory )
    setDirtyWithoutPropagating();
  }
  
  protected    :
  ValueT m_value;
};

typedef Mat11Pair HessKType;
std::ostream& operator<<( std::ostream& os, const HessKType& HessKappa );

typedef std::pair<Mat2, Mat2> ThetaHessKType;
std::ostream& operator<<( std::ostream& os, const ThetaHessKType& HessKappa );

std::ostream& operator<<( std::ostream& os, const MatX& );
/**
 * @brief Specialized template for vector quantities
 */
template<typename ElemValueT, typename AllocatorT>
class DependencyNode<std::vector<ElemValueT, AllocatorT> > : public DependencyBase
{
public:
  typedef std::vector<ElemValueT, AllocatorT> ValueT;
  
  DependencyNode( IndexType firstValidIndex, IndexType size ) :
  m_firstValidIndex( firstValidIndex ), m_size( size )
  {}
  
  virtual ~DependencyNode()
  {}
  
  const ValueT& get()
  {
    if( isDirty() || m_value.size() != m_size )
    {
      compute();
#ifdef VERBOSE_DEPENDENCY_NODE
      std::cout << "Computed " << name() << '\n';
#endif
      setClean();
    }
    
    return m_value;
  }
  
  const ElemValueT& operator[]( IndexType i )
  {
    return get()[i];
  }
  
  void set( const ValueT& value )
  {
    setDependentsDirty();
    m_value = value;
  }
  
  virtual void cleanSet( const ValueT& value )
  {
    m_value = value;
  }
  
  const ValueT& getDirty()
  {
    return m_value;
  }
  
  
  virtual void set( IndexType i, const ElemValueT& elemVal )
  {
    setDependentsDirty();
    m_value.resize( m_size );
    m_value[i] = elemVal;
  }
  
  IndexType size() const
  {
    return m_size;
  }
  
  /**
   * @brief Erase m_value and free the memory
   */
  void free()
  {
    ValueT().swap( m_value );
    // Do not propagate dirtiness, as the dependents are still valid
    // ( Free as not been called to signal that some input changed, but just to save memory )
    setDirtyWithoutPropagating();
  }
  
  /**
   * @brief Erase m_value but capacity remains untouched
   */
  void clear()
  {
    m_value.clear();
    setDirty();
  }
  
  IndexType getFirstValidIndex()
  {
    return m_firstValidIndex;
  }
  
  virtual void print( std::ostream& os )
  {
    os << name() << ":...\n";
    for( IndexType i = getFirstValidIndex(); i < size(); ++i )
    {
      os << ( *this )[i] << ' ';
    }
    os << '\n' << name() << ":^^^\n";
  }
  
  protected    :
  // Either overload compute() or elemCompute(). The second one is the "lazy way", relying on this compute() loop,
  // but note that in involves calling get() for each iteration, hence unnecessary test. If you overload compute() instead,
  // get() once for each inputs and write you own loop.
  // IMPORTANT: due to the possibility to clear the m_value, always resize it before computing.
  
  virtual void compute()
  {
    m_value.resize( m_size );
    
    for( IndexType i = m_firstValidIndex; i < size(); ++i )
    {
      m_value[i] = elemCompute( i );
    }
    setDependentsDirty();
  }
  
  virtual ElemValueT elemCompute( IndexType )
  {
    std::cerr << "Either compute() or elemCompute() method must be implemented for class " << name() << std::endl;
    
    return ElemValueT();
  }
  
  IndexType m_firstValidIndex;
  ValueT m_value;
  size_t m_size;
};

template<typename ValueT>
inline std::ostream& operator<<( std::ostream& os, DependencyNode<ValueT>& node )
{
  node.print( os );
  
  return os;
}

#endif
