//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
