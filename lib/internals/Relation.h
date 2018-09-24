/**
 * \file Relation.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: MapService.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#ifndef RELATION_H
#define RELATION_H

#include <list>
#include <map>
#include <boost/make_shared.hpp>

#include "dataclasses/I3Constants.h"
#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"

#include "ToolZ/OMKeyHash.h"
#include "ToolZ/Hitclasses.h"
#include "ToolZ/IndexMatrix.h"

#include "IceHiveZ/__SERIALIZATION.h"
static const unsigned relation_version_ = 0;

//forward declarations for serialization
#if SERIALIZATION_ENABLED
class Relation;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar,
    const Relation * t,
    const unsigned int version);

  template<class Archive>
  void load_construct_data(
    Archive & ar,
    Relation * t,
    const unsigned int version);
}};
#endif //SERIALIZATION_ENABLED

///Facilitates the access and interpretation of a indexMatrix as the connection between DOMs
class Relation {
  SET_LOGGER("Relation");
#if SERIALIZATION_ENABLED
  friend class SERIALIZATION_NS::access;

  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar,
    const Relation * t,
    const unsigned int version);

  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar,
    Relation * t,
    const unsigned int version);
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
private:
  ///the hasher that is used to adress the relation Map
  const CompactOMKeyHashServiceConstPtr hasher_;
  ///a relation map that can be evaluated with domIndexes
  indexmatrix::AsymmetricIndexMatrix_Bool relationMap_;

private: //utility classes
  /** utility for predicate
   */
  struct PredicateTranslator {
    const CompactOMKeyHashServiceConstPtr hasher_;
    const boost::function<bool (const OMKey&, const OMKey&)> predicate_;
    
    PredicateTranslator(
      const CompactOMKeyHashServiceConstPtr& hasher,
      const boost::function<bool (const OMKey&, const OMKey&)>& predicate)
    : hasher_(hasher),
      predicate_(predicate)
    {};
    
    bool operator() (const unsigned& first, const unsigned& second) {
      return predicate_(hasher_->OMKeyFromHash(first), hasher_->OMKeyFromHash(second));
    };
  }; //redicateTranslator

public: //methods
  ///Create a blank relationMap adressed by this Hasher
  Relation(
    const CompactOMKeyHashServiceConstPtr& hasher,
    const bool setall = false);
  /// constructor with predicate
  Relation(
    const CompactOMKeyHashServiceConstPtr& hasher,
    const boost::function<bool (const OMKey&, const OMKey&)> predicate);
  /// constructor by copying a relationMap
  Relation(
    const CompactOMKeyHashServiceConstPtr& hasher,
    const indexmatrix::AsymmetricIndexMatrix_Bool& relationMap);
  
public: //methods
  ///get the hasher
  CompactOMKeyHashServiceConstPtr GetHasher() const;
  ///get the relationMap
  const indexmatrix::AsymmetricIndexMatrix_Bool& GetRelationMap() const;
  
  /// if any of both are set, set the new one
  void Join(const Relation& r);
  /// if both are set, set the new one
  void Intersect(const Relation& r);
  
  /// Get the relation between two OMKeyHashes
  bool AreRelated(
    const CompactHash a,
    const CompactHash b) const;
  /// Set the relation between these OMKeyHashes
  void SetRelated(
    const CompactHash a,
    const CompactHash b,
    const bool value);
  
  ///get the relation between these DOMS  
  bool AreRelated(
    const OMKey& a,
    const OMKey& b) const;
  ///get set the relation for thess DOMs
  void SetRelated(
    const OMKey& a,
    const OMKey& b,
    const bool value);
  ///set all OMKeys to be related
  void SetAllRelated();
  ///set all OMKeys to be related
  void SetNoneRelated();
  ///set the relation by predication;
  ///NOTE assums function like object supports signature 'bool operator()(OMKey, OMKey)'
  void PredicateRelated(const boost::function<bool (const OMKey&, const OMKey&)>& callobj);
};

typedef boost::shared_ptr<Relation> RelationPtr;
typedef boost::shared_ptr<const Relation> RelationConstPtr;


#if SERIALIZATION_ENABLED
  SERIALIZATION_CLASS_VERSION(Relation, relation_version_);
#endif //SERIALIZATION_ENABLED


//=======================================================================
//======================= IMPLEMENTATION ================================
//=======================================================================


//=================== CLASS Relation =========
#if SERIALIZATION_ENABLED
template<class Archive>
void Relation::serialize(Archive & ar, const unsigned int version)
{};

// //NOTE (de)serialization overrides
namespace SERIALIZATION_NS_BASE { namespace serialization {
template<class Archive>
inline void save_construct_data(
  Archive & ar,
  const Relation * t,
  const unsigned int version)
{
  ar << SERIALIZATION_NS::make_nvp("Hasher", t->hasher_);
  ar << SERIALIZATION_NS::make_nvp("relationmap",t->relationMap_);
};

template<class Archive>
inline void load_construct_data(
  Archive & ar,
  Relation * t,
  const unsigned int version)
{
  CompactOMKeyHashServiceConstPtr hasher;
  ar >> SERIALIZATION_NS::make_nvp("Hasher", hasher); //const_cast<CompactOMKeyHashServiceConstPtr>(hasher_)
  indexmatrix::AsymmetricIndexMatrix_Bool relationMap(hasher->HashSize());
  ar >> SERIALIZATION_NS::make_nvp("relationmap",relationMap);
  
  ::new(t) Relation(hasher, relationMap);
};
}} // namespace ...
#endif //SERIALIZATION_ENABLED

inline 
CompactOMKeyHashServiceConstPtr 
Relation::GetHasher() const
  {return hasher_;};

inline
const indexmatrix::AsymmetricIndexMatrix_Bool& 
Relation::GetRelationMap() const
  {return relationMap_;};

inline
void Relation::Join(const Relation& r) //FIXME check congruence hashers
  {relationMap_ |= r.relationMap_;};
  
inline
void Relation::Intersect(const Relation& r)  //FIXME check congruence hashers
  {relationMap_ &= r.relationMap_;};

inline
bool Relation::AreRelated
  (const CompactHash a,
  const CompactHash b) const 
{return relationMap_.Get(a,b);};

inline
void Relation::SetRelated(
  const CompactHash a,
  const CompactHash b,
  const bool value)
{relationMap_.Set(a,b,value);};

inline
void Relation::SetAllRelated()
{
  for (size_t i=0; i<hasher_->HashSize(); i++) {   
    for (size_t j=i; j<hasher_->HashSize(); j++) {
      SetRelated(i, j, true);
      SetRelated(j, i, true);
    }
  }
};

inline
void Relation::SetNoneRelated()
{
  for (size_t i=0; i<hasher_->HashSize(); i++) {   
    for (size_t j=i; j<hasher_->HashSize(); j++) {
      SetRelated(i, j, false);
      SetRelated(j, i, false);
    }
  }
};

inline
void Relation::PredicateRelated(const boost::function<bool (const OMKey&, const OMKey&)>& callobj) 
{
  for (size_t i=0; i<hasher_->HashSize(); i++) {   
    for (size_t j=i; j<hasher_->HashSize(); j++) {
      SetRelated(i, j, callobj(hasher_->OMKeyFromHash(i),hasher_->OMKeyFromHash(j)));
      SetRelated(j, i, callobj(hasher_->OMKeyFromHash(j),hasher_->OMKeyFromHash(i)));
    }
  }
};

#endif //RELATION_H
