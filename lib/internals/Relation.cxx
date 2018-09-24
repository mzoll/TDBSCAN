/**
 * \file Relation.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: DistanceMapService.cxx 129048 2015-02-13 14:42:33Z mzoll $
 * \version $Revision: 129048 $
 * \date $Date: 2015-02-13 15:42:33 +0100 (fre, 13 feb 2015) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/internals/Relation.h"

#include <boost/foreach.hpp>

using namespace indexmatrix;

//=================== CLASS Relation =========

Relation::Relation(
  const CompactOMKeyHashServiceConstPtr& hasher,
  const AsymmetricIndexMatrix_Bool& relationMap)
: hasher_(hasher),
  relationMap_(relationMap)
{};

Relation::Relation(
  const CompactOMKeyHashServiceConstPtr& hasher,
  const bool setall) 
: hasher_(hasher),
  relationMap_(hasher->HashSize(), setall)
{};

Relation::Relation(
  const CompactOMKeyHashServiceConstPtr& hasher,
  const boost::function<bool (const OMKey&, const OMKey&)> predicate)
: hasher_(hasher),
  relationMap_(hasher->HashSize(), Relation::PredicateTranslator(hasher, predicate))
{};

bool Relation::AreRelated(const OMKey& a, const OMKey& b) const
{
  if (!hasher_)  
    log_fatal("no hasher set to perform this operation");
  return relationMap_.Get(hasher_->HashFromOMKey(a), hasher_->HashFromOMKey(b));
};
  
void Relation::SetRelated(const OMKey& a, const OMKey& b, const bool value) 
{
  if (!hasher_)  
    log_fatal("no hasher set to perform this operation");
  relationMap_.Set(hasher_->HashFromOMKey(a), hasher_->HashFromOMKey(b), value);
};

#if SERIALIZATION_ENABLED
  I3_SERIALIZABLE(Relation);
#endif //SERIALIZATON_ENABLED  
