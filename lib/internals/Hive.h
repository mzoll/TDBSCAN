/**
 * \file Hive.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: Hive.h 153464 2017-02-22 14:28:08Z mzoll $
 * \version $Revision: 153464 $
 * \date $Date: 2017-02-22 15:28:08 +0100 (Wed, 22 Feb 2017) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#ifndef HIVE_H
#define HIVE_H

#include "IceHiveZ/__SERIALIZATION.h"

static const unsigned hive_version_ = 0;

#include <cstdlib>
#include <vector>
#include <set>
#include <map>

#include "icetray/OMKey.h"

//forward declaration for serialization
#if SERIALIZATION_ENABLED
namespace hive {
  class StringRings;
};

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar,
    const hive::StringRings * t,
    const unsigned int version);

  template<class Archive>
  void load_construct_data(
    Archive & ar,
    hive::StringRings * t,
    const unsigned int version);
}};
#endif //SERIALIZATION_ENABLED


/** 
 * In this namespace build objects have a ring-like structure, of a central string being 
 * surounded by rings upon rings of strings
 */
namespace hive {
  /// an alias for what it is
  typedef unsigned int StringNbr;
  /// a set of strings define a ring (around a central string)
  typedef std::set<StringNbr> Ring;
  
  /** A class that stores a central string and its surrounding rings;
   * NOTE Ring at position 0 in the rings_ is the center-string
   */
  class StringRings {
  #if SERIALIZATION_ENABLED
    friend class SERIALIZATION_NS::access;
    
    template<class Archive>
    friend void SERIALIZATION_NS::save_construct_data(
      Archive & ar,
      const StringRings * t,
      const unsigned int version);

    template<class Archive>
    friend void SERIALIZATION_NS::load_construct_data(
      Archive & ar,
      StringRings * t,
      const unsigned int version);
      
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version);
  #endif //SERIALIZATION_ENABLED
  private:
    /// member: storing the index of the ring [0...n] and the contained strings;
    std::vector<Ring> rings_;

  public:
    /// Constructor
    StringRings(const StringNbr center);
      
    /** @brief constructor with arbitrary number of rings
     * @param center the center
     * @param rings vector of rings and their contained strings
     */
    StringRings(const StringNbr center, std::vector<Ring>& rings);
  
  public: //methods
    /// get the central string
    StringNbr GetCenter() const;
    
    /// get the number of Rings
    unsigned int GetNRings() const;

    /// get the strings of the ring [ringnbr]
    /// \param ringnbr number of the ring that should be returned
    Ring GetRing(const unsigned ringnbr) const;
    
    /** Sets the strings for a specific ring, number, might erase a previous existing entry
     * @param ringnbr number of the ring to set
     * @param ring strings contained on that ring
     */
    void SetRing(const unsigned ringnbr, const Ring &ring);
    
    /** @brief Add the single string [string] to the ring at [ringnbr]
     * @param string that should be added
     * @param ringnbr the ring on which to add the string to
     */
    void AddStringToRing(
      const StringNbr string,
      const unsigned ringnbr);
      
    /** @brief is this [string] on ring[ringnbr] of [center]
     * @param string the string to be found
     * @param ringnbr ring number around center the string should be located on
     * @return true, if string could be correctly located;
     * false, if not;
     */
    bool IsRingX(
      const StringNbr string,
      const unsigned ringnbr) const;
                  
    /** @brief Which ring is this on?
      * @param string the string to locate on any ring
      * @param search_depth the maximum ring that this string is searched for (does limit the amount of computation)
      * @return 0, if it is the center;
      * n, if it is on ring n;
      * -1, if string is not registered;
      */
    int WhichRing(
      const StringNbr string,
      const unsigned max_search_depth=1000) const;
    
    /// Dump the entire entry
    std::string Dump() const;
  }; //class StringRings

  
  //=========== CLASS StringRingRegister ===============
  
  typedef std::map<StringNbr, StringRings> StringRingRegister;

  
  // ================ CLASS HiveTopology ========================
  
  ///complex type that holds all magic information associated
  class HiveTopology{
  #if SERIALIZATION_ENABLED
    friend class SERIALIZATION_NS::access;
    
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;
    
    template<class Archive>
    void load(Archive & ar, const unsigned int version);
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
  #endif //SERIALIZATION_ENABLED  
  private:
    /// the Register holding all information
    StringRingRegister register_;
    
  public: //constructors
    ///constructor blank
    HiveTopology();
    /// Construct a HiveTopology from the Dump in an inputstream
    static HiveTopology FromDump(std::istream& iss);
    
  public: //methods
    /// Add a StringRing to the register
    void AddStringRing(const StringRings& sr);
    
    /** @brief Register this string [center_A] as on ring[ringnbr] on string [center_B] and vise vers
     * @param center_A the string to register
     * @param center_B the other string to register
     * @param ringnbr the ring that strings should be registered to; should not be further than number of already registered rings plus 1
     */
    void MutualAddStringToRing(const StringNbr center_A,
                               const StringNbr center_B,
                               const unsigned ringnbr);
    
    /** @brief is this [string] on ring[ringnbr] of [center]
     * @param center center string around to look
     * @param string the string to be found
     * @param ringnbr ring number around center the string should be located on
     * @return true, if string could be correctly located;
     * false, if not;
     */
    bool IsRingX (const StringNbr center,
                  const StringNbr string,
                  const unsigned ringnbr) const;
                  
    /** @brief Which ring is this on?
      * @param center the center from which the ring should be found
      * @param string the string to locate on any ring
      * @param search_depth the maximum ring that this string is searched for (does limit the amount of computation)
      * @return 0, if it is the center;
      * n, if it is on ring n;
      * -1, if string is not registered;
      */
    int WhichRing(const StringNbr center,
                  const StringNbr string,
                  const unsigned max_search_depth=1000) const;
                                
    /// Does this object hold information of this string?
    /// \param string check this string
    bool HoldsCenterString(const StringNbr string) const;
    
    /// Dump an entire hive in string format
    std::string Dump() const;
    /// Dump an entire hive in to a stream object
    std::ostream& Dump(std::ostream& oss) const;
  };
  
  typedef boost::shared_ptr<HiveTopology> HiveTopologyPtr;
  typedef boost::shared_ptr<const HiveTopology> HiveTopologyConstPtr;
  
  std::ostream& operator<< (std::ostream& oss, const HiveTopology& ht);
};

#if SERIALIZATION_ENABLED
SERIALIZATION_CLASS_VERSION(hive::StringRings, hive_version_);
  SERIALIZATION_CLASS_VERSION(hive::HiveTopology, hive_version_);
#endif //SERIALIZATION_ENABLED



//==============================================================================
//========================= IMPLEMENTATIONS ====================================
//==============================================================================

#include <boost/foreach.hpp>

//========================= CLASS StringRings ========================

#if SERIALIZATION_ENABLED
template<class Archive>
void hive::StringRings::serialize(Archive & ar, const unsigned int version)
{};

// //NOTE (de)serialization overrides
template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar,
  const hive::StringRings * t,
  const unsigned int version)
{ ar << SERIALIZATION_NS::make_nvp("rings",t->rings_); };

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar,
  hive::StringRings * t,
  const unsigned int version)
{
  std::vector<hive::Ring> rings;
  ar >> SERIALIZATION_NS::make_nvp("rings",rings);
  const hive::StringNbr center = *(rings.at(0).begin());
  
  ::new(t) hive::StringRings(center, rings);
};
#endif //SERIALIZATION_ENABLED

/// get the number of implemented Rings
inline 
unsigned int hive::StringRings::GetNRings() const
  { return rings_.size()-1; };

//=========================== CLASS HiveTopology ====================
  
#if SERIALIZATION_ENABLED
template<class Archive>
void hive::HiveTopology::serialize(Archive & ar, const unsigned int version) 
// { ar & SERIALIZATION_NS::make_nvp("StringRingRegister", register_); };
//NOTE TODO does not work on serialization/trunk@2673); Problems in the deserialization of std::maps!
// replace by save/load/split_free
{ SERIALIZATION_NS::split_member(ar, *this, version); };

template<class Archive>
void hive::HiveTopology::save(Archive & ar, const unsigned int version) const
{ 
  const size_t size = register_.size();
  ar << SERIALIZATION_NS::make_nvp("size", size);
  
  BOOST_FOREACH(const StringRingRegister::value_type& sr ,register_)
    ar << SERIALIZATION_NS::make_nvp("StringRing", sr.second); 
};

template<class Archive>
void hive::HiveTopology::load(Archive & ar, const unsigned int version) 
{
  size_t n_many;
  ar >> SERIALIZATION_NS::make_nvp("size", n_many);
  for (size_t i=0; i<n_many; i++) {
    StringRings* sr = (StringRings*)malloc(sizeof(StringRings));
    ar >> SERIALIZATION_NS::make_nvp("StringRing", sr);
    AddStringRing(*sr);
  }
};
#endif //SERIALIZATION_ENABLED
  
inline
bool hive::HiveTopology::HoldsCenterString(const StringNbr string) const
  { return bool(register_.count(string)); };
  
#endif //HIVE_H
