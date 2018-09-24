/**
 * \file Hive-lib.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: Hive.cxx 151252 2016-11-01 14:33:39Z mzoll $
 * \version $Revision: 151252 $
 * \date $Date: 2016-11-01 15:33:39 +0100 (Tue, 01 Nov 2016) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/internals/Hive.h"

#include <algorithm>

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace hive;

//================= class DOMHoneyComb ===================

StringRings::StringRings(const StringNbr center) :
  rings_()
{
  rings_.push_back(boost::assign::list_of(center));
}

StringRings::StringRings(
  const StringNbr center,
  std::vector<Ring>& rings) :
  rings_()
{  
  rings_.push_back(boost::assign::list_of(center));
  BOOST_FOREACH(const Ring& ring, rings)
    rings_.push_back(ring);
}

StringNbr StringRings::StringRings::GetCenter() const {
  return *rings_[0].begin();
}

Ring StringRings::GetRing(const unsigned ringnbr) const {
  if (rings_.size() <= ringnbr)
    return Ring();
  else
    return rings_.at(ringnbr);
}

void StringRings::SetRing(const unsigned ringnbr, const Ring &ring) {
  if (rings_.size()<=ringnbr)
    rings_.resize(ringnbr+1);
  rings_[ringnbr]=ring;
}

void StringRings::AddStringToRing(
  const StringNbr string,
  const unsigned ringnbr)
{
  if (ringnbr<GetNRings()){
    rings_[ringnbr].insert(string);
  }
  else {
    Ring newring;
    newring.insert(string);
    SetRing(ringnbr, newring);
  }
};

bool StringRings::IsRingX (
  const StringNbr string,
  const unsigned ringnbr) const
{
  return GetRing(ringnbr).count(string);
}

int StringRings::WhichRing(
  const StringNbr string, 
  const unsigned max_search_depth) const 
{
  for (unsigned ringnbr =0; ringnbr<= GetNRings() && ringnbr<=max_search_depth; ringnbr++) {
    if (GetRing(ringnbr).count(string))
      return ringnbr;
  }
  log_trace_stream("Could not locate string rings around "<< *(GetRing(0).begin()));
  return -1;
};


std::string StringRings::Dump() const {
  std::ostringstream dump;
  dump << "========= DUMPING STRING "<<*rings_[0].begin()<<"==========" << std::endl;
  for (unsigned ring = 0; ring<=GetNRings(); ring++){
    dump << std::endl << " Ring "<<ring<<" : ";
    Ring strings = GetRing(ring);
    BOOST_FOREACH(const StringNbr& s, GetRing(ring)){
      dump << s << ", ";
    }
  }
  dump << std::endl;;
  return dump.str();
}

#if SERIALIZATION_ENABLED
  I3_SERIALIZABLE(StringRings);
#endif //SERIALIZATON_ENABLED  


//============ CLASS HiveTopology ==========

HiveTopology::HiveTopology() {};

void HiveTopology::AddStringRing(const StringRings& sr) {
  if (register_.count(sr.GetCenter()))
    log_warn("entry for string already exists: overwriting");
  register_.insert(std::make_pair(sr.GetCenter(),sr));
};
  
void HiveTopology::MutualAddStringToRing(
  const StringNbr center_A,
  const StringNbr center_B,
  const unsigned ringnbr) 
{
  StringRingRegister::iterator honey_A = register_.find(center_A);
  StringRingRegister::iterator honey_B = register_.find(center_B);
  
  if (honey_A==register_.end()){
    log_info_stream("String "<< center_A <<" not registered to the register_; Add it now");
    AddStringRing(StringRings(center_A));
    honey_A = register_.find(center_A);
  }
  
  if (honey_B==register_.end() ){
    log_info_stream("String "<< center_B <<" not registered to the register_; Add it now");
    AddStringRing(StringRings(center_B));
    honey_B = register_.find(center_B);
  }
  honey_A->second.AddStringToRing(center_B, ringnbr);
  honey_B->second.AddStringToRing(center_A, ringnbr);
};

bool HiveTopology::IsRingX (
  const StringNbr center,
  const StringNbr string,
  const unsigned ringnbr) const
{
  StringRingRegister::const_iterator honey = register_.find(center);
  if (honey==register_.end()) {
    log_trace_stream("Center string '"<< center<< "' not registered in map_"<< endl);
    return false;
  }
  return honey->second.IsRingX(string, ringnbr);
}

int HiveTopology::WhichRing(const StringNbr center,
              const StringNbr string, 
              const unsigned max_search_depth) const {
  StringRingRegister::const_iterator honey = register_.find(center);
  if (honey==register_.end()) {
    log_info_stream("Could not locate center string '"<< center << "' in map_" << std::endl);
    return -1;
  }
  return honey->second.WhichRing(string, max_search_depth);
}

/// Dump an entire hive in string format
std::string HiveTopology::Dump() const {
  std::ostringstream oss;
  Dump(oss);
  return oss.str();
};

std::ostream& HiveTopology::Dump(std::ostream& oss) const {
  BOOST_FOREACH(const StringRingRegister::value_type &entry, register_) {
    oss << "S "<< entry.first << std::endl;
    for (unsigned ring = 1; ring<= entry.second.GetNRings(); ring++) {
      oss << "R"<<ring<<" ";
      BOOST_FOREACH(const StringNbr& s, entry.second.GetRing(ring)) {
        oss << s << " ";
      }
      oss << std::endl;
    }
  }
  oss << std::endl;

  return oss;
};

HiveTopology
HiveTopology::FromDump(std::istream& iss) {
  HiveTopology ht;
  
  //containers
  std::string line, word;
  unsigned centerNbr, strNbr, ringNbr;
  
  //skip forward until there is no more line starting in a hash
  while(!iss.eof()){
    getline(iss, line);
    if (line[0]!='#')
      break;
  }
  
  //read the main data
  while(!iss.eof()){
    std::istringstream buffer(line);
    buffer>>word;
    if (word=="S") {
      bool more=true;
      //new string
      buffer>>centerNbr;
      
      ht.register_.insert(std::make_pair(centerNbr,StringRings(centerNbr))); //<unsigned,StringRings>
      
      more = true;
      while(more){
        getline(iss, line);
        if (line[0]=='R') {
          std::istringstream ring_buffer(line);
          ring_buffer>>word;
          word.erase(word.begin());
          std::stringstream(word)>>ringNbr;
          Ring str_in_ring;
          while(ring_buffer>>strNbr){
            str_in_ring.insert(strNbr);
          }
          ht.register_.at(centerNbr).SetRing(ringNbr, str_in_ring);
        }
        else{
          more=false;
          break;
        }
      }
    }
  }

  return ht;
};


std::ostream& operator<< (std::ostream& oss, const HiveTopology& ht) {
  return oss<<ht.Dump();
}

#if SERIALIZATION_ENABLED
  I3_SERIALIZABLE(HiveTopology);
#endif //SERIALIZATON_ENABLED  
