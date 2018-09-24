/**
 * \file ConfigInterfaces.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: MapService.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 * 
 * Facilitates Serialization of objects to Files, Write and Read operations
 */

#ifndef CONFIGINTERFACES_H
#define CONFIGINTERFACES_H

#include "icetray/I3Logging.h"

#include <fstream>

namespace configfileinterfaces {
//----------------------------
template<class T>
void 
WriteToConfigFile(
  const T& t,
  const std::string &filename);

template<class T>
void 
ReadFromConfigFile(
  T& t,
  const std::string &filename);

template<class T>
void
ReadFromConfigFile(
  T*& t,
  const std::string &filename);


}; //fileinterfaces

//============================================================
//=========================== IMPLEMETATIONS =================  
//============================================================

//============== from/to config files ================

template<class T>
void 
configfileinterfaces::WriteToConfigFile(
  const T& t,
  const std::string &filename)
{
  log_info_stream("Writing ConnectorBlock to config-file: "<<filename);

  std::ofstream ofs(filename.c_str());
  if (!ofs.is_open())
    log_fatal_stream("problems opening file " << filename.c_str());
  ofs << t.Dump();
  ofs.close();  
};

template<class T>
void
configfileinterfaces::ReadFromConfigFile(
  T& t,
  const std::string &filename)
{
  log_info_stream("Reading HiveTopology from config-file: "<<filename);
  
  std::ifstream ifs(filename.c_str());
  if (!ifs.is_open())
    log_fatal_stream(filename.c_str()<< " not found");
  
  t = T::FromDump(ifs);
  ifs.close();
};

//ptr types
template<class T>
void
configfileinterfaces::ReadFromConfigFile(
  T*& t,
  const std::string &filename)
{
  log_info_stream("Reading HiveTopology from config-file: "<<filename);
  
  std::ifstream ifs(filename.c_str());
  if (!ifs.is_open())
    log_fatal_stream(filename.c_str()<< " not found");
  
//   ::new(t)T::FromDump(ifs);
  t = new T(T::FromDump(ifs));
  ifs.close();
};



#endif //CONFIGINTERFACES_H
