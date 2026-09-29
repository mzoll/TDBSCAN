//
// Created by mzoll on 08/08/2026.
//

#ifndef TDBSCAN_CONNECTOR_HH
#define TDBSCAN_CONNECTOR_HH

#include "connector.h"

#include <list>
#include <string>

//======================= Connector ============
namespace tdbscan {

template <class tBlib>
void ConnectorAssembly<tBlib>::addConnector(const ConnectorSingle<tBlib>* connector_ptr) {
  connectorlist_.push_back(connector_ptr);
};


template <class tBlib>
std::list<std::string> ConnectorAssembly<tBlib>::diagnose
(const tBlib& h1, const tBlib& h2) const {
  std::list<std::string> result;
  for (const auto& connector : connectorlist_) {
    if (connector->eval(h1, h2))
      result.push_back(connector->getName());
  }
  return result;
};



template <class tBlib>
bool ConnectorAssembly_AND<tBlib>::eval(const tBlib& h1, const tBlib& h2) const {
  for (const auto& connector : ConnectorAssembly<tBlib>::connectorlist_) {
    if (! connector->eval(h1, h2))
      return false;
  }
  return true;
}


template <class tBlib>
bool ConnectorAssembly_OR<tBlib>::eval(const tBlib& h1, const tBlib& h2) const {
  for (const auto& connector : ConnectorAssembly<tBlib>::connectorlist_) {
    if (connector->eval(h1, h2))
      return true;
  }
  return false;
}






}  // namespace tdbscan;

#endif //TDBSCAN_CONNECTOR_HH
