//
// Created by netsu on 08/08/2026.
//

#ifndef TDBSCAN_CONNECTOR_HH
#define TDBSCAN_CONNECTOR_HH


#include "absblib.h"

#include <list>
#include <string>
#include <memory>

//======================= Connector ============
namespace tdbscan {
  namespace detail {

    template <class tBlib>
    bool CausallyConnected(
      const tBlib& h1,
      const tBlib& h2,
      const Connector<tBlib>& connector)
    {
      if (h1.GetTime() > h2.GetTime())
        return CausallyConnected(h2, h1, connector); //recursive call to enforce time-order at this point
      return connector.Connected(h1, h2);
    }


  } // namespace tdbscan::detail


template <class tBlib>
void ConnectorBlock<tBlib>::addConnector(const ConnectorSingle<tBlib>* connector_ptr) {
  connectorlist_.push_back(connector_ptr);
};

template <class tBlib>
bool ConnectorBlock<tBlib>::eval(const tBlib& h1, const tBlib& h2) const {
  for (const auto& connector : connectorlist_) {
    if (connector->eval(h1, h2))
      return true;
  }
  return false;
}

template <class tBlib>
std::list<std::string> ConnectorBlock<tBlib>::diagnose
(const tBlib& h1, const tBlib& h2) const {
    std::list<std::string> result;
    for (const auto& connector : connectorlist_) {
      if (connector->eval(h1, h2))
        result.push_back(connector->getName());
    }
    return result;
  };



}


#endif //TDBSCAN_CONNECTOR_HH
