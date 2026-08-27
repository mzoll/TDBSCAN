//
// Created by netsu on 08/08/2026.
//

#ifndef TDBSCAN_CONNECTOR_H
#define TDBSCAN_CONNECTOR_H


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




}


#endif //TDBSCAN_CONNECTOR_H
