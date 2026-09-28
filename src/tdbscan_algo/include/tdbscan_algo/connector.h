//
// Created by mzoll on 08/08/2026.
//

#ifndef TDBSCAN_CONNECTOR_H
#define TDBSCAN_CONNECTOR_H


#include "absblib.h"

#include <list>
#include <string>
#include <memory>

//======================= Connector ============
namespace tdbscan {
template<class Blib_t>
class Connector {

  public:
  virtual ~Connector() = default;

  /** Are Hits h1 and h2 connected by being spatially and causally connected to each other?
    * @param h1
    * @param h2
    * @return true if hits are connected
    */
    [[nodiscard]]
    virtual bool eval(const Blib_t& h1, const Blib_t& h2) const = 0;
  };

  //============ CLASS ConnectorSingle ===========
  /**
   * A service which tells you if hits are (causally and topologically) connected
   */
  template<class Blib_t>
  class ConnectorSingle : public Connector<Blib_t> {
  protected: // params
    ///a unique name for this service
    const std::string name_; //className constructed
  protected: //constructor
    ConnectorSingle(const std::string& name) : name_(name) {};
  public: //methods
    [[nodiscard]]
    std::string getName() const {return name_;};
  };


  //============ CLASS ConnectorBlock ===========
  /**
   * A collection of Connectors, which can be evaluated en-block.
   */
  template <class tBlib>
class ConnectorAssembly : public Connector<tBlib> {
  public:
    typedef Connector<tBlib> Connector_t;
    typedef std::list<const Connector_t*> ConnectorList;
  protected: //property
    ///list of all connectors
    ConnectorList connectorlist_;

  public: //methods
    /// Add a Connector to the list of to be evaluated Connectors
    void addConnector (
      const ConnectorSingle<tBlib>* connector_ptr);

    /**
     * diagnose the connection for these hits, as by which connector they are voted as connected
     * @tparam HitClass
     * @param h1
     * @param h2
     * @return the names of the Connectors that see these Hits as connected
     */
    std::list<std::string> diagnose(
      const tBlib& h1,
      const tBlib& h2) const;

    //=== getters ===
    ///retrieve a connector from the ConnectorList; 0 will pass the cumulative one
    ConnectorSingle<tBlib>* getConnector (int index) const;
    ///Get the complete list of Relations
    ConnectorList getConnectorList() const;
  };


/**
 * chains Connectors together by a logical AND operation
 * @tparam tBlib
 */
  template<class tBlib>
  class ConnectorAssembly_AND : public ConnectorAssembly<tBlib> {
  public:
    bool eval(const tBlib& h1, const tBlib& h2) const;
  };

  /**
   * chains Connectors together by a logical OR operation
   * @tparam tBlib
   */
  template<class tBlib>
  class ConnectorAssembly_OR : public ConnectorAssembly<tBlib> {
  public:
    bool eval(const tBlib& h1, const tBlib& h2) const;
  };






}


#include "tdbscan_algo/connector.hh"

#endif //TDBSCAN_CONNECTOR_H
