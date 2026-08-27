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
class ConnectorBlock : public Connector<tBlib> {
  public:
    typedef ConnectorSingle<tBlib> Connector_t;
    typedef std::list< Connector_t* > ConnectorList;
  private: //property
    ///list of all connectors
    ConnectorList connectorlist_;

  public: //constructors
    /// constructor purely with a geometry and no Connectors
    ConnectorBlock();

  public: //methods
    /// Add a Connector to the list of to be evaluated Connectors
    void addConnector (
      const ConnectorSingle<tBlib>* connector_ptr);

    ///check if to Hits are connected by any of the Connectors
    bool eval(const tBlib& h1, const tBlib& h2) const;

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
    ConnectorSingle<tBlib>* getConnector (const int index) const;
    ///Get the complete list of Relations
    ConnectorList getConnectorList() const;
  };
}


#include "tdbscan_algo/connector.hh"

#endif //TDBSCAN_CONNECTOR_H
