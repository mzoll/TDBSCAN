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



class Connector {
protected: // params
  ///a unique name for this service
  const std::string name_;
protected: //constructor
  Connector(const std::string& name);
public: //methods
  [[nodiscard]] std::string GetName() const;

    /** Are Hits h1 and h2 connected by being spatially and causally connected to each other?
  * @param h1
  * @param h2
  * @return true if hits are connected
  */
  template <class Position_t>
  bool Eval (const AbsBlib<Position_t>& h1, const AbsBlib<Position_t>& h2) const;
};


/**
 * A service which tells you if hits are (causally and topologically) connected
 */
class ConnectorSingle : public Connector {
private:
  friend class ConnectorBlock;
#if SERIALIZATION_ENABLED
  friend class SERIALIZATION_NS::access;

  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar,
    const Connector * t,
    const unsigned int version);

  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar,
    Connector* t,
    const unsigned int version);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int file_version);
#endif //SERIALIZATION_ENABLED
private: // params
  ///a unique name for this service
  const std::string name_;

public:
  ///constructor
  ConnectorSingle(
    const std::string& name);

public: //methods
  std::string GetName() const;

  /** Are Hits h1 and h2 connected by being spatially and causally connected to each other?
  * @param h1
  * @param h2
  * @return true if hits are connected
  */
  template <class Hitclass>
  bool Connected (const Hitclass& h1, const Hitclass& h2) const;
};

typedef std::shared_ptr<Connector> ConnectorPtr;
typedef std::shared_ptr<const Connector> ConnectorConstPtr;

#if SERIALIZATION_ENABLED
  SERIALIZATION_CLASS_VERSION(Connector, connector_version_);
#endif //SERIALIZATION_ENABLED

//============ CLASS ConnectorBlock ===========

/**
 * A collection of Connectors, which can be evaluated en-block.
 *
 */
class ConnectorBlock : public Connector {
#if SERIALIZATION_ENABLED
  friend class SERIALIZATION_NS::access;

  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar,
    const ConnectorBlock * t,
    const unsigned int version);

  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar,
    ConnectorBlock * t,
    const unsigned int version);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
public:
  ///list of connectors
  typedef std::list<ConnectorPtr> ConnectorList;
private: //property
  ///list of all connectors
  ConnectorList connectorlist_;

public: //constructors
  /// constructor purely with a geometry and no Connectors
  ConnectorBlock(
    const std::string& name,
    const ConnectorList cons = ConnectorList()
  );

public: //methods
  /// Add a Connector and add the Relation map to the cumRel Map
  void AddConnector (
    const ConnectorPtr& connector);

  ///check if to Hits are connected by any of the Connectors
  template <class Hitclass>
  bool Connected(
    const Hitclass& h1,
    const Hitclass& h2) const;

  ///diagnose the connection for these hits, as by which connector they are voted as connected
  template <class Hitclass>
  void DiagnoseConnected(
    const Hitclass& h1,
    const Hitclass& h2) const;

  //=== getters ===
  ///retrieve a connector from the ConnectorList; 0 will pass the cumulative one
  ConnectorPtr GetConnector (const int index) const;
  ///Get the complete list of Relations
  ConnectorList GetConnectorList() const;
};

typedef std::shared_ptr<ConnectorBlock> ConnectorBlockPtr;
typedef std::shared_ptr<const ConnectorBlock> ConnectorBlockConstPtr;

}


#endif //TDBSCAN_CONNECTOR_H
