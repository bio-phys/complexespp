#ifndef REBUILDER_HPP
#define REBUILDER_HPP

#include <functional>
#include <map>
#include <memory>
#include <string>

#include "io/serializer.h"
#include "util/util.h"

namespace io {

/*
  The REBUILDER class and MACRO can be used to generate a automated rebuild
  function for classes inheriting from a parent class.

```

// You can add more template arguments here if your deserializer constructor
// expects them
class Parent : public io::RebuilderCore<Parent>,
               public io::AbstractSerializable {
 public:
  // Here we declare a pure virtual serialize function to force the child-class
  // to define one too.
  virtual std::string type() const = 0;

 protected:
  void serializeCore(io::Serializer& serializer) const {
    // Store information of child type for rebuilding. This member has
    // to be called 'type' !!!
    // the rebuild function expects this
    serializer.append(type(), "type");
    // serializer any member functions here if defined.
  }
};

class Child : Parent {
 public:
  Child(io::Deserializer& deserializer)
      : Parent(deserializer), m_member(deserializer.restore<>()) {}
  void serialize(io::Serializer& serializer) {
    serializeCore(serializer);
    serializer.append(m_member, "m_member");
  }
  // define function to encode sub type information
  static const std::string& Type() { return __PRETTY_FUNCTION; }
  std::string type() { return Type(); }
};
REBUILDER_REGISTER(Child)
// The Rebuilder for child should be registered in a last step.

// To unpack the serialized child run
auto foo = Parent::rebuild(deserializer);
// The `rebuild` method is defined by the inheritance of RebuilderCore

```

*/

template <typename Target, typename... Params>
class Rebuilder {
  // The corresponding rebuild functions
  using FunctionType = std::function<std::unique_ptr<Target>(
      io::Deserializer& deserializer, Params...)>;
  // The dictionary of type/rebuild function
  std::map<std::string, FunctionType> rebuildFunctions;
  // To keep track of the original datatype
  std::map<std::string, std::string> rebuildTypenames;

  // Use singleton
  Rebuilder() = default;

 public:
  // Access singleton
  static Rebuilder& Controller() {
    static Rebuilder singleton;
    return singleton;
  }

  // Add a type/function pair
  template <typename InstanceClass>
  bool registerType(const std::string& type) {
    auto iter = rebuildTypenames.find(type);
    if (iter == rebuildTypenames.end()) {
      rebuildFunctions.emplace(type,
                               [](io::Deserializer& deserializer,
                                  Params... params) -> std::unique_ptr<Target> {
                                 return std::make_unique<InstanceClass>(
                                     deserializer, params...);
                               });
      rebuildTypenames.emplace(type, typeid(InstanceClass).name());
    } else if (iter->second != typeid(InstanceClass).name()) {
      throw std::runtime_error(
          fmt::format("Error in rebuilder, invalid register.\n"
                      "Type: {}\n"
                      "Already register name: {}\n"
                      "New name: {}\n",
                      type, iter->second, typeid(InstanceClass).name()));
    }
    return true;
  }

  // Rebuild a given type
  std::unique_ptr<Target> rebuild(const std::string& type,
                                  io::Deserializer& deserializer,
                                  Params... params) const {
    auto iter = rebuildFunctions.find(type);
    if (iter == rebuildFunctions.end()) {
      throw std::runtime_error(
          "Error in Rebuilder, cannot find the specified type (" + type + ")");
    }
    return (*iter).second(deserializer, params...);
  }
};

template <typename ParentClass, typename... Params>
class RebuilderCore {
 public:
  // Define the class in charge
  using RebuilderManager = Rebuilder<ParentClass, Params...>;
  // rebuild wrapper
  static std::unique_ptr<ParentClass> Rebuild(io::Deserializer& deserializer,
                                              Params... params) {
    const auto type = deserializer.restore<std::string>("type");
    return RebuilderManager::Controller().rebuild(type, deserializer,
                                                  params...);
  }
};
}  // namespace io

#define REBUILDER_REGISTER(CHILD)                                \
  namespace {                                                    \
  template <class C>                                             \
  class RebuildRegistration;                                     \
                                                                 \
  template <>                                                    \
  class RebuildRegistration<CHILD> {                             \
    static const bool reg;                                       \
  };                                                             \
                                                                 \
  const bool RebuildRegistration<CHILD>::reg =                   \
      CHILD::RebuilderManager::Controller().registerType<CHILD>( \
          CHILD::Type());                                        \
  }

#endif
