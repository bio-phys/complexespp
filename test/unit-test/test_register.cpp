#include "complexes_test.h"
#include "gtest/gtest.h"

#include "io/rebuilder.h"
#include "io/serializer.h"
#include "util/util.h"

class ParentRebuild : public io::AbstractSerializable {
 public:
  ParentRebuild(io::Deserializer& deserializer) { UNUSED(deserializer); }
  virtual ~ParentRebuild() {}
  virtual std::string type() const = 0;

 protected:
  void serializeCore(io::Serializer& serializer) const {
    serializer.append(type(), "type");
  }
};

template <class T>
class ChildRebuild : public ParentRebuild {
 public:
  ChildRebuild(io::Deserializer& deserializer) : ParentRebuild(deserializer) {}
  static std::string Type() { return typeid(ChildRebuild<T>).name(); }
  std::string type() const final { return Type(); }
  void serialize(io::Serializer& serializer) const final {
    serializeCore(serializer);
  }
};

class RebuilderTest : public testing::Test {
 public:
  RebuilderTest() : deserializer(nullptr, 0) {}
  virtual ~RebuilderTest() {}

  io::Deserializer deserializer;

 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(RebuilderTest, simpleTest) {
  using MyRebuilder = io::Rebuilder<ParentRebuild>;

  MyRebuilder::Controller().registerType<ChildRebuild<float>>("float");

  auto f1 = []() {
    MyRebuilder::Controller().registerType<ChildRebuild<float>>("float");
  };
  EXPECT_NO_THROW(f1());

  auto f2 = []() {
    MyRebuilder::Controller().registerType<ChildRebuild<double>>("float");
  };
  EXPECT_THROW(f2(), std::runtime_error);

  auto f3 = [](io::Deserializer& des) {
    MyRebuilder::Controller().rebuild("do-not-exist", des);
  };
  EXPECT_THROW(f3(deserializer), std::runtime_error);

  auto f4 = [](io::Deserializer& des) {
    MyRebuilder::Controller().rebuild("float", des);
  };
  EXPECT_NO_THROW(f4(deserializer));
}

TEST_F(RebuilderTest, rebuildTest) {
  using MyRebuilder = io::Rebuilder<ParentRebuild>;

  MyRebuilder::Controller().registerType<ChildRebuild<float>>("float");

  MyRebuilder::Controller().registerType<ChildRebuild<float>>("float");
  EXPECT_EQ(MyRebuilder::Controller().rebuild("float", deserializer)->type(),
            ChildRebuild<float>::Type());

  MyRebuilder::Controller().registerType<ChildRebuild<double>>("double");
  EXPECT_EQ(MyRebuilder::Controller().rebuild("float", deserializer)->type(),
            ChildRebuild<float>::Type());
  EXPECT_EQ(MyRebuilder::Controller().rebuild("double", deserializer)->type(),
            ChildRebuild<double>::Type());
}
