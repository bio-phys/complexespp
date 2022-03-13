// Copyright (c) 2018 the complexes++ development team and contributors
// (see the file AUTHORS for the full list of names)
//
// This file is part of complexes++.
//
// complexes++ is free software: you can redistribute it and/or modify
// it under the terms of the Lesser GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// complexes++ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with complexes++.  If not, see <https://www.gnu.org/licenses/>
#ifndef IO_YAML_H
#define IO_YAML_H

#include <algorithm>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <yaml-cpp/yaml.h>

#include "util/array.h"

namespace io {

class YamlNode {
public:
  explicit YamlNode(YAML::Node node) : m_node(node), m_level("/") {}
  explicit YamlNode(const std::string &file)
      : m_node(YAML::LoadFile(file)), m_level("/") {}

  template <typename T> T as() const {
    try {
      return m_node.as<T>();
    } catch (YAML::Exception &e) {
      throw std::runtime_error(fmt::format("ERROR in YAML: Reading key {}\n"
                                           "  LIBRARY ERROR: {}\n",
                                           m_level, e.what()));
    }
  }

  YamlNode operator[](std::string key) {
    return YamlNode(m_node[key], m_level + '.' + key);
  }

  const YamlNode operator[](std::string key) const {
    return YamlNode(m_node[key], m_level + '.' + key);
  }

  YamlNode operator[](int key) {
    return YamlNode(m_node[key], m_level + '.' + std::to_string(key));
  }

  const YamlNode operator[](int key) const {
    return YamlNode(m_node[key], m_level + '.' + std::to_string(key));
  }

  template <typename T> std::vector<T> keys() const {
    auto keys = std::vector<T>(m_node.size());
    std::transform(m_node.begin(), m_node.end(), keys.begin(),
                   [](const auto &el) { return el.first.template as<T>(); });
    return keys;
  }

  bool IsDefined() const { return m_node.IsDefined(); }
  std::size_t size() const { return m_node.size(); }
  auto Type() const -> decltype(YAML::Node().Type()) { return m_node.Type(); }
  const std::string &level() const { return m_level; }

private:
  explicit YamlNode(YAML::Node node, std::string level)
      : m_node(node), m_level(level) {
    if (!m_node.IsDefined()) {
      throw std::runtime_error(
          fmt::format("ERROR in YAML: Reading key {}\n"
                      "  LIBRARY ERROR: Node not defined\n",
                      m_level));
    }
  }

  YAML::Node m_node;
  std::string m_level;
};

template <typename T> std::vector<T> parseVector(const YamlNode &node) {
  auto vec = std::vector<T>();
  vec.reserve(node.size());
  for (auto i = 0u; i < node.size(); ++i) {
    vec.push_back(node[i].as<T>());
  }
  return vec;
}

template <typename T> util::vec<T> parse3DVector(const YamlNode &node) {
  if (node.size() != 3) {
    throw std::runtime_error(fmt::format(
        "YAML 3D Vec length not equal to 3 for key = {}\n", node.level()));
  }
  return util::vec<T>(node[0].as<T>(), node[1].as<T>(), node[2].as<T>());
}

template <typename T> util::Array<T> parseArray(const YamlNode &node) {
  const auto nrow = node.size();
  const auto ncol = node[0].size();

  auto arr = util::Array<T>(nrow, ncol);

  for (auto i = 0u; i < nrow; ++i) {
    if (node[i].size() != ncol) {
      throw std::runtime_error(
          fmt::format("In YAML key = {}, on row '{}' the number of columns "
                      "doesn't match. Got '{}' expected '{}'\n",
                      node.level(), i, node[i].size(), ncol));
    }
    for (auto j = 0u; j < ncol; ++j) {
      arr(i, j) = node[i][j].as<T>();
    }
  }
  return arr;
}

template <typename T, int rows, int cols>
util::FixedArray<T, rows, cols> parseFixedArray(const YamlNode &node) {
  if (node.size() != rows) {
    throw std::runtime_error(fmt::format(
        "In YAML key = {}, number of rows ([]) in node doesn't match expected "
        "number of {}.\n",
        node.level(), node.size(), rows));
  }
  if (node[0].size() != cols) {
    throw std::runtime_error(fmt::format(
        "In YAML key = {}, number of columns ([]) in node doesn't match "
        "expected number of {}.\n",
        node.level(), node[0].size(), cols));
  }

  auto arr = util::FixedArray<T, rows, cols>();

  for (auto i = 0u; i < rows; ++i) {
    if (node[i].size() != cols) {
      throw std::runtime_error(
          fmt::format("In YAML key = {}, on row '{}' the number of columns "
                      "doesn't match. Got '{}' expected '{}'\n",
                      node.level(), i, node[i].size(), cols));
    }
    for (auto j = 0u; j < cols; ++j) {
      arr(i, j) = node[i][j].as<T>();
    }
  }
  return arr;
}

} // namespace io

#endif // IO_YAML_H
