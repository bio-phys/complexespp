// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef IO_TRAJECTORYFILE_H
#define IO_TRAJECTORYFILE_H

#include <fstream>

#include "io/serializer.h"
#include "io/xdr.h"

namespace io {

class Serialize;
class Deserializer;

std::string deduceType(const std::string& fname);

class TrajectoryFile : public AbstractSerializable {
 public:
  enum class Mode { app, in, out };

  explicit TrajectoryFile(const std::string& fname,
                          const Mode mode = Mode::app);

  XDR& xdr();
  std::ofstream& fstream();
  const std::string& fname() const noexcept;
  const std::string& type() const noexcept;
  void flush();

  // class movable
  TrajectoryFile(TrajectoryFile&& rhs) = default;
  TrajectoryFile& operator=(TrajectoryFile&& rhs) = default;
  // class not copyable
  TrajectoryFile(const TrajectoryFile& rhs) = delete;
  TrajectoryFile& operator=(const TrajectoryFile& rhs) = delete;

  void serialize(io::Serializer& serializer) const final;
  TrajectoryFile(io::Deserializer& deserializer);

 private:
  std::string m_fname;
  std::string m_type;
  Mode m_mode;
  std::unique_ptr<std::ofstream> m_fstream = nullptr;
  std::unique_ptr<XDR> m_xdr = nullptr;
};
}  // namespace io
#endif  // IO_TRAJECTORYFILE_H
