// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <fmt/format.h>
#include <xdrfile/xdrfile_trr.h>
#include <xdrfile/xdrfile_xtc.h>

#include "io/xdr.h"
#include "util/file.h"

namespace io {

XDR::XDR(const std::string& file, const std::string& mode)
    : m_file(file), m_pXdrfile(xdrfile_open(m_file.c_str(), mode.c_str())) {
  if (m_pXdrfile == nullptr) {
    throw std::invalid_argument(
        fmt::format("Couldn't open {} for writing", m_file));
  }
}

XDR::~XDR() {
  const auto ok = xdrfile_close(m_pXdrfile);
  // NEVER EVER think about throwing here instead of logging.
  if (ok != 0) {
    fmt::print(stderr, "Couldn't close xdrfile, errorCode = {}", ok);
  }
}

XDRFILE* XDR::operator&() {
  if (m_pXdrfile == nullptr) {
    throw std::invalid_argument(fmt::format("damn pointer was lost"));
  }
  return m_pXdrfile;
}

void XDR::flush() {
  const auto ok = xdrfile_close(m_pXdrfile);
  if (ok != 0) {
    throw std::invalid_argument(fmt::format("couldn't flush xdr file."));
  }
  m_pXdrfile = xdrfile_open(m_file.c_str(), "a");
  if (m_pXdrfile == nullptr) {
    throw std::invalid_argument(
        fmt::format("Couldn't open {} for writing", m_file));
  }
}

std::array<float, DIM * DIM> topologyesBoxToXDRbox(const util::rvec& box) {
  return {{static_cast<float>(box[0] / 10), 0, 0, 0,
           static_cast<float>(box[1] / 10), 0, 0, 0,
           static_cast<float>(box[2] / 10)}};
}

util::Array<float, util::ColumnMajor> collectCoordinates(
    const int nBeads, const domains::Domains& model) {
  auto x = util::Array<float, util::ColumnMajor>(nBeads, DIM);
  auto idx = 0;
  for (const auto& dom : model) {
    for (auto i = 0; i < dom->nBeads(); ++i) {
      for (auto j = 0u; j < DIM; ++j) {
        // xtc/trr use nm as length unit. We use Angstrom.
        x(idx + i, j) = static_cast<float>(dom->xyz()(i, j) / 10);
      }
    }
    idx += dom->nBeads();
  }
  return x;
}

int beadsInModel(const domains::Domains& model) {
  auto nBeads = 0;
  for (const auto& dom : model) {
    nBeads += dom->nBeads();
  }
  return nBeads;
}

void writeXTC(XDR& xdr, const domains::Domains& model, const util::rvec& box,
              const int step, const double time) {
  const auto nBeads = beadsInModel(model);
  auto x = collectCoordinates(nBeads, model);
  auto xbox = topologyesBoxToXDRbox(box);
  // This give me 3 digits of accuracy after the point. Only then it will give
  // the same results as the PDB files. This is due to the very large boxes we
  // more often tend to use in complexes. So the 'normal' precision of 1000
  // isn't enough. With these settings we are good up to 1000 \AA boxes.
  const auto prec = 10000;
  // change time from ns to ps
  const auto result =
      write_xtc(&xdr, nBeads, step, static_cast<float>(time * 1e3),
                reinterpret_cast<float(*)[3]>(xbox.data()),
                reinterpret_cast<rvec*>(x.data()), prec);
  if (result != exdrOK) {
    throw std::invalid_argument(
        fmt::format("Error writting xtc coordinates, code = {}", result));
  }
}

void writeTRR(XDR& xdr, const domains::Domains& model, const util::rvec& box,
              const int step, const double time) {
  const auto nBeads = beadsInModel(model);
  auto x = collectCoordinates(nBeads, model);
  auto xbox = topologyesBoxToXDRbox(box);

  // defaults we don't have other appropriate values for
  const auto lambda = 0.0f;
  auto velocities = nullptr;
  auto forces = nullptr;

  // change time from ns to ps
  const auto result =
      write_trr(&xdr, nBeads, step, static_cast<float>(time * 1e3), lambda,
                reinterpret_cast<float(*)[3]>(xbox.data()),
                reinterpret_cast<rvec*>(x.data()), velocities, forces);
  if (result != exdrOK) {
    throw std::invalid_argument(
        fmt::format("Error writting trr coordinates, code = {}", result));
  }
}

XDR_TYPE getType(const std::string& file) {
  const auto suffix = util::fileSuffix(file);
  if (suffix == ".trr") {
    return XDR_TYPE::trr;
  } else if (suffix == ".xtc") {
    return XDR_TYPE::xtc;
  } else {
    throw std::invalid_argument(
        fmt::format("{} <-- unkown file format.\n", file));
  }
}

class XDRResult {
  util::rArray m_x;
  bool m_ok;
  int m_error;

 public:
  XDRResult(const int natoms, const util::Array<float, util::ColumnMajor>& x_,
            const int retcode)
      : m_x(natoms, 3), m_ok(!retcode), m_error(retcode) {
    for (auto i = 0; i < natoms; ++i) {
      for (auto j = 0; j < 3; ++j) {
        m_x(i, j) = x_(i, j) * 10;  // convert units;
      }
    }
  }
  bool ok() const { return m_ok; }
  const util::rArray& x() const { return m_x; }
  int error() const { return m_error; }
};

XDRResult readXDR(XDR& xd, const int natoms, const XDR_TYPE type) {
  std::array<float, DIM * DIM> xbox;
  util::Array<float, util::ColumnMajor> x(natoms, 3);
  int ok = -1;

  int step = 0;
  float time = 0;
  float lambda = 0;
  float prec = 0;
  switch (type) {
    case XDR_TYPE::trr:
      ok = read_trr(&xd, natoms, &step, &time, &lambda,
                    reinterpret_cast<float(*)[3]>(xbox.data()),
                    reinterpret_cast<rvec*>(x.data()), nullptr, nullptr);
      break;
    case XDR_TYPE::xtc:
      ok = read_xtc(&xd, natoms, &step, &time,
                    reinterpret_cast<float(*)[3]>(xbox.data()),
                    reinterpret_cast<rvec*>(x.data()), &prec);
      break;
  }
  return XDRResult(natoms, x, ok);
}

int readNbrAtoms(const std::string& filename, const XDR_TYPE type) {
  int natoms;
  switch (type) {
    case XDR_TYPE::trr:
      read_trr_natoms(const_cast<char*>(filename.c_str()), &natoms);
      break;
    case XDR_TYPE::xtc:
      read_xtc_natoms(const_cast<char*>(filename.c_str()), &natoms);
      break;
  }
  return natoms;
}

XDRReader::XDRReader(std::shared_ptr<domains::Domains> dom,
                     const util::rvec& box, const std::string& file)
    : BaseReader(dom, box),
      m_xdr(file, "r"),
      m_type(getType(file)),
      m_bReadLast(false),
      m_lastFrame(m_numAtomsModel, 3) {
  auto res = readXDR(m_xdr, m_numAtomsModel, m_type);
  m_bReadLast = res.ok();
  m_lastFrame = res.x();
  if (m_numAtomsModel != readNbrAtoms(file, m_type)) {
    throw std::runtime_error("xdr and model don't have same number  of atoms");
  }
}

XDRReader::~XDRReader() {}

void XDRReader::readXDRFrame() {
  // load last read frame into topology
  auto frame = m_lastFrame;
  auto total = 0;
  for (auto& dom : *m_dom) {
    auto xyz = dom->xyz();
    for (auto i = 0; i < dom->nBeads(); ++i) {
      xyz(i, 0) = frame(i + total, 0);
      xyz(i, 1) = frame(i + total, 1);
      xyz(i, 2) = frame(i + total, 2);
    }
    dom->setXyz(std::move(xyz));
    total += dom->nBeads();
  }
  // see if we can read another frame
  auto res = readXDR(m_xdr, m_numAtomsModel, m_type);
  m_bReadLast = res.ok();
  m_lastFrame = res.x();
}

bool XDRReader::hasNextFrame() const {
  return m_bReadLast;
}

std::shared_ptr<domains::Domains> XDRReader::nextFrame() {
  readXDRFrame();
  return m_dom;
}
}  // namespace io
