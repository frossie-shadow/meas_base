/* 
 * LSST Data Management System
 * Copyright 2008-2016  AURA/LSST.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include <memory>

#include "pybind11/pybind11.h"

#include "lsst/meas/base/CircularApertureFlux.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyApertureFluxClass = py::class_<CircularApertureFluxAlgorithm,
                                       std::shared_ptr<CircularApertureFluxAlgorithm>, ApertureFluxAlgorithm>;

}  // <anonymous>

PYBIND11_PLUGIN(circularApertureFlux) {
    py::module mod("circularApertureFlux");
    
    PyApertureFluxClass cls(mod, "CircularApertureFluxAlgorithm");

    cls.def(py::init<CircularApertureFluxAlgorithm::Control const &, std::string const &,
                     afw::table::Schema &, daf::base::PropertySet &>(),
            "ctrl"_a, "name"_a, "schema"_a, "metadata"_a);

    cls.def("measure", &CircularApertureFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
