/*
 * LSST Data Management System
 * Copyright 2008-2017  AURA/LSST.
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

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <memory>

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/base/LocalBackground.h"
#include "lsst/meas/base/FluxUtilities.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

using PyAlgorithm = py::class_<LocalBackgroundAlgorithm, std::shared_ptr<LocalBackgroundAlgorithm>,
                               SimpleAlgorithm>;
using PyControl = py::class_<LocalBackgroundControl>;
using PyTransform = py::class_<LocalBackgroundTransform, std::shared_ptr<LocalBackgroundTransform>,
                               BaseTransform>;


PyControl declareControl(py::module &mod) {
    PyControl cls(mod, "LocalBackgroundControl");

    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, badMaskPlanes);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, annulusStart);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, annulusStop);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, bgRej);
    LSST_DECLARE_CONTROL_FIELD(cls, LocalBackgroundControl, bgIter);

    cls.def(py::init<>());

    return cls;
}


PyAlgorithm declareAlgorithm(py::module &mod) {
    PyAlgorithm cls(mod, "LocalBackgroundAlgorithm");

    cls.attr("FAILURE") = py::cast(LocalBackgroundAlgorithm::FAILURE);
    cls.attr("NO_GOOD_PIXELS") = py::cast(LocalBackgroundAlgorithm::NO_GOOD_PIXELS);

    cls.def(py::init<LocalBackgroundAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def(py::init<LocalBackgroundAlgorithm::Control const &, std::string const &, afw::table::Schema &,
            std::string const &>(),
            "ctrl"_a, "name"_a, "schema"_a, "logName"_a);
    return cls;
}

PyTransform declareTransform(py::module &mod) {
    PyTransform cls(mod, "LocalBackgroundTransform");

    cls.def(py::init<LocalBackgroundTransform::Control const &, std::string const &,
                     afw::table::SchemaMapper &>(),
            "ctrl"_a, "name"_a, "mapper"_a);

    return cls;
}

}  // <anonymous>

PYBIND11_PLUGIN(localBackground) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base.algorithm");
    py::module::import("lsst.meas.base.flagHandler");
    py::module::import("lsst.meas.base.fluxUtilities");
    py::module::import("lsst.meas.base.transform");

    py::module mod("localBackground");

    auto clsControl = declareControl(mod);
    auto clsAlgorithm = declareAlgorithm(mod);
    auto clsTransform = declareTransform(mod);

    clsAlgorithm.attr("Control") = clsControl;
    clsTransform.attr("Control") = clsControl;

    python::declareAlgorithm<LocalBackgroundAlgorithm, LocalBackgroundControl, LocalBackgroundTransform>(
            clsAlgorithm, clsControl, clsTransform);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
