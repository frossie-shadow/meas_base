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
#include "pybind11/stl.h"

#include "lsst/utils/python.h"
#include "lsst/meas/base/FlagHandler.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace base {

namespace {

void declareFlagDefinition(py::module &mod) {
    py::class_<FlagDefinition, std::shared_ptr<FlagDefinition>> cls(mod, "FlagDefinition");

    cls.def(py::init<>());
    cls.def(py::init<std::string, std::string, std::size_t>(), "name"_a, "doc"_a,
            "number"_a = FlagDefinition::number_undefined);

    cls.def_readwrite("name", &FlagDefinition::name);
    cls.def_readwrite("doc", &FlagDefinition::doc);
    cls.def_readwrite("number", &FlagDefinition::number);

    cls.def("__eq__", &FlagDefinition::operator==, py::is_operator());
    cls.def("__ne__", &FlagDefinition::operator!=, py::is_operator());
}

void declareFlagDefinitionList(py::module &mod) {
    py::class_<FlagDefinitionList, std::shared_ptr<FlagDefinitionList>> cls(mod, "FlagDefinitionList");

    cls.def(py::init<>());
    cls.def(py::init<std::initializer_list<FlagDefinition> const &>());

    cls.def("__getitem__",
            [](FlagDefinitionList const &self, int i) {
                auto cind = lsst::utils::python::cppIndex(self.size(), i);
                return self[cind];
            },
            py::is_operator());
    cls.def("__len__", &FlagDefinitionList::size);

    cls.def("getEmptyList", &FlagDefinitionList::getEmptyList);
    cls.def("getDefinition",
            (FlagDefinition(FlagDefinitionList::*)(std::size_t) const) & FlagDefinitionList::getDefinition,
            "index"_a);
    cls.def("getDefinition", (FlagDefinition(FlagDefinitionList::*)(std::string const &) const) &
                                     FlagDefinitionList::getDefinition,
            "name"_a);
    cls.def("hasDefinition", &FlagDefinitionList::hasDefinition, "name"_a);
    cls.def("addFailureFlag", &FlagDefinitionList::addFailureFlag, "doc"_a = "General Failure Flag");
    cls.def("add", &FlagDefinitionList::add, "name"_a, "doc"_a);
}

void declareFlagHandler(py::module &mod) {
    py::class_<FlagHandler, std::shared_ptr<FlagHandler>> cls(mod, "FlagHandler");

    cls.def(py::init<>());
    cls.def(py::init<afw::table::SubSchema const &, FlagDefinitionList const &, FlagDefinitionList const &>(),
            "s"_a, "flagDefs"_a, "exclDefs"_a = FlagDefinitionList::getEmptyList());

    cls.def_static("getFailureFlagName", &FlagHandler::getFailureFlagName);
    cls.def_static("addFields", &FlagHandler::addFields, "schema"_a, "prefix"_a, "flagDefs"_a,
                   "exclDefs"_a=FlagDefinitionList::getEmptyList());

    cls.def("getFlagNumber", &FlagHandler::getFlagNumber, "flagName"_a);
    cls.def("getFlagName", &FlagHandler::getFlagName, "i"_a);
    cls.def("getValue", (bool (FlagHandler::*)(afw::table::BaseRecord const &, std::size_t) const) &
                                FlagHandler::getValue,
            "record"_a, "i"_a);
    cls.def("getValue", (bool (FlagHandler::*)(afw::table::BaseRecord const &, std::string const &) const) &
                                FlagHandler::getValue,
            "record"_a, "flagName"_a);
    cls.def("setValue", (void (FlagHandler::*)(afw::table::BaseRecord &, std::size_t, bool) const) &
                                FlagHandler::setValue,
            "record"_a, "i"_a, "value"_a);
    cls.def("setValue", (void (FlagHandler::*)(afw::table::BaseRecord &, std::string const &, bool) const) &
                                FlagHandler::setValue,
            "record"_a, "flagName"_a, "value"_a);
    cls.def("getFailureFlagNumber", &FlagHandler::getFailureFlagNumber);
    cls.def("handleFailure", &FlagHandler::handleFailure, "record"_a, "error"_a = nullptr);
}

}  // <anonymous>

PYBIND11_PLUGIN(flagHandler) {
    py::module mod("flagHandler");

    declareFlagDefinition(mod);
    declareFlagDefinitionList(mod);
    declareFlagHandler(mod);

    return mod.ptr();
}

}  // base
}  // meas
}  // lsst
