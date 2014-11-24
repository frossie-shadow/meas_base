
// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

// Miscellaneous algorithms

%include "lsst/meas/base/ApertureFlux.i"

%feature("notabstract") lsst::meas::base::CircularApertureFluxAlgorithm;
%include "lsst/meas/base/CircularApertureFlux.h"

%feature("notabstract") lsst::meas::base::PsfFluxAlgorithm;
%include "lsst/meas/base/PsfFlux.h"

%feature("notabstract") lsst::meas::base::SincFluxAlgorithm;
%include "lsst/meas/base/SincFlux.h"

%feature("notabstract") lsst::meas::base::NaiveFluxAlgorithm;
%include "lsst/meas/base/NaiveFlux.h"

%feature("notabstract") lsst::meas::base::PeakLikelihoodFluxAlgorithm;
%include "lsst/meas/base/PeakLikelihoodFlux.h"

%feature("notabstract") lsst::meas::base::GaussianFluxAlgorithm;
%include "lsst/meas/base/GaussianFlux.h"

%feature("notabstract") lsst::meas::base::GaussianCentroidAlgorithm;
%include "lsst/meas/base/GaussianCentroid.h"

%feature("notabstract") lsst::meas::base::NaiveCentroidAlgorithm;
%include "lsst/meas/base/NaiveCentroid.h"

%feature("notabstract") lsst::meas::base::SdssCentroidAlgorithm;
%include "lsst/meas/base/SdssCentroid.h"

%feature("notabstract") lsst::meas::base::PixelFlagsAlgorithm;
%include "lsst/meas/base/PixelFlags.h"

%feature("notabstract") lsst::meas::base::SdssShapeAlgorithm;
%include "lsst/meas/base/SdssShape.h"
%template (apply) lsst::meas::base::SdssShapeAlgorithm::apply<float>;
%template (apply) lsst::meas::base::SdssShapeAlgorithm::apply<double>;