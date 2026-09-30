// SPDX-License-Identifier: MIT
#include "ignis/nozzle/NozzleGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "ignis/core/Exceptions.hpp"

namespace ignis {
namespace {
constexpr double kDeg = constants::pi / 180.0;
}

std::string toString(DivergentType t) {
  return t == DivergentType::kConical ? "conical" : "bell";
}

DivergentType divergentTypeFromString(const std::string& s) {
  if (s == "conical" || s == "cone") return DivergentType::kConical;
  if (s == "bell" || s == "parabolic" || s == "rao") return DivergentType::kBell;
  throw ConfigError("nozzle divergent type must be 'conical' or 'bell', got '" + s + "'");
}

std::string toString(NozzleSegment s) {
  switch (s) {
    case NozzleSegment::kChamber: return "chamber";
    case NozzleSegment::kChamberFillet: return "chamber_fillet";
    case NozzleSegment::kConvergingCone: return "converging_cone";
    case NozzleSegment::kThroatUpstream: return "throat_upstream";
    case NozzleSegment::kThroatDownstream: return "throat_downstream";
    case NozzleSegment::kDivergent: return "divergent";
  }
  return "?";
}

/// Analytic description of one contour piece.
struct NozzleGeometry::Segment {
  enum class Kind { kLine, kArc, kBezier } kind = Kind::kLine;
  NozzleSegment label = NozzleSegment::kChamber;
  double x0 = 0.0, x1 = 0.0;      // axial extent
  // line: r = r0 + slope * (x - x0)
  double r0 = 0.0, slope = 0.0;
  // arc: centre (xc, rc), radius R.  sign = +1 when the wall lies above the
  // centre (the convex chamber-to-cone fillet) and -1 when it lies below (the
  // concave throat arcs).
  double xc = 0.0, rc = 0.0, R = 0.0, sign = 1.0;
  // quadratic Bezier control points
  double bx0 = 0.0, br0 = 0.0, bxq = 0.0, brq = 0.0, bx1 = 0.0, br1 = 0.0;

  double radius(double x) const {
    switch (kind) {
      case Kind::kLine:
        return r0 + slope * (x - x0);
      case Kind::kArc: {
        const double dx = x - xc;
        const double disc = R * R - dx * dx;
        const double s = std::sqrt(std::max(disc, 0.0));
        return rc + sign * s;
      }
      case Kind::kBezier: {
        const double t = bezierParam(x);
        const double mt = 1.0 - t;
        return mt * mt * br0 + 2.0 * t * mt * brq + t * t * br1;
      }
    }
    return 0.0;
  }

  double slopeAt(double x) const {
    switch (kind) {
      case Kind::kLine:
        return slope;
      case Kind::kArc: {
        const double dx = x - xc;
        const double s = std::sqrt(std::max(R * R - dx * dx, 1e-300));
        return -sign * dx / s;
      }
      case Kind::kBezier: {
        const double t = bezierParam(x);
        const double dr = 2.0 * (1.0 - t) * (brq - br0) + 2.0 * t * (br1 - brq);
        const double dx = 2.0 * (1.0 - t) * (bxq - bx0) + 2.0 * t * (bx1 - bxq);
        return (dx != 0.0) ? dr / dx : 0.0;
      }
    }
    return 0.0;
  }

  /// Invert the monotone quadratic x(t) of the Bezier.
  double bezierParam(double x) const {
    const double a = bx0 - 2.0 * bxq + bx1;
    const double b = 2.0 * (bxq - bx0);
    const double c = bx0 - x;
    if (std::abs(a) < 1e-14 * std::max(1.0, std::abs(b))) {
      if (b == 0.0) return 0.0;
      return std::min(1.0, std::max(0.0, -c / b));
    }
    const double disc = b * b - 4.0 * a * c;
    const double sq = std::sqrt(std::max(disc, 0.0));
    const double t1 = (-b + sq) / (2.0 * a);
    const double t2 = (-b - sq) / (2.0 * a);
    const double t = (t1 >= -1e-9 && t1 <= 1.0 + 1e-9) ? t1 : t2;
    return std::min(1.0, std::max(0.0, t));
  }
};

NozzleGeometry::NozzleGeometry() = default;
NozzleGeometry::~NozzleGeometry() = default;
NozzleGeometry::NozzleGeometry(const NozzleGeometry&) = default;
NozzleGeometry& NozzleGeometry::operator=(const NozzleGeometry&) = default;
NozzleGeometry::NozzleGeometry(NozzleGeometry&&) noexcept = default;
NozzleGeometry& NozzleGeometry::operator=(NozzleGeometry&&) noexcept = default;

NozzleGeometry NozzleGeometry::build(const NozzleGeometrySpec& spec) {
  // A supplied contour wins over the analytic parameterisation, so that every
  // caller of build() -- the engine, the nozzle app, the transient -- picks up
  // imported geometry without knowing it exists.
  if (!spec.contour.empty()) return fromContour(spec.contour, spec);

  NozzleGeometry g;
  g.spec_ = spec;

  // --- throat and chamber sizes ----------------------------------------
  if (spec.throat_radius > 0.0 && spec.throat_area > 0.0)
    throw ConfigError("nozzle: give either throat_radius or throat_area, not both");
  if (spec.throat_radius > 0.0) {
    g.rt_ = spec.throat_radius;
  } else if (spec.throat_area > 0.0) {
    g.rt_ = std::sqrt(spec.throat_area / constants::pi);
  } else {
    throw ConfigError("nozzle: one of throat_radius or throat_area must be positive");
  }
  g.at_ = constants::pi * g.rt_ * g.rt_;

  if (spec.chamber_radius > 0.0 && spec.contraction_ratio > 0.0)
    throw ConfigError("nozzle: give either chamber_radius or contraction_ratio, not both");
  if (spec.chamber_radius > 0.0) {
    g.rc_ = spec.chamber_radius;
  } else if (spec.contraction_ratio > 0.0) {
    if (spec.contraction_ratio <= 1.0)
      throw ConfigError("nozzle: contraction_ratio must exceed 1");
    g.rc_ = g.rt_ * std::sqrt(spec.contraction_ratio);
  } else {
    throw ConfigError("nozzle: one of chamber_radius or contraction_ratio must be positive");
  }
  if (g.rc_ <= g.rt_)
    throw ConfigError("nozzle: chamber radius must exceed the throat radius");
  if (!(spec.expansion_ratio > 1.0))
    throw ConfigError("nozzle: expansion_ratio must exceed 1, got " +
                      std::to_string(spec.expansion_ratio));
  g.re_ = g.rt_ * std::sqrt(spec.expansion_ratio);
  if (!(spec.chamber_length >= 0.0))
    throw ConfigError("nozzle: chamber_length must be non-negative");

  const double beta = spec.converging_half_angle * kDeg;
  if (!(beta > 2.0 * kDeg && beta < 80.0 * kDeg))
    throw ConfigError("nozzle: converging_half_angle must lie in (2, 80) degrees");
  const double R1 = spec.chamber_fillet_ratio * g.rc_;
  const double Ru = spec.throat_upstream_ratio * g.rt_;
  const double Rd = spec.throat_downstream_ratio * g.rt_;
  if (!(R1 > 0.0) || !(Ru > 0.0) || !(Rd > 0.0))
    throw ConfigError("nozzle: fillet and throat curvature ratios must be positive");

  const double rA = g.rc_ - R1 * (1.0 - std::cos(beta));
  const double rB = g.rt_ + Ru * (1.0 - std::cos(beta));
  if (!(rA > rB)) {
    std::ostringstream os;
    os << "nozzle: the converging section does not fit -- the chamber fillet ends at r = "
       << rA << " m but the throat upstream arc starts at r = " << rB
       << " m. Reduce chamber_fillet_ratio (" << spec.chamber_fillet_ratio
       << ") or throat_upstream_ratio (" << spec.throat_upstream_ratio
       << "), or increase the contraction ratio.";
    throw ConfigError(os.str());
  }

  const double xA = spec.chamber_length + R1 * std::sin(beta);
  const double l_cone = (rA - rB) / std::tan(beta);
  const double xB = xA + l_cone;
  g.xt_ = xB + Ru * std::sin(beta);

  // --- divergent --------------------------------------------------------
  const double t15 = 15.0 * kDeg;
  g.l15_ = (g.rt_ * (std::sqrt(spec.expansion_ratio) - 1.0) +
            Rd * (1.0 / std::cos(t15) - 1.0)) / std::tan(t15);

  double theta_n, theta_e, x_end;
  if (spec.divergent == DivergentType::kConical) {
    theta_n = spec.cone_half_angle * kDeg;
    theta_e = theta_n;
    if (!(theta_n > 1.0 * kDeg && theta_n < 45.0 * kDeg))
      throw ConfigError("nozzle: cone_half_angle must lie in (1, 45) degrees");
  } else {
    if (!(spec.bell_length_fraction > 0.3 && spec.bell_length_fraction <= 1.2))
      throw ConfigError("nozzle: bell_length_fraction must lie in (0.3, 1.2]");
    // Default angle schedule when the user does not supply one.  These are the
    // commonly quoted values for an 80 %-length Rao bell; they are inputs, not
    // a reproduction of Rao's charts (see docs/limitations.md).
    theta_n = (spec.bell_initial_angle > 0.0 ? spec.bell_initial_angle : 33.0) * kDeg;
    theta_e = (spec.bell_exit_angle > 0.0 ? spec.bell_exit_angle : 8.0) * kDeg;
    if (!(theta_n > theta_e))
      throw ConfigError("nozzle: bell_initial_angle must exceed bell_exit_angle");
    if (!(theta_n < 60.0 * kDeg) || !(theta_e >= 0.0))
      throw ConfigError("nozzle: bell angles out of range (theta_n < 60 deg, theta_e >= 0)");
  }

  const double xN = g.xt_ + Rd * std::sin(theta_n);
  const double rN = g.rt_ + Rd * (1.0 - std::cos(theta_n));
  if (!(rN < g.re_))
    throw ConfigError("nozzle: the throat downstream arc already exceeds the exit radius; "
                      "reduce throat_downstream_ratio or the initial divergence angle");

  Segment div;
  div.label = NozzleSegment::kDivergent;
  if (spec.divergent == DivergentType::kConical) {
    div.kind = Segment::Kind::kLine;
    div.x0 = xN;
    div.r0 = rN;
    div.slope = std::tan(theta_n);
    x_end = xN + (g.re_ - rN) / div.slope;
    div.x1 = x_end;
  } else {
    x_end = g.xt_ + spec.bell_length_fraction * g.l15_;
    if (!(x_end > xN + 1e-9 * g.rt_))
      throw ConfigError("nozzle: bell_length_fraction is too small for the chosen "
                        "throat downstream arc and initial angle");
    const double tn = std::tan(theta_n), te = std::tan(theta_e);
    const double xQ = (g.re_ - rN + xN * tn - x_end * te) / (tn - te);
    const double rQ = rN + (xQ - xN) * tn;
    if (!(xQ > xN && xQ < x_end)) {
      std::ostringstream os;
      os << "nozzle: the bell control point falls outside the divergent (x_N = " << xN
         << " m, x_Q = " << xQ << " m, x_E = " << x_end
         << " m). Increase bell_length_fraction or reduce bell_initial_angle.";
      throw ConfigError(os.str());
    }
    if (!(rQ >= rN))
      throw ConfigError("nozzle: the bell control point lies below the throat arc exit");
    div.kind = Segment::Kind::kBezier;
    div.x0 = xN;
    div.x1 = x_end;
    div.bx0 = xN;  div.br0 = rN;
    div.bxq = xQ;  div.brq = rQ;
    div.bx1 = x_end; div.br1 = g.re_;
  }
  g.xe_ = x_end;
  g.theta_n_ = theta_n;
  g.theta_e_ = theta_e;

  // --- assemble segments ------------------------------------------------
  auto line = [](NozzleSegment lab, double x0, double x1, double r0, double slope) {
    Segment s;
    s.kind = Segment::Kind::kLine;
    s.label = lab; s.x0 = x0; s.x1 = x1; s.r0 = r0; s.slope = slope;
    return s;
  };
  auto arc = [](NozzleSegment lab, double x0, double x1, double xc, double rc, double R,
                double sign) {
    Segment s;
    s.kind = Segment::Kind::kArc;
    s.label = lab; s.x0 = x0; s.x1 = x1; s.xc = xc; s.rc = rc; s.R = R; s.sign = sign;
    return s;
  };

  if (spec.chamber_length > 0.0)
    g.segments_.push_back(line(NozzleSegment::kChamber, 0.0, spec.chamber_length, g.rc_, 0.0));
  // Chamber fillet: x = Lc + R1 sin(t), r = Rc - R1 (1 - cos t), so the centre
  // sits at (Lc, Rc - R1) and the wall lies *above* it.
  g.segments_.push_back(arc(NozzleSegment::kChamberFillet, spec.chamber_length, xA,
                            spec.chamber_length, g.rc_ - R1, R1, +1.0));
  if (l_cone > 0.0)
    g.segments_.push_back(line(NozzleSegment::kConvergingCone, xA, xB, rA, -std::tan(beta)));
  g.segments_.push_back(arc(NozzleSegment::kThroatUpstream, xB, g.xt_, g.xt_,
                            g.rt_ + Ru, Ru, -1.0));
  g.segments_.push_back(arc(NozzleSegment::kThroatDownstream, g.xt_, xN, g.xt_,
                            g.rt_ + Rd, Rd, -1.0));
  g.segments_.push_back(div);

  // --- Bartz throat curvature radius -----------------------------------
  if (spec.bartz_curvature == "upstream")      g.bartz_rc_ = Ru;
  else if (spec.bartz_curvature == "downstream") g.bartz_rc_ = Rd;
  else if (spec.bartz_curvature == "mean")     g.bartz_rc_ = 0.5 * (Ru + Rd);
  else throw ConfigError("nozzle: bartz_curvature must be 'upstream', 'downstream' or 'mean'");

  // --- sample -----------------------------------------------------------
  if (spec.num_stations < 21)
    throw ConfigError("nozzle: num_stations must be at least 21");
  double total = 0.0;
  for (const auto& s : g.segments_) total += (s.x1 - s.x0);
  std::vector<double> xs;
  for (const auto& s : g.segments_) {
    const double frac = (s.x1 - s.x0) / total;
    int n = std::max(6, static_cast<int>(std::lround(frac * spec.num_stations)));
    for (int i = 0; i <= n; ++i) xs.push_back(s.x0 + (s.x1 - s.x0) * i / n);
  }
  xs.push_back(g.xt_);
  xs.push_back(g.xe_);
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end(),
                       [&](double a, double b) { return std::abs(a - b) < 1e-12 * g.xe_; }),
           xs.end());

  g.stations_.reserve(xs.size());
  for (double x : xs) {
    const Segment& s = g.segmentFor(x);
    ContourPoint pnt;
    pnt.x = x;
    pnt.r = s.radius(x);
    pnt.area = constants::pi * pnt.r * pnt.r;
    pnt.area_ratio = pnt.area / g.at_;
    pnt.drdx = s.slopeAt(x);
    pnt.segment = s.label;
    g.stations_.push_back(pnt);
  }

  g.finalise();
  return g;
}

namespace {

/// Least-squares circle through (x, r) samples, Kasa's algebraic form.
///
/// Minimises sum of (x^2 + r^2 + D x + E r + F)^2, which is linear in
/// (D, E, F).  The centre is (-D/2, -E/2) and the radius follows.  It is
/// biased for short arcs, but the throat region is a well-sampled arc of
/// substantial angular extent and the alternative -- a proper geometric fit --
/// would need an iteration whose failure modes are harder to report.
/// Returns false when the normal equations are too ill-conditioned to trust,
/// which is what happens when the samples are collinear, i.e. a sharp throat.
bool fitCircle(const std::vector<ContourSample>& pts, double& radius_out) {
  const std::size_t n = pts.size();
  if (n < 4) return false;
  double sx = 0, sr = 0, sxx = 0, srr = 0, sxr = 0;
  double sxz = 0, srz = 0, sz = 0;
  for (const auto& p : pts) {
    const double z = p.x * p.x + p.r * p.r;
    sx += p.x; sr += p.r; sxx += p.x * p.x; srr += p.r * p.r; sxr += p.x * p.r;
    sxz += p.x * z; srz += p.r * z; sz += z;
  }
  const double dn = static_cast<double>(n);
  // Normal equations for (D, E, F), centred to keep the conditioning sane.
  const double a11 = sxx - sx * sx / dn;
  const double a12 = sxr - sx * sr / dn;
  const double a22 = srr - sr * sr / dn;
  const double b1 = 0.5 * (sxz - sx * sz / dn);
  const double b2 = 0.5 * (srz - sr * sz / dn);
  const double det = a11 * a22 - a12 * a12;
  const double scale = std::max(a11 * a22, a12 * a12);
  if (scale <= 0.0 || std::abs(det) < 1e-12 * scale) return false;
  const double cx = (b1 * a22 - b2 * a12) / det;
  const double cr = (b2 * a11 - b1 * a12) / det;
  double sum = 0.0;
  for (const auto& p : pts) sum += std::hypot(p.x - cx, p.r - cr);
  radius_out = sum / dn;
  return std::isfinite(radius_out) && radius_out > 0.0;
}

}  // namespace

NozzleGeometry NozzleGeometry::fromContour(std::vector<ContourSample> samples,
                                           NozzleGeometrySpec spec) {
  if (samples.size() < 9)
    throw ConfigError("nozzle: a tabulated contour needs at least 9 samples, got " +
                      std::to_string(samples.size()));
  std::sort(samples.begin(), samples.end(),
            [](const ContourSample& a, const ContourSample& b) { return a.x < b.x; });
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!(samples[i].r > 0.0)) {
      std::ostringstream os;
      os << "nozzle: contour sample " << i << " has a non-positive radius ("
         << samples[i].r << " m at x = " << samples[i].x << " m)";
      throw ConfigError(os.str());
    }
    if (i > 0 && !(samples[i].x > samples[i - 1].x)) {
      std::ostringstream os;
      os << "nozzle: contour samples must be strictly increasing in x, but sample "
         << i << " is at x = " << samples[i].x << " m and sample " << i - 1
         << " is at x = " << samples[i - 1].x << " m";
      throw ConfigError(os.str());
    }
  }

  NozzleGeometry g;
  g.tabulated_ = true;

  // --- throat: the minimum of the table -----------------------------------
  // Its resolution is the file's resolution.  A contour sampled coarsely
  // through the throat gives a throat area that is slightly too large, which
  // is a property of the input and is not silently corrected here.
  std::size_t it = 0;
  for (std::size_t i = 1; i < samples.size(); ++i)
    if (samples[i].r < samples[it].r) it = i;
  if (it < 3 || it + 3 >= samples.size())
    throw ConfigError(
        "nozzle: the throat of a tabulated contour must lie inside it with room on "
        "both sides -- the minimum radius is at sample " + std::to_string(it) + " of " +
        std::to_string(samples.size()) + ", which leaves no convergent or divergent "
        "section to resolve");

  g.rt_ = samples[it].r;
  g.xt_ = samples[it].x;
  g.at_ = constants::pi * g.rt_ * g.rt_;
  g.xe_ = samples.back().x;
  g.re_ = samples.back().r;
  g.rc_ = 0.0;
  for (std::size_t i = 0; i <= it; ++i) g.rc_ = std::max(g.rc_, samples[i].r);

  // --- throat curvature, measured ------------------------------------------
  // Bartz needs the wall radius of curvature at the throat.  Take the samples
  // whose radius is within 15 % of the throat -- an arc wide enough to fit but
  // still local to the throat -- and fit a circle to them.
  std::vector<ContourSample> near;
  for (const auto& s : samples)
    if (s.r <= 1.15 * g.rt_) near.push_back(s);
  double fitted = 0.0;
  if (fitCircle(near, fitted) && fitted > 0.05 * g.rt_ && fitted < 50.0 * g.rt_) {
    g.bartz_rc_ = fitted;
    std::ostringstream os;
    os << "circle fitted to " << near.size() << " samples within 15 % of the throat radius";
    g.curvature_note_ = os.str();
  } else {
    g.bartz_rc_ = spec.throat_upstream_ratio * g.rt_;
    std::ostringstream os;
    os << "circle fit failed or was implausible on " << near.size()
       << " near-throat samples, so throat_upstream_ratio ("
       << spec.throat_upstream_ratio << ") was used instead -- Bartz's (Dt/Rc)^0.1 term "
          "rests on an assumption here, not on the contour";
    g.curvature_note_ = os.str();
  }

  // --- back-fill the spec with what was measured ---------------------------
  // Downstream code consults spec(); it must see the real geometry.  The exit
  // wall angle matters most: NozzleFlow turns it into the divergence loss, and
  // defaulting it would silently apply an 8-degree bell's loss to whatever was
  // imported.
  spec.throat_radius = g.rt_;
  spec.throat_area = 0.0;
  spec.chamber_radius = g.rc_;
  spec.contraction_ratio = 0.0;
  spec.expansion_ratio = (g.re_ * g.re_) / (g.rt_ * g.rt_);
  const std::size_t last = samples.size() - 1;
  const double exit_slope = (samples[last].r - samples[last - 1].r) /
                            (samples[last].x - samples[last - 1].x);
  spec.bell_exit_angle = std::atan(exit_slope) / kDeg;
  spec.divergent = DivergentType::kBell;
  g.theta_e_ = spec.bell_exit_angle;
  const double init_slope = (samples[it + 2].r - samples[it + 1].r) /
                            (samples[it + 2].x - samples[it + 1].x);
  g.theta_n_ = std::atan(init_slope) / kDeg;
  spec.chamber_length = samples[0].x;
  g.spec_ = spec;

  g.l15_ = (g.rt_ * (std::sqrt(spec.expansion_ratio) - 1.0) +
            spec.throat_downstream_ratio * g.rt_ * (1.0 / std::cos(15.0 * kDeg) - 1.0)) /
           std::tan(15.0 * kDeg);

  // --- the wall is the polyline --------------------------------------------
  g.segments_.reserve(samples.size() - 1);
  for (std::size_t i = 1; i < samples.size(); ++i) {
    Segment s;
    s.kind = Segment::Kind::kLine;
    s.x0 = samples[i - 1].x;
    s.x1 = samples[i].x;
    s.r0 = samples[i - 1].r;
    s.slope = (samples[i].r - samples[i - 1].r) / (s.x1 - s.x0);
    // Label by where the segment sits and what the wall is doing, which is all
    // a table supports: it has no analytic pieces to name.
    if (s.x1 <= g.xt_)
      s.label = (std::abs(s.slope) < 1e-6) ? NozzleSegment::kChamber
                                           : NozzleSegment::kConvergingCone;
    else
      s.label = NozzleSegment::kDivergent;
    g.segments_.push_back(s);
  }
  // The two segments touching the throat are its curvature, however the file
  // resolved it.
  g.segments_[it - 1].label = NozzleSegment::kThroatUpstream;
  g.segments_[it].label = NozzleSegment::kThroatDownstream;

  // --- stations are the samples themselves ---------------------------------
  // Resampling a polyline onto a uniform grid would only interpolate what is
  // already there while moving the throat off a station.
  g.stations_.reserve(samples.size());
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const Segment& s = g.segments_[std::min(i, g.segments_.size() - 1)];
    ContourPoint pnt;
    pnt.x = samples[i].x;
    pnt.r = samples[i].r;
    pnt.area = constants::pi * pnt.r * pnt.r;
    pnt.area_ratio = pnt.area / g.at_;
    pnt.drdx = s.slope;
    pnt.segment = s.label;
    g.stations_.push_back(pnt);
  }

  g.finalise();
  return g;
}

void NozzleGeometry::finalise() {
  // --- validity ---------------------------------------------------------
  // The throat must be the unique minimum of the area distribution.
  double r_min = 1e300;
  std::size_t i_min = 0;
  for (std::size_t i = 0; i < stations_.size(); ++i)
    if (stations_[i].r < r_min) { r_min = stations_[i].r; i_min = i; }
  if (std::abs(r_min - rt_) > 1e-9 * rt_) {
    std::ostringstream os;
    os << "nozzle: the minimum sampled radius (" << r_min << " m at x = "
       << stations_[i_min].x << " m) is not the throat radius (" << rt_ << " m)";
    throw ConfigError(os.str());
  }
  for (std::size_t i = 1; i < stations_.size(); ++i) {
    const auto& a = stations_[i - 1];
    const auto& b = stations_[i];
    if (b.x <= a.x) throw ConfigError("nozzle: station positions are not strictly increasing");
    const double tol = -1e-9 * rt_;
    if (b.x <= xt_ && (b.r - a.r) > -tol && (b.r - a.r) > 1e-9 * rt_) {
      std::ostringstream os;
      os << "nozzle: the converging section is not monotonically contracting -- the wall "
            "rises from " << a.r << " m at x = " << a.x << " m to " << b.r << " m at x = "
         << b.x << " m, upstream of the throat at x = " << xt_ << " m";
      throw ConfigError(os.str());
    }
    if (a.x >= xt_ && (b.r - a.r) < tol) {
      std::ostringstream os;
      os << "nozzle: the diverging section is not monotonically expanding -- the wall "
            "falls from " << a.r << " m at x = " << a.x << " m to " << b.r << " m at x = "
         << b.x << " m, downstream of the throat at x = " << xt_ << " m";
      throw ConfigError(os.str());
    }
  }

  // --- chamber volume (solid of revolution up to the throat) ------------
  double vol = 0.0;
  for (std::size_t i = 1; i < stations_.size(); ++i) {
    const auto& a = stations_[i - 1];
    const auto& b = stations_[i];
    if (b.x > xt_ + 1e-15) break;
    const double dx = b.x - a.x;
    vol += constants::pi * dx * (a.r * a.r + a.r * b.r + b.r * b.r) / 3.0;   // exact for a frustum
  }
  v_chamber_ = vol;
}

const NozzleGeometry::Segment& NozzleGeometry::segmentFor(double x) const {
  const double tol = 1e-9 * std::max(xe_, 1e-9);
  for (const auto& s : segments_)
    if (x >= s.x0 - tol && x <= s.x1 + tol) return s;
  std::ostringstream os;
  os << "nozzle: axial position " << x << " m lies outside the contour [0, " << xe_ << "] m";
  throw RangeError(os.str());
}

double NozzleGeometry::radius(double x) const { return segmentFor(x).radius(x); }

NozzleSegment NozzleGeometry::segmentAt(double x) const { return segmentFor(x).label; }

double NozzleGeometry::wettedArea(double x_from, double x_to) const {
  if (x_to < x_from) std::swap(x_from, x_to);
  const int n = 2000;
  double area = 0.0;
  double xp = x_from, rp = radius(x_from);
  for (int i = 1; i <= n; ++i) {
    const double x = x_from + (x_to - x_from) * i / n;
    const double r = radius(x);
    area += constants::pi * (rp + r) * std::hypot(x - xp, r - rp);
    xp = x;
    rp = r;
  }
  return area;
}

void NozzleGeometry::checkContinuity(double& max_radius_jump, double& max_slope_jump) const {
  max_radius_jump = 0.0;
  max_slope_jump = 0.0;
  const double eps = 1e-7 * xe_;
  for (std::size_t i = 1; i < segments_.size(); ++i) {
    const double xj = segments_[i].x0;
    const double r_left = segments_[i - 1].radius(xj);
    const double r_right = segments_[i].radius(xj);
    max_radius_jump = std::max(max_radius_jump, std::abs(r_left - r_right) / rt_);
    const double s_left = segments_[i - 1].slopeAt(std::max(xj - eps, segments_[i - 1].x0));
    const double s_right = segments_[i].slopeAt(std::min(xj + eps, segments_[i].x1));
    max_slope_jump = std::max(max_slope_jump, std::abs(s_left - s_right));
  }
}

std::string NozzleGeometry::summary() const {
  std::ostringstream os;
  os << std::fixed << std::setprecision(5);
  os << "nozzle contour (" << toString(spec_.divergent) << ")\n"
     << "  throat radius        " << rt_ * 1e3 << " mm  (At = " << at_ * 1e4 << " cm^2)\n"
     << "  chamber radius       " << rc_ * 1e3 << " mm  (contraction " << std::setprecision(3)
     << contractionRatio() << ")\n" << std::setprecision(5)
     << "  exit radius          " << re_ * 1e3 << " mm  (expansion " << std::setprecision(3)
     << expansionRatio() << ")\n" << std::setprecision(5)
     << "  throat at x          " << xt_ * 1e3 << " mm\n"
     << "  exit at x            " << xe_ * 1e3 << " mm\n"
     << "  divergent length     " << divergentLength() * 1e3 << " mm ("
     << std::setprecision(1) << 100.0 * divergentLength() / l15_ << " % of 15 deg conical)\n"
     << std::setprecision(5)
     << "  chamber volume       " << v_chamber_ * 1e6 << " cm^3\n"
     << "  characteristic L*    " << characteristicLength() * 1e3 << " mm\n"
     << "  Bartz curvature      " << bartz_rc_ * 1e3 << " mm (" << spec_.bartz_curvature << ")\n"
     << "  divergence angles    theta_n = " << std::setprecision(2) << theta_n_ / kDeg
     << " deg, theta_e = " << theta_e_ / kDeg << " deg\n"
     << "  stations             " << stations_.size();
  return os.str();
}

}  // namespace ignis
